#include "brain.h"

/* Model-specific methods:
 *     p0 defines the time-dependent pressure input
 *     rhs defines the differential equations
 */
double p0(double t) {
    double p0;
    p0 = 1.5;
    return p0;
}
void rhs(workspace *W, double t, double *y, double *p, double *dy) {
    double r, f, cb, ct;
    double pt, e, r0, q, g;
    int ir, ismc, icb, ict;
    for (int i = 0; i < W->P->n_auto; i++) {
        ir   = 4*i + W->P->i_radius;
        ismc = 4*i + W->P->i_smc;
        ict  = 4*i + W->P->i_ctissue;
        icb  = 4*i + W->P->i_cblood;

        r  = y[ir];
        g  = pow(r, 4) / W->P->l[i];
        f  = y[ismc];
        ct = y[ict];
        cb = y[icb];

        pt = 0.5 * (p[i/2] + W->P->pcap);
        q  = (p[i/2] - W->P->pcap) * g;
        e  = 1. + W->P->a5[i] * f;
        r0 = W->P->a3[i] * (1 - W->P->a4[i] * f);

        dy[ir]   = -W->P->a1[i] * e * (r / r0 - 1) + W->P->a2[i] * r * pt;
        dy[ismc] = -W->P->b1[i] * \
                   (f - 1 / (1 + exp(W->P->gamma[i] * (ct - W->P->cstar[i]))));
        dy[icb]  =  W->P->d1[i] * q * (1 - cb) + W->P->d2[i] * (ct - cb);
        dy[ict]  = -W->P->g1[i] * (ct - cb) + W->P->g2[i];
    }

}

/* Methods for external use: essentially function evaluation and Jacobian
 * updates. evaluate, solve and jacupdate both have the side effect of updating
 * conductance, pressure, and flow values in the workspace data structure,
 * so these values should not be relied to remain static
 * */
workspace * init(int argc, char **argv) {
    workspace *W;
    W = malloc(sizeof *W);

    init_geometry(W, argc, argv);
    init_subtree(W);
    init_roottree(W);
    compute_symbchol(W);
    init_problem(W);

    return W;
} 
void evaluate(workspace *W, double t, double *y, double *dy) {
    double r, l;

    // Update conductances of autoregulating vessels
    for (int i = 0; i < W->P->n_auto; i++) {
        r = y[4*i + W->P->i_radius];
        l = W->P->l[i];
        W->g[i] = pow(r, 4) / l;
    }
    // printf("%4d: Boo! ", W->rank);
    // Solve for pressure and flow
    solve(W, p0(t), W->P->pcap);
    // Evaluate the right hand side equations
    rhs(W, t, y, W->p, dy);
}
void solve(workspace *W, double pin, double pc) {
    compute_uv(W, pc);
    if (W->n_procs > 1) {
        communicate(W);
        compute_root(W, pin);
    }
    compute_sub(W, pin, pc);
}
void jacupdate(workspace *W, double t, double *y) {
    double *f;
    double eps = 1e-6;

    // Evaluate the right hand side
    f = malloc(W->P->ny * sizeof (*f ));
    evaluate(W, t, y, f);
    //rhs(W, t, y, W->p, f);

    // This order is important. dpdg depends on dgdx, and on the correct
    // value of W->w from the call to evaluate. dfdx and dfdp will modify
    // this value
    eval_dgdx(W, t, y);
    eval_dpdg(W, t, y);
    eval_dfdx(W, t, y, f, eps);
    eval_dfdp(W, t, y, f, eps); 

    if (W->isjac) cs_spfree(W->J);
    cs *Y;
    // Compute the product dfdx + dfdp dpdg dgdx

    Y = cs_multiply(W->dfdp->A, W->dpdgneg);
    W->J = cs_add(W->dfdx->A, Y, 1, -1);
    W->isjac = 1;
    cs_spfree(Y);


    free(f);
}

/* Private functions */
void init_geometry(workspace *W, int argc, char **argv) {
    // Sort out MPI configuration
    MPI_Comm_size(MPI_COMM_WORLD, &W->n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &W->rank);
    assert(is_power_of_two(W->n_procs));

    // Parse input parameters and set tree sizes
    W->N    = NDEFAULT; 
    W->Nsub = NSUBDEFAULT;
    if (argc > 1)
        W->N = atoi(argv[1]);
    if (argc > 2)
        W->Nsub = atoi(argv[2]);
    W->N0 = (int) round(log2((double) W->n_procs));
    W->Np = W->N - W->N0;

    // Configure buffer and parameters for MPI_Allgather)
    W->buf  = malloc(W->n_procs * NSYMBOLS * sizeof(*W->buf));
    W->flag = malloc(W->n_procs * sizeof (*W->flag));

}
int is_power_of_two (unsigned int x) {
      return ((x != 0) && !(x & (x - 1)));
}
void init_subtree(workspace *W) {
    // First construct the local subtree 
    W->A    = adjacency(W->Np);
    W->At   = cs_transpose(W->A, 1);
    W->G    = speye(W->A->n);
    W->g    = W->G->x;
    W->level = malloc(W->A->n * sizeof (*W->level));
    for (int i = 0; i < W->A->n; i++) 
        W->level[i] = (int) floor(log2(W->A->n - i)) + W->N0;

    // Initialise workspace variables for solving
    W->b = malloc (W->A->n * (sizeof *W->b));
    W->u = malloc (W->A->m * (sizeof *W->u));
    W->v = malloc (W->A->m * (sizeof *W->v));
    W->p = malloc (W->A->m * (sizeof *W->p));
    W->q = malloc (W->A->n * (sizeof *W->q));
    W->w = malloc (W->A->n * (sizeof *W->w));
    W->ucomm = malloc (W->n_procs * (sizeof *W->ucomm));
    W->vcomm = malloc (W->n_procs * (sizeof *W->vcomm));
    W->gcomm = malloc (W->n_procs * (sizeof *W->gcomm));
    // Initialise general purpose vectors for workspaces
    W->xn = malloc(W->A->n * sizeof(*W->xn));
    W->xm = malloc(W->A->m * sizeof(*W->xm));
}
cs *adjacency(int N) {
    int m = (1 << (N-1)) - 1;
    int n = (1 << N) - 1;
    int *Ti, *Tj;
    double *Tx;

    cs *T, *A;
    T = cs_spalloc(m, n, 3*m, 1, 1);
    if (T == NULL) return T;
    Ti = T->i; Tj = T->p; Tx = T->x;

    int k = 0;
    for (int i = 0; i < m; i++) {
        Ti[k] = i;   Ti[k+1] = i;       Ti[k+2] = i;
        Tj[k] = 2*i; Tj[k+1] = 2*i + 1; Tj[k+2] = m + i + 1;
        Tx[k] = 1.;  Tx[k+1] = 1.;      Tx[k+2] = -1.;
        k += 3;
    }
    T->nz = 3*m;
    A = cs_compress(T);
    cs_spfree(T);
    return A;
}
void compute_symbchol(workspace *W) {
    cs *X, *Y;
    X = cs_multiply(W->A, W->G);
    Y = cs_multiply(X, W->At);
    cs_spfree(X);

    W->symbchol = cs_schol(0, Y);
    cs_spfree(Y);

    X = cs_multiply(W->A0, W->G0);
    Y = cs_multiply(X, W->A0t);
    cs_spfree(X);

    W->symbchol0 = cs_schol(0, Y);
    cs_spfree(Y);
}
void init_roottree(workspace *W) {
    cs *T;
    int *ikeep;
    T = adjacency(W->N0 + 1);
    ikeep = malloc(T->m * (sizeof *ikeep));
    for (int i = 0; i < T->m; i++)
        ikeep[i] = T->n - T->m + i;
    W->A0 = subsref(T, NULL, ikeep, -1, T->m);
    cs_spfree(T);
    W->A0t = cs_transpose(W->A0, 1);
    W->G0 = speye(W->A0->n);
    W->g0 = W->G0->x;
    W->level0 = malloc(W->A->n * sizeof (*W->level0));
    for (int i = 0; i < W->A0->n; i++)
        W->level0[i] = (int) floor(log2(W->A0->n - i));
    // Initialise workspace variables for solving 
    W->b0 = malloc (W->A0->n * sizeof (*W->b0));
    W->p0 = malloc (W->A0->m * sizeof (*W->b0));
    W->q0 = malloc (W->A0->n * sizeof (*W->b0));
    W->xn0 = malloc (W->A0->n * sizeof (*W->xn0));
}
void init_problem(workspace *W) {
    /* The trees, but not conductances, have been initialised. We need to
     * set up all of the ODE parameters, and set the conductances in the
     * tree */

    // Initialise conductances to their scaled values 
    set_conductance(W, 1, 1); // set unscaled conductances
    W->P = init_params(W);
    set_conductance(W, 0, 1); // set scaled conductances
    init_jacobians(W);        // Initialise Jacobian data structures
}
void set_conductance(workspace *W, int isunscaled, int computeroot) {
    double r, l;
    if (isunscaled) {
        for (int i = 0; i < W->A->n; i++) {
            r = compute_radius(W->level[i], W->N);
            l = compute_length(W->level[i], W->N);
            W->g[i] = M_PI * pow(r, 4) / (8. * MU * l);
        }
    } else {
        for (int i = 0; i < W->A->n; i++) {
            r = compute_radius(W->level[i], W->N) / R0;
            l = compute_length(W->level[i], W->N) / L0;
            W->g[i] = pow(r, 4) / l;
        }
    } if (computeroot) {
        if (isunscaled) {
            for (int i = 0; i < W->A0->n; i++) {
                r = compute_radius(W->level0[i], W->N);
                l = compute_length(W->level0[i], W->N);
                W->g0[i] = M_PI * pow(r, 4) / (8. * MU * l);
            }
        } else {
            for (int i = 0; i < W->A0->n; i++) {
                r = compute_radius(W->level0[i], W->N) / R0;
                l = compute_length(W->level0[i], W->N) / L0;
                W->g0[i] = pow(r, 4) / l;
            }
        }
    }
}
double compute_length(int level, int n_levels) {
    double l;
    l = LRR * RMIN * (double) (1 << ((n_levels - level - 1)/ 2));
    return l;
}
double compute_radius(int level, int n_levels) {
    double r;
    r = RMIN * pow(2., ((double) (n_levels - level - 1)) / 2.);
    return (r);
}
user *init_params(workspace *W) {
    int na, i_start = 0;
    double Rstar, hstar;
    double Lmin, Acap;
    double Vb, Vt;
    double Q0;
    user *P;

    P = malloc(sizeof *P);

    // Define the indices for the state variables. At the moment we assume
    // only the tail end vertices are autoregulating.
    // TODO: Allow more generic set of autoregulating vessels 
    P->n_auto    = W->A->n - W->A->m;
    P->ny = 4*P->n_auto;
    na = P->n_auto;
    P->i_radius  = 0;
    P->i_cblood  = 1;
    P->i_smc     = 2;
    P->i_ctissue = 3;

    // Define the l vector for autoregulating vessels
    P->l = malloc(na * sizeof(*P->l));
    for (int i = 0; i < na; i++)
        P->l[i] = compute_length(W->level[i_start + i], W->N) / L0; 

    // Define a workspace for the transmural pressure 
    // Define a sparse matrix that can compute it (essentially a projection
    // of abs(A) 
    
    // Radial parameters 
    P->a1 = zerosv(na);
    P->a2 = zerosv(na);
    P->a3 = zerosv(na);
    P->a4 = zerosv(na);
    P->a5 = zerosv(na);
    for (int i = 0; i < na; i++) {
        Rstar = compute_radius(W->level[i], W->N);
        hstar = HRR * Rstar;
        P->a1[i] = E0 * T0 * Rstar / (ETA * R0);
        P->a2[i] = P0 * Rstar * T0 / (ETA * hstar);
        P->a3[i] = Rstar / R0;
        P->a4[i] = 1 - RSCALE;
        P->a5[i] = EACTIVE / EPASSIVE - 1;
    }
    Lmin = LRR * RMIN;
    Acap = pow(2*Lmin, 3) / (RCAP / 2. + 1. / AVR);
    Vt   = Acap / AVR;
    Vb   = RCAP / 2 * Acap;
    Q0   = M_PI * pow(R0, 4) * P0 / (8 * MU * L0);
    P->d1 = zerosv(na);
    P->d2 = zerosv(na);
    P->g1 = zerosv(na);
    P->g2 = zerosv(na);
    P->b1 = zerosv(na);
    P->gamma = zerosv(na);
    P->cstar = zerosv(na);

    for (int i = 0; i < na; i++) {
        P->d1[i] = T0 * Q0 / Vb;
        P->d2[i] = T0 * PW * Acap / Vb;
        P->g1[i] = T0 * PW * Acap / Vt;
        P->g2[i] = T0 * M0 / C0;
        P->b1[i] = T0 / TAUM;
        P->gamma[i] = GAMMA;
        P->cstar[i] = CSTAR;
    }
    P->pcap = PCAP / P0;
    return P;
}
void init_jacobians(workspace *W) {
    // Create the data structures for each of the Jacobians
    W->isjac = 0;
    init_dgdx(W);
    init_dpdg(W);
    init_dfdx(W);
    init_dfdp(W);
}
void init_dgdx(workspace *W) {
    cs *T;
    int *Ti, *Tj;
    double *Tx;
    int neq = 4;
    int na, ir;
    na = W->P->n_auto;
    ir = W->P->i_radius;

    T = cs_spalloc(W->A->n, W->P->ny, na, 1, 1);
    Ti = T->i; Tj = T->p; Tx = T->x;
    for (int i = 0; i < na; i++) {
        Ti[i] = i;
        Tj[i] = neq*i + ir;
        Tx[i] = 1.;
    }
    T->nz = na;
    W->dgdx = cs_compress(T);
    cs_spfree(T);
}
void init_dpdg(workspace *W) {
    // In this function we'll figure out A_A, etc, and
    //     * compute the symbolic factorisation of A_A G A_A^T

    // Compute the indices B of the nodes to remove from A
    int Nr = W->Np - W->Nsub; int m = W->A->m;
    int B0 = (1 << (W->Np - 1)) - (1 << Nr);
    int B1 = (1 << (W->Np - 1)) - (1 << (Nr - 1));

    // Construct the temporary matrix to store the projector in
    cs *T;
    int *Ti, *Tj; double *Tx;
    T = cs_spalloc(m - (B1 - B0), m, m - (B1 - B0), 1, 1);
    Ti  = T->i; Tj = T->p; Tx = T->x;

    int k = 0;
    for (int i = 0; i < B0; i++) {
        Ti[k] = k; Tj[k] = i; Tx[k++] = 1.;
    }
    for (int i = B1; i < m; i++) {
        Ti[k] = k; Tj[k] = i; Tx[k++] = 1.;
    }
    T->nz = k;
    cs *Pt;
    Pt = cs_compress(T);
    cs_spfree(T);
    
    // Compute A_A, A_A.', and the matrix Proj used to construct dp/dg
    W->A_A  = cs_multiply(Pt, W->A);
    W->A_At = cs_transpose(W->A_A, 1);
    W->Proj = cs_transpose(Pt, 1);
    cs_spfree(Pt);

    // Compute symbolic Cholesky factorisation of A_A G A_A^T
    T = cs_multiply(W->A_A, W->A_At);
    W->symbchol_reduced = cs_schol(0, T);
    cs_spfree(T);
}
void init_dfdx(workspace *W) {
    int j = 0, k = 0, na, ir, ismc, ict, icb;
    int nzloc = 9;
    int neqloc = 4;
    int *Ti, *Tj;
    double *Tx;
    na = W->P->n_auto;
    cs *T, *J;
    ir = W->P->i_radius;
    ismc = W->P->i_smc;
    ict = W->P->i_ctissue;
    icb = W->P->i_cblood;

    T = cs_spalloc(neqloc*na, neqloc*na, nzloc*na, 1, 1);
    if (T == NULL) fprintf(stderr, "failed to allocate matrix in create_dfdx");
    Ti = T->i; Tj = T->p; Tx = T->x;
    for (int i = 0; i < na; i++) {
        j = neqloc * i;
        k = nzloc * i;
        Ti[k    ] = j + ir;   Tj[k    ] = j + ir;   Tx[k    ] = 1.;
        Ti[k + 1] = j + ir;   Tj[k + 1] = j + ismc; Tx[k + 1] = 1.;
        Ti[k + 2] = j + ismc; Tj[k + 2] = j + ismc; Tx[k + 2] = 1.;
        Ti[k + 3] = j + ismc; Tj[k + 3] = j + ict;  Tx[k + 3] = 1.;
        Ti[k + 4] = j + ict;  Tj[k + 4] = j + ict;  Tx[k + 4] = 1.;
        Ti[k + 5] = j + ict;  Tj[k + 5] = j + icb;  Tx[k + 5] = 1.;
        Ti[k + 6] = j + icb;  Tj[k + 6] = j + ict;  Tx[k + 6] = 1.;
        Ti[k + 7] = j + icb;  Tj[k + 7] = j + icb;  Tx[k + 7] = 1.;
        Ti[k + 8] = j + icb;  Tj[k + 8] = j + ir;   Tx[k + 8] = 1.;
    }
    T->nz = nzloc * na;
    J = cs_compress(T);
    cs_spfree(T);
    W->dfdx = numjacinit(J);
    cs_spfree(J);
}
void init_dfdp(workspace *W) {
    int j = 0, k = 0, na, ir, icb;
    int neqloc = 4, nzloc = 2;
    int *Ti, *Tj;
    double *Tx;
    cs *T, *J;

    na = W->P->n_auto;
    ir  = W->P->i_radius;
    icb = W->P->i_cblood;

    T = cs_spalloc(neqloc*na, W->A->m, nzloc*na, 1, 1);
    if (T == NULL) fprintf(stderr, "failed to allocate matrix in create_dfdx");
    Ti = T->i; Tj = T->p; Tx = T->x;
    for (int i = 0; i < na; i++) {
        k = nzloc * i;
        j = neqloc * i;
        Ti[k    ] = j + ir;    Tj[k    ] = i/2; Tx[k    ] = 1.;
        Ti[k + 1] = j + icb;   Tj[k + 1] = i/2; Tx[k + 1] = 1.;
    }
    T->nz = nzloc * na;
    J = cs_compress(T);

    cs_spfree(T);
    W->dfdp = numjacinit(J);
    cs_spfree(J);

}
void compute_uv(workspace *W, double pc) {
    cs *AG, *B;
    csn *Nu;

    // Set up boundary conditions - put them in q for now 
    for (int i = 0; i < W->A->m + 1; i++)
        W->q[i] = -pc;
    for (int i = W->A->m + 1; i < W->A->n; i++)
        W->q[i] = 0;

    // Define the matrices AG = A*G, and B = A*G*A'
    AG = cs_multiply(W->A, W->G);
    B  = cs_multiply(AG, W->At);

    // Numerical Cholesky factorisation using precomputed symbolic
    Nu = cs_chol(B, W->symbchol);
    if (!Nu) printf("Numerical cholesky decomposition failed in compute_uv");
    cs_spfree(B); // B is (potentially) big, so free it

    // Define the RHS of the u equations
    for (int i = 0; i < W->A->m; i++)
        W->u[i] = 0.;                       /* u = 0 */
    cs_gaxpy(AG, W->q, W->u);               /* u = AGb + u = AGb */
    for (int i = 0; i < W->A->m; i++)
        W->u[i] = -W->u[i];                 /* u = -u = -AGb */

    cs_spfree(AG); 

    // And solve, using our computed factorisations
    cholsoln(Nu, W->symbchol, W->A->m, W->u, W->xm);

    // Define RHS of v equations
    W->v[W->A->m-1] = W->g[W->A->n-1];
    for (int i = 0; i < W->A->m - 1; i++)
        W->v[i] = 0.;

    // And solve, using our computed factorisations
    cholsoln(Nu, W->symbchol, W->A->m, W->v, W->xm);
    cs_nfree(Nu); 

    // Put the end entries into communication buffers
    W->ucomm[W->rank] = W->u[W->A->m-1];
    W->vcomm[W->rank] = W->v[W->A->m-1];
    W->gcomm[W->rank] = W->g[W->A->n-1];
}
void communicate(workspace *W) {
    double send_buf[NSYMBOLS], *recv_buf;

    recv_buf = W->buf;

    // Put local variables into the communication buffer
    send_buf[0] = W->ucomm[W->rank];
    send_buf[1] = W->vcomm[W->rank];
    send_buf[2] = W->gcomm[W->rank];
    send_buf[3] = (double) W->flag[W->rank];

    // Fill up the buffer with MPI all-to-all communication
    MPI_Allgather(send_buf, NSYMBOLS, MPI_DOUBLE, recv_buf, NSYMBOLS, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Populate u1, v1, g1
    for (int i = 0; i < W->n_procs; i++) {
        if (i != W->rank) {
            W->ucomm[i] = recv_buf[NSYMBOLS * i    ];
            W->vcomm[i] = recv_buf[NSYMBOLS * i + 1];
            W->gcomm[i] = recv_buf[NSYMBOLS * i + 2];
            W->flag[i]  = (int) recv_buf[NSYMBOLS * i + 3]; 
        }
    }
}
void compute_root(workspace *W, double pin) {
    cs *AG, *B, *X, *D;
    csn *Nu;
    double *d;

    d = W->xn0;
    for (int i = 0; i < W->A0->n; i++)
        d[i] = 0.;

    // Define the diagonal modification D to be added to AGA'
    for (int i = 0; i < W->n_procs; i++)
        d[i/2] += W->gcomm[i] * (1 - W->vcomm[i]); // Parent is simply i/2 

    D = spdiags(d, W->A0->n);

    // Define the matrices AG = A_0 G_0 and B = A_0 G_0 A_0^T + D 
    AG = cs_multiply(W->A0, W->G0);
    X  = cs_multiply(AG, W->A0t);
    B = cs_add(X, D, 1., 1.);
    cs_spfree(D); 
    cs_spfree(X);

    // Perform numerical Cholesky factorisation 
    Nu = cs_chol(B, W->symbchol0);
    if (!Nu) printf("Numerical Cholesky failed in compute_root\n");
    cs_spfree(B);

    // Define the boundary condition b (stick it in q0 for now)
    W->q0[W->A0->n-1] = pin;
    for (int i = 0; i < (W->A0->n - 1); i++)
        W->q0[i] = 0.;

    // Compute the right hand side of the p equation
    for (int i = 0; i < W->A0->m; i++) 
        W->p0[i] = 0.;                              // p0 = 0
    cs_gaxpy(AG, W->q0, W->p0);                     // p0 = AGb + p0 = AGb 
    for (int i = 0; i < W->A0->m; i++)
        W->p0[i] = -W->p0[i];                       // p = -p = -AGb */
    // p = -AGb + sum g_i u_i e_ki 
    for (int i = 0; i < W->n_procs; i++)
        W->p0[i/2] += W->gcomm[i] * W->ucomm[i];   
    cs_spfree(AG);

    // And solve, using our numerical Cholesky
    cholsoln(Nu, W->symbchol0, W->A0->m, W->p0, W->xn0);          
    // W->p0 is overwritten with p 

    // Now compute q0: W->q0 is currently p_in e_n 
    cs_gaxpy(W->A0t, W->p0, W->q0);       // q = A'p + q = A'p + pin en 
    for (int i = 0; i < W->A0->m; i++)
        W->q0[i] *= W->g0[i];
    cs_nfree(Nu);
}
void compute_sub(workspace *W, double pin, double pc) {
    double pk;
    if (W->n_procs == 1)
        pk = pin;
    else
        pk = W->p0[W->rank / 2];
    for (int i = 0; i < W->A->m; i++)
        W->p[i] = W->u[i] + pk * W->v[i];

    // figure out w - first set it to b + p0k e_n (effectively b) 
    W->w[W->A->n - 1] = pk;
    for (int i = 0; i < W->A->m + 1 ; i++)
        W->w[i] = -pc;
    for (int i = W->A->m + 1; i < (W->A->n - 1); i++)
        W->w[i] = 0;

    cs_gaxpy(W->At, W->p, W->w); // compute A'p + b + pKi en 
    for (int i = 0; i < W->A->n; i++)
        W->q[i] = W->w[i] * W->g[i]; // q = w g
}
void eval_dgdx(workspace *W, double t, double *y) {
    int neqloc = 4;
    double r, l;
    for (int i = 0; i < W->P->n_auto; i++) {
        r = y[i * neqloc + W->P->i_radius];
        l = W->P->l[i];
        W->dgdx->x[i] = 4.*pow(r, 3) / l;
    }
}
void eval_dpdg(workspace *W, double t, double *y) {
    // Free the old one
    if (W->isjac) cs_spfree(W->dpdgneg);

    // Update the conductance matrix
    double r;
    for (int i = 0; i < W->P->n_auto; i++) {
        r = y[4*i + W->P->i_radius];
        W->g[i] = pow(r, 4) / W->P->l[i];
    }

    // Form the matrix A_A G A_A.'
    cs *X, *Y;
    Y = cs_multiply(W->A_A, W->G);
    X = cs_multiply(Y, W->A_At);

    cs_spfree(Y);

    // Form the right hand side matrix
    cs *D, *B, *C;
    D = spdiags(W->w, W->A->n); // w is set by the last evaluate
    C = cs_multiply(D, W->dgdx);
    cs_spfree(D);
    B = cs_multiply(W->A_A, C);
    cs_spfree(C);

    // Solve to form dpdg_A
    cs *PA;
    PA = mldivide_chol(X, W->symbchol_reduced, B);
    W->dpdgneg = cs_multiply(W->Proj, PA);
    cs_spfree(X);
    cs_spfree(PA);
    cs_spfree(B);
}
void eval_dfdx(workspace *W, double t, double *y, double *f, double eps) {
    int i, j;
    double *y1, *h, *f1;

    y1 = malloc(W->P->ny * sizeof (*y1));
    h  = malloc(W->P->ny * sizeof (*h));
    f1 = malloc(W->P->ny * sizeof (*f1));
    
    for (int igrp = 0; igrp < W->dfdx->ng; igrp++) {
        for (int k = 0; k < W->P->ny; k++)
            y1[k] = y[k];

        for (int k = W->dfdx->r[igrp]; k < W->dfdx->r[igrp + 1]; k++) {
            j = W->dfdx->g[k];
            h[j] = eps; // * fabs(y[j]); 
            y1[j] += h[j];
        }
        rhs(W, t, y1, W->p, f1);  // f(t, y + dy, p)
        for (int k = W->dfdx->r[igrp]; k < W->dfdx->r[igrp+1]; k++) {
            j = W->dfdx->g[k];
            for (int ip = W->dfdx->A->p[j]; ip < W->dfdx->A->p[j+1]; ip++) {
                i = W->dfdx->A->i[ip];
                W->dfdx->A->x[ip] = (f1[i] - f[i]) / h[j];
            }
        }
    }
    free(y1); free(h); free(f1);
}
void eval_dfdp(workspace *W, double t, double *y, double *f, double eps) {
    int i, j;
    double *h, *p1, *f1;
    h =  malloc(W->A->m  * sizeof (*h));
    p1 = malloc(W->A->m  * sizeof (*p1));
    f1 = malloc(W->P->ny * sizeof (*f1));

    for (int igrp = 0; igrp < W->dfdp->ng; igrp++) {
        // Set p1 back to p
        for (int k = 0; k < W->A->m; k++)
            p1[k] = W->p[k];

        // Increment entry of h for every entry column in the group 
        for (int k = W->dfdp->r[igrp]; k < W->dfdp->r[igrp+1]; k++) {
            j = W->dfdp->g[k];
            h[j] = eps; // * fabs(W->p[j]);
            p1[j] += h[j];
        }

        // Evaluate the right hand side
        rhs(W, t, y, p1, f1);
        // Iterate over the columns in the gropu
        for (int k = W->dfdp->r[igrp]; k < W->dfdp->r[igrp+1]; k++) {
            j = W->dfdp->g[k];
            // and then the rows in the column
            for (int ip = W->dfdp->A->p[j]; ip < W->dfdp->A->p[j+1]; ip++) {
                i = W->dfdp->A->i[ip];
                W->dfdp->A->x[ip] = (f1[i] - f[i]) / h[j];
            }
        }
    }
    free(h); free(f1); free(p1);
}
