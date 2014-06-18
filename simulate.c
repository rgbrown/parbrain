#include "brain.h"

typedef struct odews {
    workspace *W;
    csn *N; // Newton matrix numeric factorisation
    css *S; // Newton matrix sybolic factorisation
    double *y; // Workspace variable
    double *p; // - hab ich dazugefuegt, funktioniert aber anscheinend so nicht...
    double *q; // - hab ich dazugefuegt, funktioniert aber anscheinend so nicht...
    double *f; // Workspace variable
    double gamma;
    double t0;
    double tf;
    double ftol;
    double ytol;
    int    maxits;
    int    nconv;
    int mdeclared;
    double dtwrite;
} odews;

// prototypes
void newton_matrix(odews *ws);
int lusoln(odews *ws, double *b);
css * newton_sparsity(cs *J);
int sizecheck(double *x, int n, double tol);
void back_euler(odews *ws);
void solver_init(odews *ws, int argc, char **argv);

// Main simulation program 
int main(int argc, char **argv) {
    // Initialisation step. Construct and initialise the ODE workspace
    odews *ws;
    MPI_Init(&argc, &argv);
    ws = malloc(sizeof *ws);
    int verbose = 1;

    // Problem parameters
    ws->gamma  = 1e-5; // time step  1e-5
    ws->t0     = 0.;   // initial time 0
    ws->tf     = 0.1;  // final time  10
    ws->ftol   = 1e-3; // function evaluation tolerance for Newton convergence 1e-3
    ws->ytol   = 1e-3; // relative error tolerance for Newton convergence 1e-3
    ws->nconv  = 5;    // Newton iteration threshold for Jacobian reevaluation 5
    ws->maxits = 100;   // Maximum number of Newton iterations 100
    ws->dtwrite = 5e-2; // Time step for writing to file (and screen)

    // Initialise the solver with all the bits and pieces
    solver_init(ws, argc, argv);

    // Begin computation. t0 and tf just measure elapsed simulation time
    double t0 = MPI_Wtime();
    back_euler(ws); // All of the simulation occurs in here
    double tf = MPI_Wtime();

    // Display diagnostics
    if (ws->W->rank == 0) {
        if (verbose) {
            printf("Levels: %d Subtree size: %d N procs: %d\n", ws->W->N, ws->W->Nsub, ws->W->n_procs);
            printf("Solution time:                %g seconds\n", tf - t0);
            printf("    # fevals:                 %d\n", ws->W->fevals);
            printf("    # Jacobians:              %d\n", ws->W->jacupdates);
            printf("    # feval time:             %g seconds\n", ws->W->tfeval);
            printf("    # Jacobian update time:   %g seconds\n", ws->W->tjacupdate);
            printf("    # Jacobian symbolic time: %g seconds\n", ws->W->tjacfactorize);
        } else {
            printf("%4d%4d%4d%12.4e%4d%4d%12.4e%12.4e%12.4e\n", ws->W->N, ws->W->Nsub, ws->W->n_procs, tf - t0, ws->W->fevals, ws->W->jacupdates, ws->W->tfeval, ws->W->tjacupdate, ws->W->tjacfactorize);
        }
    }
    // And clean up
    close_io(ws->W);
    MPI_Finalize();
    return 0;
}

// Fixed step Backward Euler ODE solver
void back_euler(odews *ws) {
    // Declare and initialise additional workspace variables
    double *beta, *w, *x;
    workspace *W;
    W = ws->W;
    int ny = W->nu;
    beta = zerosv(ny);
    w    = zerosv(ny);
    x    = zerosv(ny);

    double t = ws->t0;
    double tnext;

    // Initial newton_matrix computation
    newton_matrix(ws);
    int jac_needed = 0; // flag to say that Jacobian is current

    int converged = 0;
    write_data(W, t, ws->y); // Write initial data to file
    // write_vtk(W, t, ws->y, W->p, W->q);
    for (int i = 0; t < ws->tf; i++) {
        // Perform a Jacobian update if necessary. This is in sync across
        // all processors
        if (jac_needed) {
            jacupdate(W, t, ws->y); 
            jac_needed = 0;
            newton_matrix(ws);
        }

        // Copy values from previous completed timestep
        dcopy(ny, ws->y, beta);
        dcopy(ny, ws->y, w);
        tnext = t + ws->gamma;

        // Indicate that we haven't converged yet
        W->flag[W->rank] = 0;
        converged = 0;
        for (int k = 0; k < ws->maxits; k++) {  //
            evaluate(W, tnext, w, ws->f); // f = g(w)

            // evaluate also exchanges convergence information. If everyone
            // has converged, then we can stop. Everyone likewise gets the
            // request for Jacobian update
            if (all(W->flag, W->n_procs)) {
                converged = 1;
                if (k > ws->nconv)
                    jac_needed = 1;
                break; // w contains the correct value 
            }
            // Form x = w - beta - gamma g (our fcn value for Newton)
            dcopy(ny, w, x);
            daxpy(ny, -1, beta, x);
            daxpy(ny, -ws->gamma, ws->f, x);
            for (int la = 0; la < 24; la++) {
		//printf("iteration %d, state variable %2d - x: %e w: %e\n", k, la, x[la], w[la] );
	    } // TEST x[0] = radius etc. - w is the state var value, x is passed on to Newton
            W->flag[W->rank] = sizecheck(x, ny, ws->ftol); // function value size check
            lusoln(ws, x);  // solve (x is now increment)
            W->flag[W->rank] |= sizecheck(x, ny, ws->ytol); // increment size check
            daxpy(ny, -1, x, w); // update w with new value
        } // Newton loop
        if (!converged) {
            printf("Newton iteration failed to converge\n");
            exit(1);
        }
        t = tnext;
        dcopy(ny, w, ws->y); // update y values
        //if (W->rank == 0) {   // TEST
		//if (t <= 0.1) {
        if (fmod(t, ws->dtwrite) < ws->gamma) {
            write_data(W, t, ws->y); //ws->p, ws->q);
		    //write_vtk(W, t, ws->y, W->p, W->q);
            if (W->rank == 0) 
                printf("time: %e \n",t);
        }
		//}
        //}
    } // timestep loop
}
void solver_init(odews *ws, int argc, char **argv) {
    workspace *W;

    // Initialise the workspace
    ws->W = init(argc, argv); W = ws->W;
    ws->mdeclared = 0;

    // Put initial conditions in to y
    ws->y = zerosv(W->nu);
    ws->f = zerosv(W->nu);
    set_initial_conditions(W, ws->y);

    // Initial Jacobian computation
    double t0 = MPI_Wtime();
    evaluate(W, ws->t0, ws->y, ws->f);
    double tf = MPI_Wtime();
    ws->W->tfeval = tf - t0;
    t0 = MPI_Wtime();
    jacupdate(W, ws->t0, ws->y);
    double ta = MPI_Wtime();
    ws->S = newton_sparsity(W->J);
    double tb = MPI_Wtime();
    newton_matrix(ws);
    tf = MPI_Wtime();
    ws->W->tjacupdate = (tf - t0) - (tb - ta);
    ws->W->tjacfactorize = (tb - ta);
}
int sizecheck(double *x, int n, double tol) { // n - no of equ total (nblocks*nequs)
    int smallenough = 1;
    double x0[34] = 	{1,   	// 0
		 	1e-7,	// 1
			1e-4,	// 2
			1e-3,	// 3
			1e-4,	// 4
			1e-4,	// 5
			1e-3,	// 6
			1e-5, 	// 7
			1e-4,	// 8	
			1e+3,	// 9
			1e-4,	// 10
			1e-1,	// 11 *
			1.,	// 12
			1e-1,	// 13
			1e-1, 	// 14
			1e-1,	// 15 **
			1e+5,	// 16
			1.,	// 17
			1e-1,	// 18
			1e1,	// 19 *
			1.,	// 20 **
			1e-1,	// 21
			1e-1,	// 22
			1e-1,	// 23

			1e-2, 	// 24 NO pathway
			1e-2,	// 25 
			1e-1,	// 26
			1,	// 27
			1e-1,	// 28
			1e-4,	// 29 
			1e-1,	// 30
			1e-1,	// 31
			1e-1,	// 32
			1e+1};	// 33

    for (int i = 0; i < n; i++) {
 	    for (int la = 0; la < 34; la++) {
                // printf("***** tolerance check: var = %d: %e %e  %e \n", la, x[la], x0[la % 24], fabs(x[la] / x0[la % 24])); // TEST
            }
        smallenough &= (fabs(x[i] / x0[i % 34]) < tol);  // W->nequ hardcoded
        //smallenough &= (fabs(x[i]) < tol);
        //printf("%f \n", x[i]);
        if (!smallenough)
            break;
    }
    return smallenough;

}
css * newton_sparsity(cs *J) {
    // Perform symbolic analysis of the Jacobian for subsequent LU
    // factorisation
    css *S;
    int order = 0;           // 0 = natural, 1 = min deg order of A + A'
    S = cs_sqr(order, J, 0); // 0 means we're doing LU and not QR
    return S;
}
void newton_matrix(odews *ws) {
    // Create a Newton matrix from the given step gamma and Jacobian in W
    cs *M, *eye;
    if (ws->mdeclared)
        cs_nfree(ws->N);
    else
        ws->mdeclared = 1;

    eye = speye(ws->W->J->m);
    M = cs_add(eye, ws->W->J, 1, -ws->gamma);
    cs_spfree(eye);

    ws->N = cs_lu(M, ws->S, 1);
    cs_spfree(M);
}
int lusoln(odews *ws, double *b) {
    // Can only be called if newton_matrix has been called already 
    double *x;
    int n = ws->W->J->n;
    x = cs_malloc (n, sizeof (*x));
    int ok = ws->S && ws->N && x;
    if (ok) {
        cs_ipvec(ws->N->pinv, b, x, n);
        cs_lsolve(ws->N->L, x);
        cs_usolve(ws->N->U, x);
        cs_ipvec(ws->S->q, x, b, n);
    }
    cs_free (x);
    return ok;
}
