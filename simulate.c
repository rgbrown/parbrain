#include "brain.h"

typedef struct odews {
    workspace *W;
    csn *N; // Newton matrix numeric factorisation
    css *S; // Newton matrix sybolic factorisation
    double gamma;
} odews;

// prototypes
void newton_matrix(odews *ws);
int lusoln(odews *ws, double *b);
css * newton_sparsity(cs *J);
void dcopy(int n, const double *x, double *y);
void daxpy(int n, double a, const double *x, double *y);

int main(int argc, char **argv) {
    double *y;
    odews *ws;
    workspace *W;
    double wt0, wtf;
    
    MPI_Init(&argc, &argv);
    if (W->rank == 0) fprintf(stderr, "Boo! ... ");

    // Initialise the workspace 
    ws = malloc(sizeof (*ws));
    ws->W = init(argc, argv);
    W = ws->W;

    // Set up initial conditions 
    int ny = W->P->ny;
    y = repmatv(2., ny);

    // Initialise vector f for function evaluation
    double *dy;
    dy = repmatv(0., ny);

    // Perform initial evaluation and Jacobian computation
    evaluate(W, 0., y, dy);
    jacupdate(W, 0., y);

    int donewton = 1;
    if (donewton) {
        wt0 = MPI_Wtime();
        // Compute symbolic factorisation of Newton Matrix
        ws->S = newton_sparsity(W->J);

        // Get initial Newton matrix
        ws->gamma = 1e-1;
        newton_matrix(ws);

        // Do simple newton iteration
        double t0 = 0., tf = 10;
        double t = t0;

        
        double *beta, *w, *x;
        double tnext;
        beta = zerosv(ny);
        w = zerosv(ny);
        x = zerosv(ny);
        for (int i = 0; t < tf; i++) {
            //if (W->rank == 0) vecprint(y, ny);
            // w[k+1] = w[k] - M^-1 (w[k] - beta - gamma g(w[k])) 
            // where for Backward Euler, beta = y_i, g(w) = f(w, t_{i+1})

            // copy current y into w and beta (we use current value as initial
            // condition)
            dcopy(ny, y, beta);
            dcopy(ny, y, w);
            tnext = t + ws->gamma;

            for (int k = 0; k < 5; k++) {
                evaluate(W, tnext, w, dy);
                // Do a fixed number of Newton steps
                //if (W->rank == 0) vecprint(w, ny);

                // Form x = w - beta - gamma g
                dcopy(ny, w, x);
                daxpy(ny, -1, beta, x);
                daxpy(ny, -ws->gamma, dy, x);

                // solve
                lusoln(ws, x);

                // update w
                daxpy(ny, -1, x, w);
            }
            dcopy(ny, w, y);
            t += ws->gamma;
            //if (!(i % 10)) vecprint(y, ny);

        }
        wtf = MPI_Wtime();
        if (W->rank == 0) printf("Elapsed time %.2g seconds\n", wtf - wt0);
    }
    MPI_Finalize();
    return 0;
}

void dcopy(int n, const double *x, double *y) {
    for (int i = 0; i < n; i++)
        y[i] = x[i];
}

void daxpy(int n, double a, const double *x, double *y) {
    for (int i = 0; i < n; i++) 
        y[i] += a * x[i];
}


css * newton_sparsity(cs *J) {
    css *S;
    int order = 0;           // 0 = natural order, 1 = minimum degree ordering of A + A'
    S = cs_sqr(order, J, 0); // 0 means we're doing LU and not QR
    return S;
}

void newton_matrix(odews *ws) {
    // Create a Newton matrix from the given step gamma and Jacobian in W
    cs *M, *eye;

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
