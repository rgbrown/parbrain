#include "brain.h"

int main(int argc, char **argv)
{
    workspace *W;
    MPI_Init(&argc, &argv);

    W = init(argc, argv);

    // Test function evaluatios
    double *u;
    double *du;

    u = malloc (W->nu * sizeof(*u));
    du = malloc (W->nu * sizeof(*du));
    for (int i = 0; i < W->nblocks; i++) {
        nvu_ics(u + W->neq*i, W->x[i], W->y[i], W->nvu);
    }
    if (W->rank == 0) {
        for (int i = 0; i < W->nblocks; i++) {
            for (int j = 0; j < W->neq; j++) {
                printf("%-16f", u[W->neq*i + j]);
            }
            printf("\n");
        }
    }

    evaluate(W, 0, u, du);

    MPI_Finalize();
    return 0;
}

