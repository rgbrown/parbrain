#include "brain.h"

int main(int argc, char **argv)
{
    workspace *W;
    MPI_Init(&argc, &argv);

    W = init(argc, argv);
    if (W->rank == 0) {
        spy(W->A0);
        printf("\n");
        spy(W->A);

        printf("N:    %d\n", W->N);
        printf("Nsub: %d\n", W->Nsub);
        printf("N0:   %d\n", W->N0);
        printf("Np:   %d\n", W->Np);
    }

    // Test function evaluatios
    double *u;
    double *du;
    u = malloc (W->nu * sizeof(*u));
    du = malloc (W->nu * sizeof(*du));
    for (int i = 0; i < W->nblocks; i++) {
        nvu_ics(u + W->neq*i,0., 0., W->nvu);
    }
    evaluate(W, 0, u, du);
    if (W->rank == 0) {
        for (int i = 0; i < W->nblocks; i++) {
            printf("%16f %16f\n", W->x[i], W->y[i]);
        }
    }

    MPI_Finalize();
    return 0;
}
