#include "brain.h"

int main(void)
{
    cs *A, *At, *AAt;
    css *S;
    csn *N;
    double **B;

    printf("\nTesting adjacency\n");
    A = adjacency(4);
    B = sparse2dense(A);
    denseprint(B, A->m, A->n);
    densefree(B, A->n);

    printf("\nTesting A*A'\n");
    At = cs_transpose(A, 1);
    AAt = cs_multiply(A, At);
    B = sparse2dense(AAt);
    denseprint(B, AAt->m, AAt->n);
    densefree(B, AAt->n);

    printf("\nTesting Cholesky with no reordering\n "
           "\t(should have same NZ pattern in LT portion as A*A')\n");
    S = cs_schol(0, AAt);
    N = cs_chol(AAt, S);
    B = sparse2dense(N->L);
    denseprint(B, N->L->m, N->L->n);

    /* Obsolete code 
    printf("\nTesting create_tree\n");
    E = malloc(sizeof (*E));
    E->N = 7;
    E->N_root = 2;
    E->N_parallel = 4;
    E->N_sub = 2;
    create_tree(E);
    printf("\nA = ");
    sparseprint(E->A);
    printf("\nG = ");
    sparseprint(E->G);
    printf("level: \n");
    for (i = 0; i < E->A->n; i++)
        printf("\t%d\n", E->level[i]);

    printf("\nTesting create_reduced_tree\n");
    create_reduced_tree(E);
    printf("\nA_A = ");
    sparseprint(E->A_A);
    */

    return 0;
}
