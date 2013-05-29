#include <stdio.h>
#include "matops.h"

int main(void)
{
    cs *T, *A, *B;
    int m = 100, n = 3;
    T = cs_load(stdin);
    A = cs_compress(T);
    B = blkdiag(A, m, n);
    printf("B is %d x %d\n", B->m, B->n);
    sparseprint(B);


    // Test the block diagonal feature

//    cs *T, *L, *Lt, *A, *X, *B;
//    css *S;
//    T = cs_load(stdin);
//    L = cs_compress(T);
//    Lt = cs_transpose(L, 1);
//    A = cs_multiply(L, Lt);
//    cs_spfree(T);
//    sparseprint(A);
//    S = cs_schol(0, A);
//    B = speye(A->n);
//    X = mldivide_chol(A, S, B);
//    sparseprint(X);
}
