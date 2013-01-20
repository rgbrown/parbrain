#include <stdio.h>
#include "matops.h"

int main(void)
{
    cs *T, *L, *Lt, *A, *X, *B;
    css *S;
    T = cs_load(stdin);
    L = cs_compress(T);
    Lt = cs_transpose(L, 1);
    A = cs_multiply(L, Lt);
    cs_spfree(T);
    sparseprint(A);
    S = cs_schol(0, A);
    B = speye(A->n);
    X = mldivide_chol(A, S, B);
    sparseprint(X);
}
