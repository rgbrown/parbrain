#include <stdio.h>
#include "matops.h"

int main(void)
{
    int j;
    int m = 2, n = 5;
    int ni = 3, nj = 2;
    int Ii[] = {0, 2, 4};
    int Jj[] = {1, 2};
    uint32_t i;
    double a;
    double *x, *y, *z;
    double **B;
    numjac *N;
    cs *A, *C, *T;

    /* Test vector construction routines */
    x = onesv(n);
    y = repmatv(-3.0, n);
    z = zerosv(n);

    printf("Testing onesv, repmatv, and zerosv:\n");
    for (i = 0; i < n; i++)
        printf("\t%6.2f %6.2f %6.2f\n", x[i], y[i], z[i]);

    /* Test vector utility routines */
    printf("Testing mean:\n");
    printf("\tx= ");
    for (i = 0; i < n; i++) 
    {
        x[i] = (double) i + 1;
        printf("%6.2f", x[i]);
    }
    a = mean(x, n);
    printf("\n\tmean(x) = %6.3f\n", a);

    printf("\nTesting dxxy: y = repmatv(-3., n), x = 1:n, y = x .* y\n");
    dxxy(n, x, y);
    printf("\ty = ");
    for (i = 0; i < n; i++) printf("%-6.2f", y[i]);
    printf("\n");

    /* Test matrix construction routines */
    printf("\nTesting spzeros:\n");
    A = spzeros(m, n);
    sparseprint(A);
    cs_spfree(A);

    printf("\nTesting speye:\n");
    A = speye(n);
    sparseprint(A);
    cs_spfree(A);


    printf("\nTesting spdiags\n");
    for (i = 0; i < n; i++)
        x[i] = sqrt(i);
    A = spdiags(x, n);
    sparseprint(A);

    printf("\nTesting sparse2dense\n");
    B = sparse2dense(A);
    denseprint(B, n, n);
    cs_spfree(A);

    printf("\nTesting dense2sparse\n");
    A = dense2sparse(B, n, n);
    sparseprint(A);
    densefree(B, n);

    printf("\nTesting horzcat\n");
    C = horzcat(A, A);
    sparseprint(C);
    cs_spfree(C);

    printf("\nTesting vertcat\n");
    C = vertcat(A, A);
    sparseprint(C);
    cs_spfree(C);

    printf("\nTesting subsref\n");
    cs_spfree(A);
    A = speye(6);
    printf("\t A = eye(6)\n");
    sparseprint(A);
    printf("\n\t A([1, 3, 5], :)\n");
    C = subsref(A, Ii, NULL, ni, -1);
    sparseprint(C);
    cs_spfree(C);
    printf("\n\t A(:, [2 3])\n");
    C = subsref(A, NULL, Jj, -1, nj);
    sparseprint(C);
    cs_spfree(C);
    printf("\n\t A([1 3 5], [2 3])\n");
    C = subsref(A, Ii, Jj, ni, nj);
    sparseprint(C);
    cs_spfree(C);
    cs_spfree(A);

    printf("\nTesting inverse Mortons\n");
    for (i = 0; i < 16; i++)
    {
        printf("%4d: %4d %4d\n", i, imortonx(i), imortony(i));
    }

    printf("\nTesting sparse grouping\n");
    T = cs_load(stdin);
    A = cs_compress(T);
    printf("\tT = \n"); sparseprint(A);
    N = numjacinit(A);
    printf(  "\tGrouping : g = ");
    for (j = 0; j < A->n; j++)
        printf("%d ", N->g[j]);
    printf("\n\t           r = ");
    for (j = 0; j < N->ng + 1; j++)
        printf("%d ", N->r[j]);
    printf("\n");
    printf("\t#groups: %d\n", N->ng);

    printf("\nTesting mldivide_chol\n");


    return 0;
}
