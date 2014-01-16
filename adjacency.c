#include <stdio.h>
#include <cs.h>
#include <math.h>
#include "matops.h"
#include "brain.h"

// OLD VERSION
//cs *adjacency(int N) {
//    int m = (1 << (N-1)) - 1;
//    int n = (1 << N) - 1;
//    int *Ti, *Tj;
//    double *Tx;
//
//    cs *T, *A;
//    T = cs_spalloc(m, n, 3*m, 1, 1);
//    if (T == NULL) return T;
//    Ti = T->i; Tj = T->p; Tx = T->x;
//
//    int k = 0;
//    for (int i = 0; i < m; i++) {
//        Ti[k] = i;   Ti[k+1] = i;       Ti[k+2] = i;
//        Tj[k] = 2*i; Tj[k+1] = 2*i + 1; Tj[k+2] = m + i + 1;
//        Tx[k] = 1.;  Tx[k+1] = 1.;      Tx[k+2] = -1.;
//        k += 3;
//    }
//    T->nz = 3*m;
//    A = cs_compress(T);
//    cs_spfree(T);
//    return A;
//}

cs * adjacency(int N) {
    // Initialise the sparse matrix for filling
    
    cs *A, *T;
    int *Ti, *Tj;
    double *Tx;
    int m, n;
    m = (1 << (N-1)) - 1;
    n = (1 << N) - 1;
    T = cs_spalloc(m, n, 3*m, 1, 1);
    Ti = T->i; Tj = T->p; Tx = T->x;

    // Set the size of the lowest level grid
    int ncols = 1 << ((N-1)/2); 
    int nrows = 1 << ((N-1)/2 + (N-1)%2);
    
    int a, k = 0;
    int k1, k2, row = 0, col = (1 << (N-1)); 
    int xbranch = 0;
    for (int L = N - 1; L > 0; L--) {
        a = (1 << N) - (1 << (L+1));
        //b = (1 << N) - (1 << (L  ));
        //c = (1 << N) - (1 << (L-1));

        if (xbranch) {
            for (int j = 0; j < ncols; j+=2) {
                for (int i = 0; i < nrows; i++) {
                    k1 = a + i + j*nrows;
                    k2 = a + i + (j+1)*nrows;
                    Ti[k] = row; Tj[k] = k1; Tx[k++] = 1;
                    Ti[k] = row; Tj[k] = k2; Tx[k++] = 1;
                    Ti[k] = row++; Tj[k] = col++; Tx[k++] = -1;
                }
            }
            ncols /= 2;
        } 
        else {
            for (int j = 0; j < ncols; j++) {
                for (int i = 0; i < nrows; i+=2) {
                    k1 = a + i + j*nrows;
                    k2 = k1 + 1;
                    Ti[k] = row; Tj[k] = k1; Tx[k++] = 1;
                    Ti[k] = row; Tj[k] = k2; Tx[k++] = 1;
                    Ti[k] = row++; Tj[k] = col++; Tx[k++] = -1;
                }
            }
            nrows /= 2;
        }
        xbranch = !xbranch;
    } // L loop: from bottom level up to the top of the tree (internal nodes)
    T->nz = k;
    A = cs_compress(T);
    cs_spfree(T);
    cs_print(A,1);
    return A;
}








