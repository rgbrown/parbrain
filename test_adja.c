#include <stdio.h>
#include <math.h>
#include "matops.h"
#include "brain.h"


    int m = (1 << (N-1)) - 1;
    int n = (1 << N) - 1;

    // Set the size of the lowest level grid
    int ncols = 1 << ((N-1)/2); 
    int nrows = 1 << ((N-1)/2 + (N-1)%2);
    
    int a, k = 0;
    int k1, k2, row = 0, col = (1 << (N-1)); 
    int xbranch = 0;
    // int offset = nnodes + 1; // because we have to add the leaf nodes
    int x, y, z, h1, h2;  // TODO: noch ordentlich zu arrays machen!?
    for (int i = 0; i < col; i++) {
	// x[i] = i / 2;
	y[i] = i; //offset + i; // offset can be added later
        h1[i] = 0;
        h2[i] = 0;
    }

    for (int L = N - 1; L > 0; L--) {
        a = (1 << N) - (1 << (L+1));

        if (xbranch) {
            for (int j = 0; j < ncols; j+=2) {
                for (int i = 0; i < nrows; i++) {
                    k1 = a + i + j*nrows;
                    k2 = a + i + (j+1)*nrows;
                    x[k1] = row; x[k2] = row; y[col] = row; 
		    h1[col] = k1; // oder was anderes statt k1 & k2?
		    h2[col] = k2; // s.o.
		    row++; col++;
                }
            }
            ncols /= 2;
        } 
        else {
            for (int j = 0; j < ncols; j++) {
                for (int i = 0; i < nrows; i+=2) {
                    k1 = a + i + j*nrows;
                    k2 = k1 + 1;
		    x[k1] = row; x[k2] = row; y[col] = row; 
		    h1[col] = k1; // s.o.
		    h2[col] = k2;
		    row++; col++; 
                }
            }
            nrows /= 2;
        }
        xbranch = !xbranch;
    } // L loop: from bottom level up to the top of the tree (internal nodes)
    for (int i = 0; i < nnodes; i++) {
        z1[i] = y[h1[i]];
        z2[i] = y[h2[i]];
    }


