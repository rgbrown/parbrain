#include <stdio.h>
#include <cs.h>
#include <math.h>
#include "matops.h"

int *mx, *my, *dx, *dy, *level, *offset;
void ind2xy(int N, int k, int *x, int *y);
int xy2ind(int N, int x, int y, int L);

int main(int argc, char **argv) {
    // Parse inputs
    int N;
    N = argc <= 1 ? 5 : atoi(argv[1]);

    // Compute lookup tables
    int n = (1 << N) - 1;
    level = malloc(n * sizeof(*level));
    for (int k = 0; k < n; k++) {
        level[k] = (int) ceil(log2((double) ((1 << N) - k))) - 1;
    }

    mx = malloc( (N-1) * sizeof (*mx));
    my = malloc( (N-1) * sizeof (*my));
    dx = malloc( (N-1) * sizeof (*dx));
    dy = malloc( (N-1) * sizeof (*dy));
    offset = malloc( (N-1) * sizeof (*offset));

    printf("mx      my      dx      dy      offset\n");
    for (int L = 0; L < N; L++){
        mx[L] = 1 << ((int) floor(((double) L)/2));
        my[L] = 1 << ((int) ceil(((double) L)/2));
        dx[L] = 1 << ((N+1)/2 - (int) (floor(((double) L) / 2)));
        dy[L] = 1 << ((N+1)/2 - (int) (ceil(((double) L) / 2)));
        offset[L] = (1 << N) - (1 << (L+1));
        printf("%-8d%-8d%-8d%-8d%-8d\n", mx[L], my[L], dx[L], dy[L], offset[L]);
    }

    printf("Checking spatial stuff\nk   x   y\n");
    int x, y;
    for (int k = 0; k < n; k++) {
        ind2xy(N, k, &x, &y);
        printf("%-4d%-4d%-4d\n", k, x, y);
    }


    // Now we gotta do the algorithm
    cs *A, *T;
    int *Ti, *Tj;
    double *Tx;
    int m = (1 << (N-1)) - 1;
    T = cs_spalloc(m, n, 3*m, 1, 1);
    Ti = T->i; Tj = T->p; Tx = T->x;
    int xd1, yd1, kd1, xd2, yd2, kd2;
    int idx = 0, row = 0;
    for (int L = N-2; L >= 0; L--) {
        for (int k = offset[L]; k < offset[L] + (1 << L); k ++) {
            ind2xy(N, k, &x, &y);
            if ((N - L) % 2 == 0) {
                xd1 = x - dx[L+1] / 2; yd1 = y;
                xd2 = x + dx[L+1] / 2; yd2 = y;
            } else {
                yd1 = y - dy[L+1] / 2; xd1 = x;
                yd2 = y + dy[L+1] / 2; xd2 = x;
            }
            printf("k = %d, kd1 = %d, kd2 = %d\n", k, kd1, kd2);
            kd1 = xy2ind(N, xd1, yd1, L+1);
            kd2 = xy2ind(N, xd2, yd2, L+1);
            Ti[idx] = row;   Tj[idx] = k;   Tx[idx++] = -1;
            Ti[idx] = row;   Tj[idx] = kd1; Tx[idx++] = 1;
            Ti[idx] = row++; Tj[idx] = kd2; Tx[idx++] = 1;
            
        }
    }
    T->nz = 3*m;
    A = cs_compress(T);
    cs_spfree(T);
    printf("Boo!\n");
    sparseprint(A);

    free(level); free(mx); free(my); free(dx); free(dy); free(offset);

    printf("n: %d\n", n);
    return 0;
}

void ind2xy(int N, int k, int *x, int *y) {
    int L = level[k];
    k = k - offset[L];
    int i = k / my[L];
    int j = k % my[L];
    *x = i * dx[L] - (mx[L] - 1) * dx[L] / 2;
    *y = j * dy[L] - (my[L] - 1) * dy[L] / 2;
}

int xy2ind(int N, int x, int y, int L) {
    int i = (x + (mx[L]-1) * dx[L] / 2) / dx[L];
    int j = (y + (my[L]-1) * dy[L] / 2) / dy[L];
    printf("mx: %d  my: %d  x: %d   y: %d   i: %d   j: %d \n", mx[L], my[L], x, y, i, j);

    int k = j + i*my[L] + offset[L];
    return k;
}



