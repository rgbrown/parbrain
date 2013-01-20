#ifndef _MATOPS_H
#define _MATOPS_H

#include <cs.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>

typedef struct numjac
{
    cs *A;  /* Jacobian */
    int *g; /* Group */
    int *r; /* Group boundaries */
    int ng; /* Number of groups */
} numjac;

/* Vector routines */
double * onesv(int n);
double * repmatv(double a, int n);
double * zerosv(int n);
double * copyv(double *x, int n);

void dxxy (int n, const double *x, double *y);
double mean(double *x, int n);

/* Sparse matrix creation routines */
cs * spdiags(const double *x, int n);
cs * speye(int n);
cs * spzeros(int m, int n);

/* Matrix manipulation routines */
cs * vertcat(const cs *A, const cs *B);
cs * horzcat(const cs *A, const cs *B);
void cholsoln(csn *Nu, css *S, int n, double *b, double *x);

/* All dense matrices are in COLUMN major format, A[j] is the j-th column
 * of A */
double ** sparse2dense(const cs *A);
cs * dense2sparse(double **A, int m, int n);
void densefree(double **A, int n);
void denseprint(double **A, int m, int n);
void sparseprint(cs *A);
void vecprint(double *v, int n);
cs * matcopy(const cs *A);
cs * subsref(const cs *A, int *Ii, int *Jj, int ni, int nj);

/* Inverse morton transform to get spatial ordering from Morton order */
uint32_t imorton_odd(uint32_t x);
uint32_t imortonx(uint32_t x);
uint32_t imortony(uint32_t x);

/* Numerical Jacobian code */
numjac *numjacinit(cs *A);
cs * mldivide_chol(cs *A, css *S, cs *B);
#endif
