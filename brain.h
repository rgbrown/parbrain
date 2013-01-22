#ifndef BRAIN_H
#define BRAIN_H
#include "matops.h"
#include <mpi.h>
#include <stdio.h>
#include <assert.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
extern const int NDEFAULT   ;
extern const int NSUBDEFAULT;
extern const int NSYMBOLS   ;
extern const double RMIN    ; 
extern const double R0      ; 
extern const double L0      ; 
extern const double LRR     ; 
extern const double HRR     ; 
extern const double PROOT   ; 
extern const double PCAP    ; 
extern const double MU      ; 
extern const double RSCALE  ;
extern const double P0      ;
extern const double E0      ;
extern const double EPASSIVE;
extern const double EACTIVE ;
extern const double ETA     ; 
extern const double GAMMA   ; 
extern const double CSTAR   ;
extern const double TAUM    ; 
extern const double PW      ; 
extern const double AVR     ; 
extern const double RCAP    ; 
extern const double CIN     ; 
extern const double C0      ;
extern const double MNORMAL ;
extern const double M0      ; 
extern const double T0      ; 

/* Define the workspace for model-specific parameters, like RHS equations,
 * etc */
typedef struct user {
    int n_auto;
    int ny;
    int i_radius;
    int i_smc;
    int i_ctissue;
    int i_cblood;
    double *l;
    double *a1, *a2, *a3, *a4, *a5;
    double *b1, *cstar, *gamma;
    double *d1, *d2;
    double *g1, *g2;
    double pcap;
} user;

/* Define the workspace that we need to carry around. This is the base
 * workspace for computing the RHS and Jacobian, and performing the solves
 * involved in a Newton step, and as such contains pretty much everything
 * */
typedef struct workspace {
    /* First, the pieces required to describe the parallel tree problem */
    cs      *A;     /* Local adjacency matrix */
    cs      *At;
    cs      *G;     /* Conductance matrix */
    int     *level; /* Level of each vessel in A */
    double  *g;     /* Conductance vector, one per vessel */

    /* Root subtree */
    cs      *A0;    /* Root adjacency matrix */
    cs      *A0t;    
    int     *level0;/* Level of each vessel in A0 */
    cs      *G0;    /* Root conductance matrix */
    double  *g0;    /* Conductance vector */

    /* Vectors to hold local pressures and flows and solutions */
    css     *symbchol;   /* symbolic Cholesky factorisation */
    css     *symbchol0;
    double  *b;     /* RHS vector */
    double  *b0;    /* RHS vector for root */
    double  *p;     /* Pressure at each node (one per row of A) */
    double  *q;     /* Flow through each vessel */
    double  *w;     /* Pressure drop over each vessel */
    double  *p0;    /* Pressures at root */
    double  *q0;    /* Flow at root */
    double  *u;     /* Intermediate variable for parallel flow computation */
    double  *v;     /* Intermediate variable for parallel flow computation */
    double  *ucomm; /* Communication vector */
    double  *vcomm; /* Communication vector */
    double  *gcomm; /* Communication vector */
    double  *xn;    /* extra n-sized workspace */
    double  *xm;    /* extra m-sized workspace */
    double  *xn0;   /* extra n0-sized workspace */

    /* MPI Information */
    int     rank;
    int     n_procs;
    double *buf;    /* Communication buffer */

    /* Geometrical information */
    int     N;      /* Total number of levels */
    int     Nsub;   /* Subtree size for blk-diagonal Jacobian */
    int     N0;     /* Number of root levels */
    int     Np;     /* Number of levels for local subtree */

    /* Jacobian information */
    int    isjac;
    numjac *dfdx;
    numjac *dfdp;
    cs     *dgdx;
    cs     *dpdgneg;
    cs     *A_A;
    cs     *A_At;
    cs     *Proj;
    css    *symbchol_reduced;
    cs     *J;
    cs     *M;

    int *flag; // buffer for communicating convergence information

    /* Model-specific stuff */
    user   *P;
} workspace;

/* Model-specific methods:
 *     p0 defines the time-dependent pressure input
 *     rhs defines the differential equations
 */
double p0(double t);
void rhs(workspace *W, double t, double *y, double *p, double *dy);

/* Methods for external use: essentially function evaluation and Jacobian
 * updates. evaluate, solve and jacupdate both have the side effect of updating
 * conductance, pressure, and flow values in the workspace data structure,
 * so these values should not be relied to remain static
 * */
workspace * init(int argc, char **argv);
void evaluate(workspace *W, double t, double *y, double *dy);
void solve(workspace *W, double pin, double pc);
void jacupdate(workspace *W, double t, double *y);

/* Internal methods: shouldn't be used in code outside brain.c */
void init_geometry(workspace *W, int argc, char **argv);
int is_power_of_two(unsigned int x);
void init_subtree(workspace *W);
cs *adjacency(int N);
void compute_symbchol(workspace *W);
void init_roottree(workspace *W);
void init_problem(workspace *W);
void set_conductance(workspace *W, int isscaled, int computeroot) ;
double compute_length(int level, int n_levels) ;
double compute_radius(int level, int n_levels) ;
user *init_params(workspace *W);
void init_jacobians(workspace *W);
void init_dfdx(workspace *W);
void init_dfdp(workspace *W);
void init_dgdx(workspace *W);
void init_dpdg(workspace *W);

void compute_uv(workspace *W, double pc);
void communicate(workspace *W);
void compute_root(workspace *W, double pin);
void compute_sub(workspace *W, double pin, double pc);

void eval_dfdx(workspace *W, double t, double *y, double *f, double eps);
void eval_dfdp(workspace *W, double t, double *y, double *f, double eps);
void eval_dgdx(workspace *W, double t, double *y);
void eval_dpdg(workspace *W, double t, double *y);
cs * mldivide_chol(cs *A, css *S, cs *B);


#endif
