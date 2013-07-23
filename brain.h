#ifndef BRAIN_H
#define BRAIN_H
#include "matops.h"
#include "nvu.h"
#include <mpi.h>
#include <stdio.h>
#include <assert.h>

#define FILENAMESIZE 128

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
extern const double PROOT   ;
extern const double P0      ;
extern const double PCAP    ;
extern const double MU      ;


// Define the workspace that we need to carry around. This is the base
// workspace for computing the RHS and Jacobian, and performing the solves
// involved in a Newton step, and as such contains pretty much everything
typedef struct workspace {
    // First, the pieces required to describe the parallel tree problem 
    cs      *A;     // Local (per node) adjacency matrix
    cs      *At;    // transpose of A
    cs      *G;     // Conductance matrix 
    int     *level; // Level of each vessel in A
    double  *g;     // Conductance vector, one per vessel. will be a pointer to data in G
    int     nblocks;// Number of nvu blocks
    int     nu;     // Number of equations total
    int     neq;    // Number of equations per block

    // Root subtree
    cs      *A0;    // Root adjacency matrix
    cs      *A0t;    
    int     *level0;// Level of each vessel in A0
    cs      *G0;    // Root conductance matrix
    double  *g0;    // Conductance vector

    // Vectors to hold local pressures and flows and solutions */
    css     *symbchol;  // symbolic Cholesky factorisation 
    css     *symbchol0; // factorisation of root
    double  *l;     // Scaled vessel lengths for regulating vessels 
    double  *b;     // RHS vector 
    double  *b0;    // RHS vector for root 
    double  *p;     // Pressure at each node (one per row of A)
    double  *q;     // Flow through each vessel 
    double  *w;     // Pressure drop over each vessel
    double  *p0;    // Pressures at root 
    double  *q0;    // Flow at root 
    double  *u;     // Intermediate variable for parallel flow computation
    double  *v;     // Intermediate variable for parallel flow computation 
    double  *ucomm; // Communication vector 
    double  *vcomm; // Communication vector 
    double  *gcomm; // Communication vector 
    double  *x;     // x coordinates
    double  *y;     // y coordinates
    double  *xn;    // extra n-sized workspace 
    double  *xm;    // extra m-sized workspace
    double  *xn0;   // extra n0-sized workspace

    // MPI Information 
    int     rank;
    int     n_procs;
    double  *buf;    // Communication buffer
    char    *outfilename;
    MPI_File outfile;

    // Geometrical information
    int     N;      // Total number of levels */
    int     Nsub;   // Subtree size for blk-diagonal Jacobian */
    int     N0;     // Number of root levels */
    int     Np;     // Number of levels for local subtree */

    // Jacobian information 
    int    isjac;
    numjac *dfdx;   // derivatives of DEs with respect to state
    numjac *dfdp;   // derivatives of DEs with respect to pressure
    cs     *dgdx;   // derivative of conductance with respect to state
    cs     *dpdgneg;
    cs     *A_A;
    cs     *A_At;
    cs     *Proj;
    css    *symbchol_reduced;
    cs     *J;
    cs     *M;

    int *flag;      // buffer for communicating convergence information
    int fevals;
    int jacupdates;
    double tfeval;
    double tjacupdate;
    double tjacfactorize;

    /* Model-specific stuff */
    nvu_workspace   *nvu;
} workspace;

/* Methods for external use: essentially function evaluation and Jacobian
 * updates. evaluate, solve and jacupdate both have the side effect of updating
 * conductance, pressure, and flow values in the workspace data structure,
 * so these values should not be relied to remain static
 * */
// defined in adjacency.c
cs *adjacency(int N);

// defined in brain.c
workspace *init(int argc, char **argv);
void    evaluate(workspace *W, double t, double *y, double *dy);
void    solve(workspace *W, double pin, double pc);
void    jacupdate(workspace *W, double t, double *y);
double  p0(double t);

/* Internal methods: shouldn't be used in code outside brain.c */
void    init_parallel(workspace *W, int argc, char **argv);
void    init_io(workspace *W);
void    close_io(workspace *W);
void    write_data(workspace *W, double t, double *y);
int     is_power_of_two(unsigned int x);
void    init_subtree(workspace *W);
//cs     *adjacency(int N);
void    compute_symbchol(workspace *W);
void    init_roottree(workspace *W);
void    set_spatial_coordinates(workspace *W);
void    set_conductance(workspace *W, int isscaled, int computeroot) ;
void    set_length(workspace *W);                  
double  compute_length(int level, int n_levels) ;
double  compute_radius(int level, int n_levels) ;
void    init_jacobians(workspace *W);
void    init_dfdx(workspace *W);
void    init_dfdp(workspace *W);
void    init_dgdx(workspace *W);
void    init_dpdg(workspace *W);

void    compute_uv(workspace *W, double pc);
void    communicate(workspace *W);
void    compute_root(workspace *W, double pin);
void    compute_sub(workspace *W, double pin, double pc);

void    eval_dfdx(workspace *W, double t, double *y, double *f, double eps);
void    eval_dfdp(workspace *W, double t, double *y, double *f, double eps);
void    eval_dgdx(workspace *W, double t, double *y);
void    eval_dpdg(workspace *W, double t, double *y);
void    rhs(workspace *W, double t, double *y, double *p, double *dy);
cs     *mldivide_chol(cs *A, css *S, cs *B);

#endif
