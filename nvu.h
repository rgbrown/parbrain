#ifndef NVU_H
#define NVU_H
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include <cs.h>

// Constants we need from brain. We don't need to explicitly include
// brain.h though
extern const double RMIN;
extern const double R0;
extern const double L0;
extern const double LRR;
extern const double MU;      
extern const double P0;
extern const double PROOT;
extern const double PCAP;

// nvu_workspace gets created once per node. So the parameters therein are
// generic for each nvu. Spatial inhomogeneity should be implemented by
// making the RHS explicitly dependent on the spatial coordinates
typedef struct nvu_workspace {
    // Mandatory fields (accessed externally). Initialised in nvu_init
    int neq;
    cs *dfdp_pattern; // neq * 1 matrix indicating dependence on p
    cs *dfdx_pattern; // neq * neq matrix indicating Jacobian structure of nvu 

    // Parameters that are the same for each NVU. These are model dependent
    double a1, a2, a3, a4, a5;
    double b1, d1, d2, g1, g2;
    double l;
    double gamma, cstar;
    double pcap;

} nvu_workspace;

// Initialisation routine. Gets called once before simulation
nvu_workspace* nvu_init(void); 

// Right hand side routine for one block
void           nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *w);

// Tidy up routine. Free anything allocated in nvu_init here
void           *nvu_free(nvu_workspace *w); 

// Time-varying input pressure function
double          nvu_p0(double t);     

// Initial conditions
void            nvu_ics(double *u0, double x, double y, nvu_workspace *w);


#endif
