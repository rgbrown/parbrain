#include "nvu.h"

// Constants
static const double HRR        = 0.1   ;  // Nondimensional (thickness to radius ratio)
static const double RSCALE     = 0.6   ;  // Dimensionless
static const double E0         = 66e3  ;  // Pa
static const double EPASSIVE   = 66e3  ;  // Pa
static const double EACTIVE    = 233e3 ;  // Pa
static const double ETA        = 2.8e2 ;  // Pa s
static const double GAMMA      = 5     ;  // nondim
static const double CSTAR      = 3.5   ;  // nondim
static const double TAUM       = 5     ;  // s
static const double PW         = 5e-5  ;  // m s^-1
static const double AVR        = 5e3   ;  // m^2 m^-3
static const double RCAP       = 5e-6  ;  // m
static const double CIN        = 2.2e-2;  // mol m^-3
static const double C0         = 2.2e-2;  // mol m^-3
static const double MNORMAL    = 4e-3  ;  // mol m^-3 s^-1
static const double M0         = 4e-3  ;  // mol m^-3 s^-1
static const double T0         = 1     ;  // s

static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere
static const int i_cblood  = 1;
static const int i_smc     = 2;
static const int i_ctissue = 3;

// nvu_init: this user-supplied function does any precomputation required
// for the model
nvu_workspace *nvu_init(void) {
    nvu_workspace *w;

    // Specify the sparsity patterns of the equations with respect to
    // pressure and state variables here. An equation is considered to be
    // dependent on pressure if it contains any pressure (transmural or
    // drop) or flow term
    int dfdp_pattern[] = {1,1,0,0}; // neq * 1 column vector
    int dfdx_pattern[] = {1,1,0,0, 0,1,0,1, 1,0,1,0, 0,1,1,1}; // column major, neq*neq

    // Initialise the workspace
    w = malloc(sizeof *w);
    w->neq = 4;

    // Construct sparse matrices containing the sparsity patterns
    // TODO: modify dense2sparse so we can just use two calls to that,
    // rather than this messy code
    //
    // If you just define the integer arrays dfdp_pattern and dfdx_pattern
    // as above, you can leave the following two blocks as they are.
    cs *T;
    T = cs_spalloc(w->neq, 1, 1, 1, 1);
    for (int i = 0; i < w->neq; i++) {
        if (dfdp_pattern[i]) 
            cs_entry(T, i, 0, 1.);
    }
    w->dfdp_pattern = cs_compress(T);
    cs_spfree(T);

    T = cs_spalloc(w->neq, w->neq, 1, 1, 1);
    for (int j = 0; j < w->neq; j++) {
        for (int i = 0; i < w->neq; i++) {
            if (dfdx_pattern[w->neq*j + i])
                cs_entry(T, i, j, 1.);
        }
    }
    w->dfdx_pattern = cs_compress(T);
    cs_spfree(T);


    double Rstar = R0;
    double hstar = HRR * Rstar;
    w->a1 = E0 * T0 * Rstar / (ETA * R0);
    w->a2 = P0 * Rstar * T0 / (ETA * hstar);
    w->a3 = Rstar / R0;
    w->a4 = 1 - RSCALE;
    w->a5 = EACTIVE / EPASSIVE - 1;

    double Lmin = LRR * RMIN;
    double Acap = pow(2*Lmin, 3) / (RCAP / 2. + 1. / AVR);
    double Vt   = Acap / AVR;
    double Vb   = RCAP / 2 * Acap;
    double Q0   = M_PI * pow(R0, 4) * P0 / (8 * MU * L0);

    w->d1 = T0 * Q0 / Vb;
    w->d2 = T0 * PW * Acap / Vb;
    w->g1 = T0 * PW * Acap / Vt;
    w->g2 = T0 * M0 / C0;
    w->b1 = T0 / TAUM;
    w->l  = 1; // normalised away
    w->gamma = GAMMA;
    w->cstar = CSTAR;
    w->pcap  = PCAP / P0;
    
    return w;
}
// This frees the nvu_workspace structure. If you have allocated any memory
// within this struct, here is where you free it
void *nvu_free(nvu_workspace *w) {
    cs_spfree(w->dfdp_pattern);
    cs_spfree(w->dfdx_pattern);
    free(w);
    return w;
}

// right hand side evaluation function. 
//      t       time, 
//      x,y     spatial coordinates, 
//      p       the pressure at the top of the vessel, 
//      u       state variables, the first of which is the vessel radius
//      du      output vector, in the same order (already allocated)
void nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *w) {
    double r, f, cb, ct;
    double pt, e, r0, q, g;
    
    r  = u[i_radius];
    g  = pow(r, 4) / w->l;
    f  = u[i_smc];
    ct = u[i_ctissue];
    cb = u[i_cblood];

    pt = 0.5 * (p + w->pcap);
    q  = (p - w->pcap) * g;
    e  = 1. + w->a5 * f;
    r0 = w->a3 * (1. - w->a4 * f);

    du[i_radius]  = -w->a1 * e * (r / r0 - 1.) + w->a2 * r * pt;
    du[i_smc]     = -w->b1 * (f - 1 / (1 + exp(w->gamma * (ct - w->cstar))));
    du[i_cblood]  =  w->d1 * q * (1 - cb) + w->d2 * (ct - cb);
    du[i_ctissue] = -w->g1 * (ct - cb) + w->g2;
}


// Time-varying pressure at the root of the tree. 1 is nominal value. If
// you want to work in unscaled units, make sure you divide by P0
// afterwards
double nvu_p0(double t) {
    double p0 = 1.5 * 8000 / P0;
    return p0;
}

// Initial conditions. If you want spatial inhomegeneity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *w) {
    u0[i_radius] = 1.;
    u0[i_smc] = 1.;
    u0[i_ctissue] = 1.;
    u0[i_cblood] = 1.;
}
