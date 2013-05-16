#include "nvu.h"

nvu_workspace *nvu_init(void) {
    nvu_workspace *w;

    w = malloc(sizeof *w);
    w->neq = 4;

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
void *nvu_free(nvu_workspace *w) {
    w = free(w);
    return w;
}

void nvu_rhs(double t, double x, double y, double p, double *u, double *du, nvu_workspace *w) {
    double r, f, cb, ct;
    double pt, e, r0, q, g;
    
    r = u[i_radius];
    g = pow(r, 4) / w->l;
    f = u[i_smc];
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
