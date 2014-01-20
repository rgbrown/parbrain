#include "nvu.h"
// R0, PCAP

// Constants
static const double HRR        = 0.1   ;  // Nondimensional (thickness to radius ratio)
static const double RSCALE     = 0.6   ;  // Dimensionless
static const double E0         = 66e3  ;  // Pa
static const double EPASSIVE   = 66e3  ;  // Pa
static const double EACTIVE    = 233e3 ;  // Pa
static const double ETA        = 2.8e2 ;  // Pa s
//static const double GAMMA      = 5     ;  // nondim
//static const double CSTAR      = 3.5   ;  // nondim
//static const double TAUM       = 5     ;  // s
//static const double PW         = 5e-5  ;  // m s^-1
//static const double AVR        = 5e3   ;  // m^2 m^-3
//static const double RCAP       = 5e-6  ;  // m
//static const double CIN        = 2.2e-2;  // mol m^-3
//static const double C0         = 2.2e-2;  // mol m^-3
//static const double MNORMAL    = 4e-3  ;  // mol m^-3 s^-1
//static const double M0         = 4e-3  ;  // mol m^-3 s^-1
static const double T0         = 1     ;  // s

static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere
static const int ca_i      = 1; // or use enum
static const int ca_sr_i   = 2;
static const int v_i       = 3;
static const int w_i       = 4;
static const int ip3_i     = 5;

static const int ca_j      = 6;
static const int ca_er_j   = 7;
static const int v_j       = 8;
static const int ip3_j     = 9;

static const int Mp        = 10;
static const int AMp       = 11;
static const int AM        = 12;

//static const int i_cblood  = 1;
//static const int i_smc     = 2;
//static const int i_ctissue = 3;

// nvu_init: this user-supplied function does any precomputation required
// for the model
nvu_workspace *nvu_init(void) {
    nvu_workspace *w;

    // Specify the sparsity patterns of the equations with respect to
    // pressure and state variables here. An equation is considered to be
    // dependent on pressure if it contains any pressure (transmural or
    // drop) or flow term
    int dfdp_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,0,0}; // neq * 1 column vector

    int dfdx_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,1,1,1,1,0,1,0,0,0,1,1,1,
                          0,1,1,0,0,0,0,0,0,0,0,0,0,
                          0,1,0,1,0,0,0,0,1,0,0,0,0,
                          0,0,0,1,1,0,0,0,0,0,0,0,0,
                          0,1,0,0,0,1,0,0,0,1,0,0,0,
                          0,1,0,0,0,0,1,1,1,0,0,0,0,
                          0,0,0,0,0,0,1,1,0,0,0,0,0,
                          0,0,0,1,0,0,1,0,1,0,0,0,0,
                          0,0,0,0,0,1,1,0,0,1,0,0,0,
                          0,0,0,0,0,0,0,0,0,0,1,1,0,
                          0,0,0,0,0,0,0,0,0,0,1,1,1,
                          0,0,0,0,0,0,0,0,0,0,1,1,1}; // column major, neq*neq

    // Initialise the workspace
    w = malloc(sizeof *w);
    w->neq = 13;

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

    //double Lmin = LRR * RMIN;
    //double Acap = pow(2*Lmin, 3) / (RCAP / 2. + 1. / AVR);
    //double Vt   = Acap / AVR;
    //double Vb   = RCAP / 2 * Acap;
    //double Q0   = M_PI * pow(R0, 4) * P0 / (8 * MU * L0);

    //w->d1 = T0 * Q0 / Vb;
    //w->d2 = T0 * PW * Acap / Vb;
    //w->g1 = T0 * PW * Acap / Vt;
    //w->g2 = T0 * M0 / C0;
    //w->b1 = T0 / TAUM;
    //w->l  = 1; // normalised away
    //w->gamma = GAMMA;
    //w->cstar = CSTAR;
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

// SMC constants:
    
    const double Fmax_i		= 0.23;		// (microM/s)
    const double Kr_i 		= 1; 		// (microM) ; Half saturation constant for agonist-dependent Ca$^{2+}$ entry
    const double G_Ca		= 0.00129;	// (microM/mV/s)
    const double v_Ca1		= 100;		// (mV)
    const double v_Ca2		= -24;		// (mV)
    const double R_Ca		= 8.5;		// (mV)
    const double G_NaCa		= 0.00316;	// (microM/mV/s)
    const double c_NaCa		= 0.5;		// (microM)
    const double v_NaCa		= -30;
    const double B_i		= 2.025;
    const double cb_i		= 1;
    const double C_i		= 55;
    const double sc_i		= 2;
    const double cc_i		= 0.9;
    const double D_i		= 0.24;
    const double vd_i		= -100;
    const double Rd_i		= 250;
    const double L_i		= 0.025;
    const double gam		= 1970;
    const double F_NaK		= 0.0432;
    const double G_Cl		= 0.00134;
    const double v_Cl		= -25;
    const double G_K		= 0.00446;
    const double vK_i		= -94;
    const double lam 		= 45;
    const double c_w		= 0;
    const double bet		= 0.13;
    const double v_Ca3		= -27;
    const double R_K		= 12;
    const double k_i		= 0.1;
    
    const double Fmax_j		= 0.23;		// [microM/s]
    const double Kr_j		= 1;
    const double B_j 		= 0.5;
    const double cb_j		= 1;
    const double C_j		= 5;
    const double sc_j		= 2;
    const double cc_j		= 0.9;
    const double D_j		= 0.24;
    const double L_j		= 0.025;
    const double G_cat 		= 0.66e-3; //!
    const double E_Ca		= 50;
    const double m3cat		= -0.18; //-6.18;		// changed value!
    const double m4cat 		= 0.37;
    const double JO_j 		= 0.029; //constant Ca influx (EC)
    const double C_m 		= 25.8;
    const double G_tot		= 6927;
    const double vK_j 		= -80;
    const double a1			= 53.3;
    const double a2			= 53.3;
    const double b			= -80.8;
    const double c 			= -6.4; //-0.4; 		// changed value!
    const double m3b		= 1.32e-3;
    const double m4b		= 0.3;
    const double m3s		= -0.28;
    const double m4s		= 0.389;
    const double G_R		= 955;
    const double v_rest		= -31.1;
    const double k_j		= 0.1;
    
    const double J_PLC 		= 0.4; //0.1;
    
    const double g_hat      = 0.5;
    const double p_hat      = 0.05;
    const double p_hatIP3   = 0.05;

    //// Myosin crossbridge model
    const double C_Hillmann = 1;
    const double K2_c        = 0.5 * C_Hillmann;
    const double K3_c        = 0.4 * C_Hillmann;
    const double K4_c        = 0.1 * C_Hillmann;
    const double K5_c        = 0.5 * C_Hillmann;
    const double K7_c        = 0.1 * C_Hillmann;
    const double gam_cross   = 17 * C_Hillmann;

    //// Kelvin Voigt
    const double P_r         = 4000;  // Pa
    const double rb_r        = 20e-6; // m
    const double h0_r        = 3e-6;  // m
    const double R0pas_r     = 20e-6;
    const double R0act_r     = 12e-6;
    const double Eact_r      = 233e3;
    const double Epas_r      = 66e3;
    const double nu_r        = 1e4;

    //TODO: initialise variables!
    double r, f; //, cb, ct;
    double pt, e, r0; //, q, g;
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i, state_ca_j, state_ca_er_j, state_v_j, state_ip3_j;
    double state_Mp, state_AMp, state_AM;
    double flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i, flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i, flu_degrad_i;
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j;
    double flu_K1_c, flu_K6_c, flu_M, flu_h_r, flu_F_r, flu_E_r, flu_R0_r;

    r  = u[i_radius];
    //g  = pow(r, 4) / w->l;
    //f  = u[i_smc];
    //ct = u[i_ctissue];
    //cb = u[i_cblood];

    state_ca_i    = u[ca_i];
    state_ca_sr_i = u[ca_sr_i];
    state_v_i     = u[v_i];
    state_w_i     = u[w_i];
    state_ip3_i   = u[ip3_i];
    
    state_ca_j    = u[ca_j];
    state_ca_er_j = u[ca_er_j];
    state_v_j     = u[v_j];
    state_ip3_j   = u[ip3_j];
    
    state_Mp      = u[Mp];
    state_AMp     = u[AMp];
    state_AM      = u[AM];

// Fluxes:
    //q  = (p - w->pcap) * g;
    f = state_AMp + state_AM;
    e  = 1. + w->a5 * f;
    r0 = w->a3 * (1. - w->a4 * f);
    // pressure:
    pt = 0.5 * (p + w->pcap);

    // SMC
    flu_v_cpl_i		    = - g_hat * ( state_v_i - state_v_j );
    flu_c_cpl_i         = - p_hat * ( state_ca_i - state_ca_j);
    flu_I_cpl_i         = - p_hatIP3 * ( state_ip3_i- state_ip3_j );
    flu_rho_i		    = 1; //pow((K_d + state_ca_i),2) / ( pow((K_d + state_ca_i),2) + ( K_d * B_T ) );
    flu_ip3_i		    = Fmax_i *  pow(state_ip3_i,2) / ( pow(Kr_i,2) + pow(state_ip3_i,2) );
    flu_SRuptake_i      = B_i * ( pow(state_ca_i,2) ) / ( pow(state_ca_i,2) + pow(cb_i,2) );
    flu_CICR_i		    = C_i *  ( pow(state_ca_sr_i,2) ) / ( pow(sc_i,2) + pow(state_ca_sr_i,2) ) *  ( pow(state_ca_i,4) ) / ( pow(cc_i,4) + pow(state_ca_i,4) );
    flu_extrusion_i	    = D_i * state_ca_i * (1 + ( state_v_i - vd_i ) / ( Rd_i ) );
    flu_leak_i 		    = L_i * state_ca_sr_i;
    flu_VOCC_i		    = G_Ca * ( state_v_i - v_Ca1) / ( 1 + exp( - ( state_v_i - v_Ca2 ) / ( R_Ca ) ) );
    flu_NaCa_i		    = G_NaCa * ( state_ca_i ) * ( state_v_i - v_NaCa ) / ( state_ca_i + c_NaCa ) ;
    flu_NaK_i		    = F_NaK;
    flu_Cl_i		    = G_Cl * (state_v_i - v_Cl);
    flu_K_i			    = G_K * state_w_i * ( state_v_i - vK_i );
    flu_Kactivation_i   = pow((state_ca_i + c_w),2) / ( pow((state_ca_i + c_w),2) + bet*exp(-(state_v_i - v_Ca3)/R_K) );
    flu_degrad_i	    = k_i * state_ip3_i;
    
    flu_K1_c            = gam_cross * pow(state_ca_i,3);
    flu_K6_c            = flu_K1_c;
    flu_M               = 1 - state_Mp - state_AM - state_AMp;
    //flu_h_r             = 0.1 * state_r;
    flu_h_r             = 0.1 * r;

    flu_F_r = state_AMp + state_AM;
    flu_E_r = Epas_r + flu_F_r * (Eact_r -Epas_r);
    flu_R0_r = R0pas_r + flu_F_r * R0pas_r*(0.6 - 1);



    // EC:
    flu_v_cpl_j			= - g_hat * ( state_v_j - state_v_i );
    flu_c_cpl_j			= - p_hat * ( state_ca_j - state_ca_i );
    flu_I_cpl_j			= - p_hatIP3 * ( state_ip3_j - state_ip3_i );
    flu_rho_j 			= 1; //pow(( K_d + state_ca_j ),2) / ( pow(( K_d + state_ca_j ),2) + ( K_d * B_T ) );
    flu_O_j 			= JO_j;
    flu_ip3_j			= Fmax_j * ( pow(state_ip3_j,2) ) / ( pow(Kr_j,2) + pow(state_ip3_j,2) );
    flu_ERuptake_j      = B_j * ( pow(state_ca_j,2) ) / ( pow(state_ca_j,2) + pow(cb_j,2) );
    flu_CICR_j			= C_j *  ( pow(state_ca_er_j,2) ) / ( pow(sc_j,2) + pow(state_ca_er_j,2) ) *  ( pow(state_ca_j,4) ) / ( pow(cc_j,4) + pow(state_ca_j,4) );
    flu_extrusion_j     = D_j * state_ca_j;
    flu_leak_j          = L_j * state_ca_er_j;
    
    flu_cation_j 		= G_cat * ( E_Ca - state_v_j) * 0.5 * ( 1 + tanh( ( log10( state_ca_j ) - m3cat )  /  m4cat  ) );
    //flu_cation_j      = G_cat * 0.5 * ( E_Ca - state_v_j ) * ( 1 + tanh( ( log10( state_ca_j ) - m3cat ) /  m4cat  ) );
    
    flu_BKCa_j 			= 0.4/2 * ( 1 + tanh( ( (  log10(state_ca_j) - c) * ( state_v_j - b ) - a1 ) / ( m3b* pow(( state_v_j + a2 * ( log10( state_ca_j ) - c ) - b),2) + m4b ) ) );
    
    flu_SKCa_j 			= 0.6/2 * ( 1 + tanh( ( log10(state_ca_j) - m3s ) /  m4s ));
    flu_K_j 			= G_tot * ( state_v_j - vK_j ) * ( flu_BKCa_j + flu_SKCa_j ); // Reihenfolge!!
    flu_R_j 			= G_R * ( state_v_j - v_rest);
    flu_degrad_j 		= k_j * state_ip3_j;
    


// Differential Equations:
    // Radius:
    du[i_radius]  = -w->a1 * e * (r / r0 - 1.) + w->a2 * r * pt;
    //du[i_smc]     = -w->b1 * (f - 1 / (1 + exp(w->gamma * (ct - w->cstar))));
    //du[i_cblood]  =  w->d1 * q * (1 - cb) + w->d2 * (ct - cb);
    //du[i_ctissue] = -w->g1 * (ct - cb) + w->g2 + exp(-(x*x + y*y)/0.01);
    //du[i_radius]  = R0pas_r/nu_r *(state_r * P_r / flu_h_r - flu_E_r * (state_r - flu_R0_r) / flu_R0_r); //

    //SMC:
    du[ca_i]      = flu_c_cpl_i + flu_rho_i * ( flu_ip3_i - flu_SRuptake_i + flu_CICR_i - flu_extrusion_i + flu_leak_i - flu_VOCC_i + flu_NaCa_i ) ;
    du[ca_sr_i]   = flu_SRuptake_i - flu_CICR_i - flu_leak_i ;
    du[v_i]       = flu_v_cpl_i + gam * ( - flu_NaK_i - flu_Cl_i - 2*flu_VOCC_i - flu_NaCa_i - flu_K_i ) ;
    du[w_i]       = lam * (flu_Kactivation_i - state_w_i ) ;
    du[ip3_i]     = flu_I_cpl_i - flu_degrad_i ;
    
    //EC:
    du[ca_j]      = flu_c_cpl_j + flu_rho_j * ( flu_ip3_j - flu_ERuptake_j + flu_CICR_j - flu_extrusion_j + flu_leak_j + flu_cation_j + flu_O_j ) ;
    du[ca_er_j]   = flu_ERuptake_j - flu_CICR_j - flu_leak_j ;
    du[v_j]       = flu_v_cpl_j - 1/C_m * ( flu_K_j + flu_R_j ) ;
    du[ip3_j]     = flu_I_cpl_j + J_PLC - flu_degrad_j ;

    du[Mp]        = K4_c * state_AMp + flu_K1_c * flu_M - (K2_c + K3_c) * state_Mp;
    du[AMp]       = K3_c * state_Mp + flu_K6_c *state_AM  - (K4_c + K5_c) * state_AMp;  // K7_c was corrected to K4_c
    du[AM]        = K5_c * state_AMp - (K7_c + flu_K6_c) * state_AM;
}



// Time-varying pressure at the root of the tree. 1 is nominal value. If
// you want to work in unscaled units, make sure you multiply by P0
// afterwards
double nvu_p0(double t) {
    double p0 = (0.5 * sin(t) + 1) * 8000 / P0; //1.5 * 8000 / P0;// sin(t); // * 8000 / P0; 
    return p0;
}

// Initial conditions. If you want spatial inhomegenity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *w) {
    u0[i_radius] = 1.*(0.5*sin(x)+1);
    //u0[i_smc] = 1.;
    //u0[i_ctissue] = 1.;
    //u0[i_cblood] = 1.;
    u0[ca_i]      = 0.26;
    u0[ca_sr_i]   = 1;
    u0[v_i]       = -40;
    u0[w_i]       = 0.1;
    u0[ip3_i]     = 0.5;
    
    u0[ca_j]      = 1;
    u0[ca_er_j]   = 0.6;
    u0[v_j]       = -60;
    u0[ip3_j]     = 1.5;
    
    u0[Mp]        = 0.25;
    u0[AMp]       = 0.25;
    u0[AM]        = 0.25;
}
