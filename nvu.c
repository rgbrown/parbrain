// stretch
#include "nvu.h"

// Constants
static const double HRR        = 0.1   ;  // Nondimensional (thickness to radius ratio)
static const double RSCALE     = 0.6   ;  // Dimensionless
static const double E0         = 66e3  ;  // Pa
static const double EPASSIVE   = 66e3  ;  // Pa
static const double EACTIVE    = 233e3 ;  // Pa
static const double ETA        = 2.8e2 ;  // Pa s
static const double T0         = 1     ;  // s

static const int i_radius  = 0; // radius has to be 0, this is assumed elsewhere
// AC
static const int R_k       = 1;
static const int N_Na_k    = 2;
static const int N_K_k     = 3;
static const int N_HCO3_k  = 4;
static const int N_Cl_k    = 5;
static const int N_Na_s    = 6;
static const int N_K_s     = 7;
static const int N_HCO3_s  = 8;
static const int K_p       = 9;
static const int w_k       = 10;
// SMC
static const int ca_i      = 11; // or use enum
static const int ca_sr_i   = 12;
static const int v_i       = 13;
static const int w_i       = 14;
static const int ip3_i     = 15;
static const int K_i       = 16;
// EC
static const int ca_j      = 17;
static const int ca_er_j   = 18;
static const int v_j       = 19;
static const int ip3_j     = 20;
// Mech
static const int Mp        = 21;
static const int AMp       = 22;
static const int AM        = 23;
// NO pathway
static const int NOi       = 24;
static const int NOj       = 25;
static const int NOn       = 26;
static const int cGMP      = 27;
static const int eNOS      = 28;
static const int nNOS      = 29;
static const int ca_n      = 30;
static const int E_b       = 31;
static const int E_6c      = 32;
static const int E_5c      = 33;





// nvu_init: this user-supplied function does any precomputation required
// for the model
nvu_workspace *nvu_init(void) {
    nvu_workspace *w;
//TODO: Sparsity patterns!
    // Specify the sparsity patterns of the equations with respect to
    // pressure and state variables here. An equation is considered to be
    // dependent on pressure if it contains any pressure (transmural or
    // drop) or flow term
    int dfdp_pattern[] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // neq * 1 column vector

    int dfdx_pattern[] =     {1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
    		              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // column major, neq*neq

    // Initialise the workspace
    w = malloc(sizeof *w);
    w->neq = 34;  // TODO: make sure you have the right no of DEs!

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
    w->pcap  = PCAP / P0;
    w->l  = 1; // normalised away

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
// general constants:
	const double Farad       = 96500         ;// [C mol-1]      Faradays constant
	const double R_gas       = 8.315         ;// [J mol-1K-1]
	const double Temp        = 300           ;// [K]
        const double unitcon     = 1e3        ;// [-]            Factor to convert equations to another unit

// NE & AC constants:
    const double L_p         = 2.1e-9        ;// [m uM-1s-1]
    const double R_tot       = 8.79e-8       ;// [m]
    const double X_k         = 12.41e-3      ;// [uMm]
    const double z_Na        = 1             ;// [-]
    const double z_K         = 1             ;// [-]
    const double z_Cl        = -1            ;// [-]
    const double z_NBC       = -1            ;// [-]
    const double g_K_k       = 40            ;// [ohm-1m-2]
    const double g_KCC1_k    = 1e-2          ;// [ohm-1m-2]
    const double g_NBC_k     = 7.57e-1       ;// [ohm-1m-2]
    const double g_Cl_k      = 8.797e-1      ;// [ohm-1m-2]
    const double g_NKCC1_k   = 5.54e-2       ;// [ohm-1m-2]
    const double g_Na_k      = 1.314         ;// [ohm-1m-2]
    const double J_NaK_max   = 1.42e-3       ;// [uMm s-1]
    const double K_Na_k      = 10e3          ;// [uM]
    const double K_K_s       = 1.5e3         ;// [uM]
    const double k_C         = 7.35e-5       ;// [muM s-1]

// BK channel constants:
    const double A_ef_k      = 3.7e-9        ; // m2       Area of an endfoot of an astrocyte, equal to Area astrocyte at synaptic cleft
    const double v_6         = 22e-3         ; // V        Constant
    const double v_4         = 14.5e-3       ; // V        A measure of the spread of the distrubution
    const double psi_w       = 2.664         ; // s-1      A characteristic time
    const double G_BK_k      = 4.3e3         ; // pS      Constant estimation based on Ermentrout
    const double g_BK_k      = G_BK_k * 1e-12 / A_ef_k ;// ohm-1m-2  Specific capacitance of the BK-Channel in units of Ostby
    const double VR_pa       = 0.001       ; // [-]       The estimated volume ratio of perivascular space to astrocyte: Model estimation
    const double VR_ps       = 0.001       ; // [-]       The estimated volume ratio of perivascular space to SMC: Model Estimation


// SMC constants:
    const double F_il = 7.5e2            ;//[-] scalingsfactor to fit the experimental data of Filosa
    const double z_1 =4.5                ;//[-] parameter fitted on experimental data of Filosa
    const double z_2 =-1.12e2            ;//[-] parameter fitted on experimental data of Filosa
    const double z_3 =4.2e-1             ;//[-] parameter fitted on experimental data of Filosa
    const double z_4 =-1.26e1            ;//[-] parameter fitted on experimental data of Filosa
    const double z_5 =-7.4e-2; //[-] parameter fitted on experimental data of Filosa
    const double Fmax_i		= 0.23;		// (microM/s)
    const double Kr_i 		= 1; 		// (microM) ; Half saturation constant for agonist-dependent Ca entry
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
    const double gam		= 1970; // mVmicroM-1 The change in membrane potential by a scaling factor
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
//    const double K_d         = 1;            // = 1e3 nM Gonzalez
//    const double B_T         = 100;          // = 1e5 nM Gonzalez
//    const double Kinf_i      = 1e5;          // 100 mM K+ concentration in SMC

    const double G_stretch   = 0.0061;       // uM mV-1 s-1   (stretch activated channels)
    const double P_str       = 30;
    const double Esac        = -18;          // mV
    const double alpha1      = 0.0074;
    const double sig0        = 500;



// EC constants:
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
    const double c 			= -0.4; //-6.4; 		// changed value!
    const double m3b		= 1.32e-3;
    const double m4b		= 0.3;
    const double m3s		= -0.28;
    const double m4s		= 0.389;
    const double G_R		= 955;
    const double v_rest		= -31.1;
    const double k_j		= 0.1;

    const double J_PLC 		= 0.18; //0.18

    const double g_hat      = 0.5;
    const double p_hat      = 0.05;
    const double p_hatIP3   = 0.05;
    const double C_Hillmann = 1;
    const double K2_c        = 0.5 * C_Hillmann;
    const double K3_c        = 0.4 * C_Hillmann;
    const double K4_c        = 0.1 * C_Hillmann;
    const double K5_c        = 0.5 * C_Hillmann;
    const double K7_c        = 0.1 * C_Hillmann;
    const double gam_cross   = 17 * C_Hillmann;

// NO pathway **********

    const double LArg = 100;
    const double v_spine = 8e-8;
    const double k_ex        = 1600;
    const double Ca_rest     = 0.1; 
    const double lambda      = 20; 
    const double V_maxNOS    = 25e-3; 
    const double V_nNOS      = 1.435;
    const double K_actNOS    = 9.27e-2; 
    const double D_NO 		= 3300;	
    const double dist_ni		= 50;
    const double k_O2        = 9.6e-6; 
    const double On          = 200;    
    const double v_n         = -40;
    const double G_M         = 46;   
    const double dist_ji		= 3.75;	
    const double tau_ni      = pow(dist_ni,2)/(2*D_NO); 
    const double tau_ji      = pow(dist_ji,2)/(2*D_NO); 
    const double P_Ca_P_M    = 3.6; 
    const double Ca_ex       = 2e3;  
    const double M           = 1.3e5;   
    const double betA        = 650 ;
    const double betB        = 2800 ;  
    const double Q1          = 1.9e5; 
    const double Q2          = 2.1e5; 
    const double Q3          = 0.4e5;
    const double Q4          = 0.26e5; 
    const double CaM_thresh = Ca_rest/( Ca_rest*(Q1 + 2*Q1*Q2*Ca_rest + 3*Q1*Q2*Q3* pow(Ca_rest,2) + 4*Q1*Q2*Q3*Q4* pow(Ca_rest,3))/ (1 + Q1*Ca_rest + Q1*Q2* pow(Ca_rest,2) + Q1*Q2*Q3* pow(Ca_rest,3) + Q1*Q2*Q3*Q4* pow(Ca_rest,4)) ); // 2.7764e-5; 
    const double Oj          = 200; 
    const double K_dis       = 9e-2;  
    const double K_eNOS      = 4.5e-1; 
    const double mu2         = 0.0167; 
    const double g_max       = 0.3;    	
    const double alp         = 2;    
    const double W_0         = 1.4;        
    const double delt_wss    = 2.86 ;  
    const double V_eNOS      = 0.24;   
    const double k_dno       = 0.01;  
    const double k1          = 2e3 ;  
    const double k2          = 0.1;   
    const double k3          = 3;  
    const double k_1         = 100;  
    const double V_max_sGC   = 0.8520;  //\muM s{-1}; (for m = 2)
    const double k_pde       = 0.0195;// s{-1} (for m = 2)
    const double C_4         = 0.011; // [s{-1} microM{-2}] (note: the changing units are correct!) (for m = 2)
    const double K_m_pde     = 2;           		// [microM]
    const double R_Kfit      = 30.8; //55;                  // [mV]
    const double V_cGMP      = 66.9; // 68 ;               // [mV]
    const double V_NO        = 100;                 // [mV]
    const double V_b         = 283.7; // 215; (Hannah)                 // [mV]
    const double K_m_cGMP    = 0.55;       		// [microM]
    const double K_m_NO      = 0.2;     		    // [microM]
    const double k_mlcp_b    = 0.0086;         // [s{-1}]
    const double k_mlcp_c    = 0.0327;          //[s{-1}]
    const double K_m_mlcp    = 5.5;        		// [microM]
    const double k_mlck      = 1180;
    const double v_Ca        = -27; // cGMP and NO dependent
    const double c_wi        = 0; // translation factor for Ca dependence of KCa channel activation sigmoidal [microM] 
    const double bet_i       = 0.13; // translation factor for membrane potential dependence of KCa channel activation sigmoidal [microM2] 
    const double m =2;
    const double tau_w       = 1; // will be dynamic later


//TODO: initialise variables!
    double state_r, f; //, cb, ct;
    double pt, e, r0, q, g;

    double state_R_k,   state_N_Na_k, state_N_K_k, state_N_HCO3_k, state_N_Cl_k, state_N_Na_s, state_N_K_s, state_N_HCO3_s, state_K_p, state_w_k; // AC state
    double state_ca_i, state_ca_sr_i, state_v_i, state_w_i, state_ip3_i, state_K_i; // SMC state
    double state_ca_j, state_ca_er_j, state_v_j, state_ip3_j; // EC state
    double state_Mp, state_AM, state_AMp; // Mech state
    double state_ca_n, state_nNOS, state_NOn, state_NOi, state_E_b, state_E_6c, state_E_5c, state_cGMP, state_eNOS, state_NOj; // NO pathway state


    double flu_R_s, flu_N_Cl_s, flu_Na_k, flu_K_k, flu_HCO3_k, flu_Cl_k, flu_Na_s, flu_K_s, flu_HCO3_s, flu_Cl_s, flu_E_Na_k, flu_E_K_k, flu_E_Cl_k, flu_E_NBC_k, flu_E_BK_k, flu_J_NaK_k, flu_v_k, flu_J_KCC1_k, flu_J_NBC_k, flu_J_NKCC1_k, flu_J_Na_k, flu_J_K_k, flu_J_BK_k, flu_w_inf, flu_phi_w; // AC fluxes
    double flu_M, flu_E_K_i, flu_h_r, flu_v_cpl_i, flu_c_cpl_i, flu_I_cpl_i, flu_rho_i, flu_ip3_i, flu_SRuptake_i,
    flu_CICR_i, flu_extrusion_i, flu_leak_i, flu_VOCC_i, flu_NaCa_i, flu_NaK_i, flu_Cl_i, flu_K_i, flu_Kactivation_i,
    flu_degrad_i, flu_v_KIR_i, flu_G_KIR_i, flu_J_KIR_i, flu_J_stretch_i; // SMC fluxes
    double flu_v_cpl_j, flu_c_cpl_j, flu_I_cpl_j, flu_rho_j, flu_O_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j,
    flu_extrusion_j, flu_leak_j, flu_cation_j, flu_BKCa_j, flu_SKCa_j, flu_K_j, flu_R_j, flu_degrad_j, flu_J_stretch_j; // EC fluxes
    double flu_K1_c, flu_K6_c; // Mech fluxes   
    double flu_P_NR2AO, flu_P_NR2BO, flu_openProbTerm, flu_I_Ca, flu_phi_N, flu_dphi_N, flu_N, flu_CaM, flu_W_tau_w, flu_F_tau_w, flu_k4, flu_R_cGMP1, flu_R_NO, flu_v_Ca3, flu_P_O, flu_R_cGMP2, flu_K2_c, flu_K5_c, flu_kmlcp;    // NO pathway fluxes


// BEGIN PARSER
// State Variables:
    state_r  = u[i_radius];

    state_ca_i    = u[ca_i];
    state_ca_sr_i = u[ca_sr_i];
    state_v_i     = u[v_i];
    state_w_i     = u[w_i];
    state_ip3_i   = u[ip3_i];
    state_K_i     = u[K_i];

    state_R_k     = u[R_k];
    state_N_Na_k  = u[N_Na_k];
    state_N_K_k   = u[N_K_k];
    state_N_HCO3_k= u[N_HCO3_k];
    state_N_Cl_k  = u[N_Cl_k];
    state_N_Na_s  = u[N_Na_s];
    state_N_K_s   = u[N_K_s];
    state_N_HCO3_s= u[N_HCO3_s];
    state_K_p     = u[K_p];
    state_w_k     = u[w_k];


    state_ca_j    = u[ca_j];
    state_ca_er_j = u[ca_er_j];
    state_v_j     = u[v_j];
    state_ip3_j   = u[ip3_j];

    state_Mp      = u[Mp];
    state_AMp     = u[AMp];
    state_AM      = u[AM];

// NO pathway
    state_ca_n    = u[ca_n];
    state_nNOS    = u[nNOS];
    state_NOn     = u[NOn];         
    state_NOi     = u[NOi];
    state_E_b     = u[E_b];
    state_E_6c    = u[E_6c];
    state_E_5c    = u[E_5c];
    state_cGMP    = u[cGMP];
    state_eNOS    = u[eNOS];
    state_NOj     = u[NOj];





// end


// Fluxes:
    f = state_AMp + state_AM;
    g  = pow(state_r, 4) / w->l;
    // pressure:
    pt = 0.5 * (p + w->pcap);
    q  = (p - w->pcap) * g; //brauchen wir das noch?
    e  = 1. + w->a5 * f;
    r0 = w->a3 * (1. - w->a4 * f);
    // AC fluxes
    flu_R_s            = R_tot - state_R_k;                            //
    flu_N_Cl_s         = state_N_Na_s + state_N_K_s - state_N_HCO3_s;  //
    flu_Na_k           = state_N_Na_k / state_R_k;                     //
    flu_K_k            = state_N_K_k / state_R_k;                      //
    flu_HCO3_k         = state_N_HCO3_k / state_R_k;                   //
    flu_Cl_k           = state_N_Cl_k / state_R_k;                     //
    flu_Na_s           = state_N_Na_s / flu_R_s;                       //
    flu_K_s            = state_N_K_s / flu_R_s;                        //
    flu_HCO3_s         = state_N_HCO3_s / flu_R_s;                     //
    flu_Cl_s           = flu_N_Cl_s / flu_R_s;                         //
    flu_E_Na_k         = (R_gas * Temp) / (z_Na * Farad) * log(flu_Na_s / flu_Na_k);    // V
    flu_E_K_k          = (R_gas * Temp) / (z_K  * Farad) * log(flu_K_s / flu_K_k );     // V
    flu_E_Cl_k         = (R_gas * Temp) / (z_Cl * Farad) * log(flu_Cl_s / flu_Cl_k);    // V
    flu_E_NBC_k        = (R_gas * Temp) / (z_NBC* Farad) * log((flu_Na_s * pow(flu_HCO3_s,2))/(flu_Na_k * pow(flu_HCO3_k,2)));     // V
    flu_E_BK_k         = (R_gas * Temp) / (z_K  * Farad) * log(state_K_p / flu_K_k);   // V
    //flu_J_NaK_k        = J_NaK_max * Hill(flu_Na_k, K_Na_k, 1.5) * Hill(flu_K_s, K_K_s, 1);    // uMm s-1
    flu_J_NaK_k        = J_NaK_max * ( pow(flu_Na_k,1.5) / ( pow(flu_Na_k,1.5) + pow(K_Na_k,1.5) ) ) * ( flu_K_s / (flu_K_s + K_K_s) );    // uMm s-1  
    flu_v_k            = ( g_Na_k * flu_E_Na_k + g_K_k * flu_E_K_k + g_Cl_k  * flu_E_Cl_k + g_NBC_k * flu_E_NBC_k - flu_J_NaK_k * Farad/unitcon + g_BK_k * state_w_k * flu_E_BK_k) / (g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k * state_w_k);  // V
    flu_J_KCC1_k       = 0 * (R_gas * Temp * g_KCC1_k) / (pow(Farad,2)) * log(((flu_K_s) * (flu_Cl_s))/((flu_K_k)*(flu_Cl_k))) * unitcon;   //uMm s-1
    flu_J_NBC_k        = g_NBC_k / Farad * ((flu_v_k) - (flu_E_NBC_k))*unitcon;       //uMm s-1
    flu_J_NKCC1_k      = 0 * (g_NKCC1_k * R_gas * Temp) / (pow(Farad,2))  * log(((flu_K_s) * (flu_Na_s) * pow(flu_Cl_s,2)) /((flu_K_k) * (flu_Na_k) * pow(flu_Cl_k,2)))*unitcon;        //uMm s-1
    flu_J_Na_k   = g_Na_k / Farad * (flu_v_k - flu_E_Na_k) * unitcon;              //uMm s-1
    flu_J_K_k    = g_K_k  / Farad * ((flu_v_k) - (flu_E_K_k )) * unitcon;          //uMm s-1
    flu_J_BK_k   = g_BK_k / Farad * state_w_k * (flu_v_k - flu_E_BK_k) * unitcon;  //uMm s-1
    flu_w_inf    = 0.5 * (1+tanh(((flu_v_k)+v_6)/v_4));                            //[-]
    flu_phi_w    = psi_w * cosh(((flu_v_k)+v_6)/(2*v_4));                          //s-1

    // SMC fluxes
    flu_M               = 1 - state_Mp - state_AM - state_AMp;
    flu_E_K_i           = ( R_gas * Temp ) / ( z_K  * Farad ) * unitcon * log( state_K_p / state_K_i );
    flu_h_r             = 0.1 * state_r;
    flu_v_cpl_i		    = - g_hat * ( state_v_i - state_v_j );
    flu_c_cpl_i         = - p_hat * ( state_ca_i - state_ca_j );
    flu_I_cpl_i         = - p_hatIP3 * ( state_ip3_i - state_ip3_j );
    flu_rho_i		    = 1;
    flu_ip3_i		    = Fmax_i *  pow(state_ip3_i,2) / ( pow(Kr_i,2) + pow(state_ip3_i,2) );
    flu_SRuptake_i      = B_i * pow(state_ca_i,2) / ( pow(state_ca_i,2) + pow(cb_i,2) );
    flu_CICR_i		    = C_i * pow(state_ca_sr_i,2) / ( pow(sc_i,2) + pow(state_ca_sr_i,2) ) *  ( pow(state_ca_i,4) ) / ( pow(cc_i,4) + pow(state_ca_i,4) );
    flu_extrusion_i	    = D_i * state_ca_i * (1 + ( state_v_i - vd_i ) / Rd_i );
    flu_leak_i 		    = L_i * state_ca_sr_i;
    flu_VOCC_i		    = G_Ca * ( state_v_i - v_Ca1 ) / ( 1 + exp( - ( state_v_i - v_Ca2 ) / ( R_Ca ) ) );
    flu_NaCa_i		    = G_NaCa * state_ca_i * ( state_v_i - v_NaCa ) / ( state_ca_i + c_NaCa ) ;
    flu_NaK_i		    = F_NaK;
    flu_Cl_i		    = G_Cl * (state_v_i - v_Cl);
    flu_K_i			    = G_K * state_w_i * ( state_v_i - vK_i );
    //flu_Kactivation_i   = pow((state_ca_i + c_w),2) / ( pow((state_ca_i + c_w),2) + bet*exp(-(state_v_i - v_Ca3)/R_K) );  // see NO pathway!
    flu_degrad_i	    = k_i * state_ip3_i;

    flu_J_stretch_i     = 0; //G_stretch/(1+exp(-alpha1*(P_str*state_r / flu_h_r - sig0))) * (state_v_i - Esac); // *** careful! flu_h_r is non-dimensionalised! + dfdt matrix!

    flu_v_KIR_i    = z_1 * state_K_p / unitcon + z_2;                                  // mV           state_K_p,
    flu_G_KIR_i    = exp( z_5 * state_v_i + z_3 * state_K_p / unitcon + z_4 );        // pS pF-1 =s-1  state_v_i, state_K_p
    flu_J_KIR_i    = F_il/gam * (flu_G_KIR_i) * (state_v_i-(flu_v_KIR_i));            // mV s-1 //     state_v_i, state_K_p
    


    // EC fluxes
    flu_v_cpl_j		= - g_hat * ( state_v_j - state_v_i );
    flu_c_cpl_j		= - p_hat * ( state_ca_j - state_ca_i );
    flu_I_cpl_j		= - p_hatIP3 * ( state_ip3_j - state_ip3_i );
    flu_rho_j 		= 1;
    flu_O_j 		= JO_j;
    flu_ip3_j		= Fmax_j * ( pow(state_ip3_j,2) ) / ( pow(Kr_j,2) + pow(state_ip3_j,2) );
    flu_ERuptake_j      = B_j * ( pow(state_ca_j,2) ) / ( pow(state_ca_j,2) + pow(cb_j,2) );
    flu_CICR_j			= C_j *  ( pow(state_ca_er_j,2) ) / ( pow(sc_j,2) + pow(state_ca_er_j,2) ) *  ( pow(state_ca_j,4) ) / ( pow(cc_j,4) + pow(state_ca_j,4) );
    flu_extrusion_j     = D_j * state_ca_j;
    flu_leak_j          = L_j * state_ca_er_j;

    flu_cation_j 	= G_cat * ( E_Ca - state_v_j) * 0.5 * ( 1 + tanh( ( log10( state_ca_j ) - m3cat )  /  m4cat  ) );

    flu_BKCa_j 		= 0.2 * ( 1 + tanh( ( (  log10(state_ca_j) - c) * ( state_v_j - b ) - a1 ) / ( m3b* pow(( state_v_j + a2 * ( log10( state_ca_j ) - c ) - b),2) + m4b ) ) );

    flu_SKCa_j 		= 0.3 * ( 1 + tanh( ( log10(state_ca_j) - m3s ) /  m4s ));
    flu_K_j 		= G_tot * ( state_v_j - vK_j ) * ( flu_BKCa_j + flu_SKCa_j ); // Reihenfolge!!
    flu_R_j 		= G_R * ( state_v_j - v_rest);
    flu_degrad_j 	= k_j * state_ip3_j;

    flu_J_stretch_j       = 0;// G_stretch / (1 + exp(-alpha1*(P_str * state_r / flu_h_r - sig0))) * (state_v_j - Esac);//careful!flu_h_r is non-dimensionalised! + dfdtmatrix!

// Mech fluxes
    flu_K1_c       = gam_cross * pow(state_ca_i,3);
    flu_K6_c       = flu_K1_c;

// NO pathway fluxes

   flu_P_NR2AO         = nvu_Glu(t,x,y)/(betA+nvu_Glu(t,x,y)); 
   flu_P_NR2BO         = nvu_Glu(t,x,y)/(betB+nvu_Glu(t,x,y));
   flu_openProbTerm    = 0.63 * flu_P_NR2AO + 11 * flu_P_NR2BO;
   flu_I_Ca            = (-4*v_n*G_M*P_Ca_P_M*(Ca_ex/M))/(1+exp(-0.08*(v_n+20)))*(exp(2*1e-3*v_n*Farad/(R_gas*Temp)))/(1-exp(2*1e-3*v_n*Farad/(R_gas*Temp)))*(0.63*flu_P_NR2AO+11*flu_P_NR2BO); 
   flu_phi_N           = 1 + Q1*state_ca_n + Q1*Q2*pow(state_ca_n,2) + Q1*Q2*Q3*pow(state_ca_n,3) + Q1*Q2*Q3*Q4*pow(state_ca_n,4);        // (102)
   flu_dphi_N          = Q1 + 2*Q1*Q2*state_ca_n + 3*Q1*Q2*Q3*pow(state_ca_n,2) + 4*Q1*Q2*Q3*Q4*pow(state_ca_n,3);            // == d(phi_N)/d(ind.Ca_n) ; (part of 101)
   flu_N               = (state_ca_n/flu_phi_N)*flu_dphi_N;                                                   // number of Ca2+ bound per calmodulin ; (101)
   flu_CaM             = state_ca_n/flu_N;                                      // concentration of calmodulin / calcium complexes ; (100)            

   flu_W_tau_w         = W_0*pow((tau_w + sqrt(16*pow(delt_wss,2)+pow(tau_w,2))-4*delt_wss),2) / (tau_w+sqrt(16*pow(delt_wss,2)+pow(tau_w,2))) ; 
   flu_F_tau_w         = (1/(1+alp*exp(-flu_W_tau_w)))-(1/(1+alp)); // -(1/(1+alp)) was added to get no NO at 0 wss (!)

   flu_k4             = C_4*pow(state_cGMP,m);
   flu_R_cGMP1        = (pow(state_cGMP,2))/(pow(state_cGMP,2)+pow(K_m_cGMP,2));
   flu_R_NO           = (state_NOi/(state_NOi+K_m_NO)) ;
   flu_v_Ca3          = -45*log10(state_ca_i-0.0838) + 223.276*flu_R_cGMP1 - 292.700*flu_R_NO - 198.55;
   flu_P_O            = pow((state_ca_i + c_wi ),2)/( pow((state_ca_i + c_wi ),2) + bet_i*exp(-(state_v_i - flu_v_Ca3) / (R_Kfit)) );
   flu_R_cGMP2        = (pow(state_cGMP,2))/(pow(state_cGMP,2)+pow(K_m_mlcp,2));
   flu_K2_c           = 16.88*k_mlcp_b+16.88*k_mlcp_c*flu_R_cGMP2;  // 17.64 / 16.75 errechnet sich aus dem Shift, um K2_c = 0.5 bei der baseline zu bekommen - muss vllt noch geaendert werden! 
   flu_K5_c           = flu_K2_c;
   flu_kmlcp          = k_mlcp_b + k_mlcp_c * flu_R_cGMP2;
   flu_Kactivation_i  = 0.075*(1+tanh((state_cGMP-9.7)))+ (pow((state_ca_i + c_w ),2) / ( pow((state_ca_i + c_w),2) + bet*exp(-(state_v_i - v_Ca3)/R_K) ));




// Differential Equations:
    du[i_radius]  = -w->a1 * e * (state_r / r0 - 1.) + w->a2 * state_r * pt; //Radius

    //AC:
    du[ R_k     ] = L_p * (flu_Na_k + flu_K_k + flu_Cl_k + flu_HCO3_k - flu_Na_s - flu_K_s - flu_Cl_s - flu_HCO3_s + X_k / state_R_k);  // m s-1
    du[ N_Na_k  ] = -flu_J_Na_k - 3 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_NBC_k;    // uMm s-1
    du[ N_K_k   ] = -flu_J_K_k + 2 * flu_J_NaK_k + flu_J_NKCC1_k + flu_J_KCC1_k -flu_J_BK_k; // uMm s-1
    du[ N_HCO3_k] = 2 * flu_J_NBC_k;                                                // uMm s-1
    du[ N_Cl_k  ] = du[ N_Na_k] + du[ N_K_k] - du[ N_HCO3_k];                       // uMm s-1, modified equation compared to the one of Ostby  //
    du[ N_Na_s  ] = - k_C * 0 - du[ N_Na_k];                                        // uMm s-1
    du[ N_K_s   ] = k_C * 0 - du[ N_K_k] ;                                          // uMm s-1
    du[ N_HCO3_s] = - du[ N_HCO3_k];                                                // uMm s-1
    du[ K_p     ] = flu_J_BK_k / (VR_pa * state_R_k) + flu_J_KIR_i / VR_ps;         // uM s-1
    du[ w_k     ] = flu_phi_w * (flu_w_inf - state_w_k);                            // s-1

    //SMC:
    du[ ca_i    ] = flu_c_cpl_i + flu_rho_i * ( flu_ip3_i - flu_SRuptake_i + flu_CICR_i - flu_extrusion_i + flu_leak_i - flu_VOCC_i + flu_NaCa_i + 0.1* flu_J_stretch_i);
    du[ ca_sr_i ] = flu_SRuptake_i - flu_CICR_i - flu_leak_i ;
    du[ v_i     ] = flu_v_cpl_i + gam * ( - flu_NaK_i - flu_Cl_i - 2 * flu_VOCC_i - flu_NaCa_i - flu_K_i - flu_J_stretch_i - flu_J_KIR_i );
    du[ w_i     ] = lam * (flu_Kactivation_i - state_w_i ) ;
    du[ ip3_i   ] = flu_I_cpl_i - flu_degrad_i ;          // **
    du[ K_i     ] = - flu_J_KIR_i - flu_K_i + flu_NaK_i;                                            // uM s-1


    //EC:
    du[ca_j     ] = flu_c_cpl_j + flu_rho_j * ( flu_ip3_j - flu_ERuptake_j + flu_CICR_j - flu_extrusion_j + flu_leak_j + flu_cation_j + flu_O_j + flu_J_stretch_j ) ;
    du[ca_er_j  ] = flu_ERuptake_j - flu_CICR_j - flu_leak_j ;
    du[v_j      ] = flu_v_cpl_j - 1/C_m * ( flu_K_j + flu_R_j ) ;
    du[ip3_j    ] = flu_I_cpl_j + J_PLC - flu_degrad_j ;  // **

    // Mech:
    du[ Mp   ] = K4_c * state_AMp + flu_K1_c * flu_M - (K2_c + K3_c) * state_Mp;
    du[ AMp  ] = K3_c * state_Mp + flu_K6_c * state_AM - (K4_c + K5_c) * state_AMp;
    du[ AM   ] = K5_c * state_AMp - ( K7_c + flu_K6_c ) * state_AM;

    // ***********NO pathway***********
    // NE:           
    du[ca_n]       = (((flu_I_Ca)/(2*Farad*v_spine) - (k_ex * (state_ca_n - Ca_rest)))/(1+lambda)) ;     //\muM
    du[nNOS]       = V_maxNOS * flu_CaM / (K_actNOS + flu_CaM) - mu2 * state_nNOS ;                  //\muM
    du[NOn]        = V_nNOS * LArg / ( pow((pow(LArg,2)+1),0.5) ) * (state_nNOS) - ((state_NOn-state_NOi)/tau_ni) - (k_O2* pow(state_NOn,2) * On); 
                         
    // SMC:              
    du[NOi]        = (state_NOn - state_NOi) / tau_ni + (state_NOj - state_NOi) / tau_ji - k_dno * state_NOi ; 
    du[E_b]        = -k1 * state_E_b * state_NOi + k_1 * state_E_6c + flu_k4 * state_E_5c ;     
    du[E_6c]       = k1 * state_E_b * state_NOi - k_1 * state_E_6c - k2 * state_E_6c - k3 * state_E_6c * state_NOi ;
    du[E_5c]       = k3 * state_E_6c * state_NOi + k2 * state_E_6c - flu_k4 * state_E_5c ;
    du[cGMP]       = V_max_sGC * state_E_5c - (k_pde * pow(state_cGMP,2) ) / (K_m_pde + state_cGMP );
       
    // EC:
    du[eNOS]       = ((K_dis * state_ca_j) / (K_eNOS + state_ca_j)) - mu2 * state_eNOS + g_max * flu_F_tau_w ;      // (104)
    du[NOj]        = V_eNOS*(LArg)/(pow((pow(LArg,2)+1),0.5))*(state_eNOS) - (state_NOj-state_NOi)/tau_ji - k_O2*pow((state_NOj),2)*Oj - state_NOj*4*3300/(pow(25,2));




    
   // printf("** %e** %e** %e** %e** %e** %e** %e** %e** %e** %e ** %e *** \n", L_p, flu_Na_k, flu_K_k, flu_Cl_k, flu_HCO3_k, flu_Na_s, flu_K_s, flu_Cl_s, flu_HCO3_s, X_k, state_R_k);

//flu_c_cpl_j, flu_rho_j, flu_ip3_j, flu_ERuptake_j, flu_CICR_j, flu_extrusion_j, flu_leak_j, flu_cation_j, flu_O_j, flu_J_stretch_j );
}

// Time-varying pressure at the root of the tree. 1 is nominal value. If
// you want to work in unscaled units, make sure you *multiply* by P0
// afterwards
double nvu_p0(double t) {
    //double p0 = 1. * 8000 / P0; // 8000 Pa   original: 1.5 * 8000 / P0;
    //double p0 = (0.5 * sin(t) + 1) * 8000 / P0; //
    double p0 = 1.5 * 8000 / P0;
return p0;
}

// Time-varying Glu input signal
double nvu_Glu(double t, double x, double y) {
    double Glu_min = 0;    
    double Glu_max = 1846; // uM - one vesicle (Santucci)
    double t_up   = 50;
    double t_down = 800;
    double blocks_activated_x = 4;
    double blocks_activated_x = 4;
    double Glu = ((0.5 + 0.5 *(tanh(100000*(x-0.0004*blocks_activated_x)+1))) *(0.5 + 0.5 *(tanh(100000*(y-0.0004*blocks_activated_y)+1))))        *          ((Glu_min + (Glu_max - Glu_min) / 2.0 * (1 + tanh(10*(t - t_up))) + (Glu_min - Glu_max) / 2.0 * (1 + tanh(10*(t - t_down)))));
    return Glu;
}


// Initial conditions. If you want spatial inhomegeneity, make it a
// function of the coordinates x and y. u0 is already allocated, you just
// need to fill in the entries
void nvu_ics(double *u0, double x, double y, nvu_workspace *w) {

    u0[i_radius]  = 1;   //0

    u0[R_k]       = 6.1e-8;                    //1
    u0[N_Na_k]    = 9.9796e-4;                 //2
    u0[N_K_k]     = 5.52782e-3;                //3
    u0[N_HCO3_k]  = 0.58804e-3;                //4
    u0[N_Cl_k]    = 0.32879e-3;                //5
    u0[N_Na_s]    = 4.301041e-3;               //6
    u0[N_K_s]     = 0.0807e-3;                 //7
    u0[N_HCO3_s]  = 0.432552e-3;               //8 
    u0[K_p]       = 3e3;                       //9
    u0[w_k]       = 0.1815e-3;                 //10

    u0[ca_i]      = 0.1;                       //11
    u0[ca_sr_i]   = 0.1;                       //12
    u0[v_i]       = -40;                       //13
    u0[w_i]       = 0.1;                       //14
    u0[ip3_i]     = 0.1;                       //15

    u0[K_i]       = 1e5;                       //16

    u0[ca_j]      = 1;                         //17
    u0[ca_er_j]   = 0.6;                       //18
    u0[v_j]       = -60;                       //19
    u0[ip3_j]     = 1.5;                       //20

    u0[Mp]        = 0.25;                      //21
    u0[AMp]       = 0.25;                      //22
    u0[AM]        = 0.25;                      //23
// NO pathway*****************
    u0[NOi]       = 0.07;                    //24
    u0[NOj]       = 0.07;                    //25
    u0[NOn]       = 0.2;                    //26
    u0[cGMP]      = 8.8;                    //27
    u0[eNOS]      = 3.8;                    //28
    u0[nNOS]      = 0.3;                    //29
    u0[ca_n]      = 0.1;                    //30
    u0[E_b]       = 0.3;                    //31  Eb+E6c+E5c=1 !
    u0[E_6c]      = 0.2;                    //32  
    u0[E_5c]      = 0.5;                    //33

}

//double Hill(double conc, double Hconst, double power) {
//    double hill = pow(conc,power) / ( pow(conc,power) + pow(Hconst,power) ) ;
//    return hill;
//}


