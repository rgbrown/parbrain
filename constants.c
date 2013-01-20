// Constants for brain simulation. Declared extern in brain.h 

const int NDEFAULT      = 9;
const int NSUBDEFAULT   = 3;
const int NSYMBOLS      = 3;

const double RMIN       = 10e-6 ;  // m
const double R0         = 10e-6 ;  // m (for nondimensionalising)
const double L0         = 200e-6;  // m (for nondimensionalising)
const double LRR        = 20    ;  // Nondimensional (length to radius ratio)
const double HRR        = 0.1   ;  // Nondimensional (thickness to radius ratio)
const double PROOT      = 8000  ;  // Pa (root pressure)
const double PCAP       = 4000  ;  // Pa (capillary bed pressure)
const double MU         = 3.5e-3;  // Pa s (Viscosity)
const double RSCALE     = 0.6   ;  // Dimensionless
const double P0         = 8000  ;  // Pa (scaling factor for nondim)
const double E0         = 66e3  ;  // Pa
const double EPASSIVE   = 66e3  ;  // Pa
const double EACTIVE    = 233e3 ;  // Pa
const double ETA        = 2.8e2 ;  // Pa s
const double GAMMA      = 5     ;  // nondim
const double CSTAR      = 3.5   ;  // nondim
const double TAUM       = 5     ;  // s
const double PW         = 5e-5  ;  // m s^-1
const double AVR        = 5e3   ;  // m^2 m^-3
const double RCAP       = 5e-6  ;  // m
const double CIN        = 2.2e-2;  // mol m^-3
const double C0         = 2.2e-2;
const double MNORMAL    = 4e-3  ;  // mol m^-3 s^-1
const double M0         = 4e-3  ;
const double T0         = 1     ;  // s
