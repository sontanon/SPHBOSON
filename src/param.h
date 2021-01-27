#include <libconfig.h>

// String length.
#define MAX_STR_LEN 256

// Number of grid variables.
#define GNUM 3

#ifdef MAIN_FILE
/* CONFIG FILE */
config_t cfg;

/* GRID */
double dr      	= 0.0625;
MKL_INT NrInterior 	= 128;
MKL_INT dim        	= 132;
MKL_INT ghost      	= 2;
MKL_INT order 	= 4;

/* SCALAR FIELD PARAMETERS */
double m 	= 1.0;
double phi0 	= 0.2;
double w0 	= 0.7;
MKL_INT w_idx 	= 396;
MKL_INT fixedPhi = 1;

/* INITIAL DATA */
MKL_INT readInitialData 	= 0;
const char *log_alpha_i = NULL;
const char *log_psi_i 	= NULL;
const char *phi_i 	= NULL;
const char *w_i 	= NULL;
MKL_INT dimInitial	= 0;
MKL_INT ghost_i		= 2;
MKL_INT order_i 	= 4;
double dr_i		= 1.0;
double sigma		= 1.0;

/* SCALE INITIAL DATA */
double scale_u0 = 1.0;
double scale_u1 = 1.0;
double scale_u2 = 1.0;
double scale_u3 = 1.0;
double *u_seed = NULL;

/* NEXT SCALE ADVANCE */
double scale_next = 1.0;

/* SOLVER PARAMETERS */
MKL_INT solverType		= 1;
MKL_INT localSolver		= 1;
double epsilon		= 1E-5;
MKL_INT maxNewtonIter 	= 100;
double lambda0 		= 1.0E-3;
double lambdaMin 	= 1.0E-10;
MKL_INT useLowRank		= 1;

/* INITIAL GUESS CHECK */
MKL_INT max_initial_guess_checks = 8;
double norm_f0_target		 = 1.0E-05;

/* AUXILIARY ARRAYS FOR DERIVATIVES */
double *Dr_u	= NULL;
double *Drr_u	= NULL;

/* OUTPUT */
char	work_dirname[MAX_STR_LEN] = { 0 };
char initial_dirname[MAX_STR_LEN] = { 0 };
char   final_dirname[MAX_STR_LEN] = { 0 };

/* SWEEP CONTROL */
MKL_INT hwl_min = 10;
MKL_INT hwl_max = 100;
double w_max  = 1.0;
double w_min  = 0.0;
double w_step = 0.0;

/* ANALYSIS PARAMETERS */
double M_ADM;
double M_KOM;
double phi_max;
MKL_INT hwl_res;
#else
/* CONFIG FILE */
extern config_t cfg;

/* GRID */
extern double dr;
extern MKL_INT NrInterior;
extern MKL_INT dim;
extern MKL_INT ghost;
extern MKL_INT order;

/* SCALAR FIELD PARAMETERS */
extern double m;
extern double phi0;
extern double w0;
extern MKL_INT w_idx;
extern MKL_INT fixedPhi;

/* INITIAL DATA */
extern MKL_INT readInitialData;
extern const char *log_alpha_i;
extern const char *log_psi_i;
extern const char *phi_i;
extern const char *w_i;
extern MKL_INT dimInitial;
extern MKL_INT ghost_i;
extern MKL_INT order_i;
extern double dr_i;
extern double sigma;

/* SCALE INITIAL DATA */
extern double scale_u0;
extern double scale_u1;
extern double scale_u2;
extern double scale_u3;
extern double *u_seed;

/* NEXT SCALE ADVANCE */
extern double scale_next;

/* SOLVER PARAMETERS */
extern MKL_INT solverType;
extern MKL_INT localSolver;
extern double epsilon;
extern MKL_INT maxNewtonIter;
extern double lambda0;
extern double lambdaMin;
extern MKL_INT useLowRank;

/* INITIAL GUESS CHECK */
extern MKL_INT max_initial_guess_checks;
extern double norm_f0_target;

/* AUXILIARY ARRAYS FOR DERIVATIVES. */
extern double *Dr_u;
extern double *Drr_u;

/* OUTPUT */
extern char	work_dirname[MAX_STR_LEN];
extern char initial_dirname[MAX_STR_LEN];
extern char   final_dirname[MAX_STR_LEN];

/* SWEEP CONTROL */
extern MKL_INT hwl_min;
extern MKL_INT hwl_max;
extern double w_max;
extern double w_min;
extern double w_step;

/* ANALYSIS PARAMETERS */
extern double M_ADM;
extern double M_KOM;
extern double phi_max;
extern MKL_INT hwl_res;
#endif
