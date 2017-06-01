#ifndef TOOLS_H
#define TOOLS_H
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>

//#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
//#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

double gettimeofday_sec(void);
void error(const char error_text[]);

//========================================================================
// Linear Algebra
//========================================================================
void matvecmul(double **A, double *B, double *AB, long m, long n);
void matvecmul(double **A, double *B, adouble *AB, long m, long n);
void matvecmul(double **A, adouble *B, adouble *AB, long m, long n);
void matvecmul(adouble **A, adouble *B, adouble *AB, long m, long n);
void matmul(double **A, double **B, double **AB, long m, long n, long l);
void matmul(adouble **A, adouble **B, adouble **AB, long m, long n, long l);
void transpose(double **A, double **A_t, int m, int n);
void transpose(adouble **A, adouble **A_t, int m, int n);

void copyVector(double *source, double *target, int n);
void copyVector(const double *source, double *target, int n);
void copyMatrix(double **source, double **target, int m, int n);
void copy3darray(double ***source, double ***target, int m, int n, int l);
void logVector(double *source, double *target, int size);

/*========================================================================
  Memory Allocation for Arrays
========================================================================*/
int *ivector(long n);
unsigned int *uivector(long n);
int **imatrix(long nrow, long ncol);
int ***iarray3d(long n1, long n2, long n3);

double *dvector(long n);
double **dmatrix(long nrow, long ncol);
double ***darray3d(long n1, long n2, long n3);

adouble *advector(long n);
adouble **admatrix(long nrow, long ncol);
adouble ***adarray3d(long n1, long n2, long n3);

/*========================================================================
 Memory De-allocation for Arrays
 ========================================================================*/
void free_ivector(int *v);
void free_uivector(unsigned int *v);
void free_imatrix(int **m);
void free_iarray3d(int ***m);

void free_dvector(double *v);
void free_dmatrix(double **m);
void free_darray3d(double ***m);

void free_advector(adouble *v);
void free_admatrix(adouble **m);
void free_adarray3d(adouble ***m);

/*========================================================================
 Global variables
 ========================================================================*/

struct DAT{
  double
  **ms,           // [n_period][n_product+1]: market share
  *ms0,
  **xi,           // [n_period][n_product]: unobserved demand shock
  *delta,         // [n_obs]: mean utility vector
  **valfun;
};

struct PAR{
  double
  *theta1,		// [n_theta1]: linear parameters
  *theta2,    // [n_theta2]: nonlinear parameters
  *sd1,
  *sd2,
  *xi,
  **valfun,
  *g;
};

struct PAR_AD{
  adouble
  *theta1,		// [n_theta1]: linear parameters
  *theta2,   // [n_theta2]: nonlinear parameters
  *xi,
  **valfun,
  *g;
};

struct INCV{
  double
  varPsi,   // variance of psi (error in autoregression of inclusive values)
  *w,       // [n_obs]: inclusive values
  *wGrid,   // [n_grid]: discretized grid points
  *lambda,  // [2]: coefficient vector of autoregression of inclusive values
  *y, *wy, *wlambda,  // supplemental variables to compute least squares estiamtes
  **w_lag, **w_lag_t, **ww, **ww_inv; // of lambda & var_psi
};

struct INCV_AD{
  adouble
  varPsi,
  *w,
  *wGrid,
  *lambda,
  *y, *wy, *wlambda,
  **w_lag, **w_lag_t, **ww, **ww_inv;
};

struct RND{
  double
  **nu;  // [n_person][n_theta2]: simulation draws for unobserved consumer heterogeneity
};

struct VAR{
  double
  **X,        // [n_obs][2]: covariates
  **Z,        // [n_obs][n_inst]: instruments
  **Z_t,      // [n_inst][n_obs]: tranposed Z
  **invPHI,   // [n_inst][n_inst]: inverse of PHI = Z_t' Z
  **XZPHIZXZ,
  tol_blp,
  tol_bellman,
  
  **num,
  *value0,        // [n_period]: option value of no purchase option
  *Xtheta1,
  *e,         // [n_obs]: residual of least squares regression of delta wrt X
  *Ze,
  *PHIZe,
  **weightMat;
};

struct VAR_AD{
  adouble
  **num,
  *value0,        // [n_period]: option value of no purchase option
  *Xtheta1,
  *e,         // [n_obs]: residual of least squares regression of delta wrt X
  *Ze,
  *PHIZe,
  **weightMat;
  
};

struct FXP{
  double
  *delta,
  *expDeltaIn,           // [n_obs]: input delta
  *expDeltaOut,
  **valuefunIn,       // [n_person][n_grid]: input value function
  **valuefunOut,
  ***expMu;
};

struct FXP_AD{
  adouble
  *expDeltaIn,
  **expMu,
  **valuefunOut;
};

struct PRED{
  double
  **ms,
  *survivalRate;
};

struct PRED_AD{
  adouble
  **ms,           // [n_period][n_product]: predicted market share
  *survivalRate;  // [n_period]: participation rate of consumers
};

struct POST{
  double
  **theta,    // posterior means of theta1 & theta2
  **stddev,   // posterior standard deviations
  *obj;
};

struct DIAG{
  unsigned long long int
  lipschitzIter,
  blpIter,
  bellmanIter;
  
  int
  FC_evals,
  GA_evals,
  H_evals,
  HV_evals,
  iters,
  major_iters,
  minor_iters;
  
  double
  lipschitz;
};

struct USERSTRUCT{
  short int
  tagJacobian,
  tagHessian;
  
  int
  mcID,
  seed,
  nnzHconst,
  *jacIndexConsInt,
  *jacIndexVarsInt,
  *hessIndexRows,
  *hessIndexCols,
  optionsJac[4];
  
  unsigned int
  *jacIndexCons,
  *jacIndexVars;
  
  double
  tic,
  *xcopy,
  *jac,
  **hessConstraint;
  
  adouble
  constraint,
  *c_ad;
  
  struct DAT dat;     // store dat
  struct FXP fxp;     // store fixed points
  struct PAR par;     // store accepted parameters
  struct RND rnd;     // store random numbers
  struct VAR var;     // store frequently used variables
  struct DIAG diag;   // store diagnostics including Lipschitz constat, BLP/Bellman iterations
  struct INCV incv;   // store inclusive value related variables
  struct PRED pred;
  
  struct PAR_AD par_ad;
  struct FXP_AD fxp_ad;
  struct VAR_AD var_ad;
  struct INCV_AD incv_ad;
  struct PRED_AD pred_ad;
};

#endif