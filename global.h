struct DAT{
  // market share, prices, unobserved shocks
  double
  **ms,           // [n_period][n_product+1]: market share
  **xi,           // [n_period][n_product]: unobserved demand shock
  *delta;         // [n_obs]: mean utility vector
};

struct PAR{
	// demand and price-equation parameters
  double
  *theta1,		// [n_theta1]: linear parameters
  *theta2,    // [n_theta2]: nonlinear parameters
  *sd1,
  *sd2;
};

struct INCV{
  double
  varPsi,  // variance of psi (error in autoregression of inclusive values)
  *lambda,  // [2]: coefficient vector of autoregression of inclusive values
  *wGrid,  // [n_grid]: discretized grid points
  *w,       // [n_obs]: inclusive values
  
  *y, *wy, *wlambda,  // supplemental variables to compute least squares estiamtes
  **w_lag, **w_lag_t, **ww, **ww_inv; // of lambda & var_psi

};

struct RND{
	double
  **nu;  // [n_person][n_theta2]: simulation draws for unobserved consumer heterogeneity
};

struct VAR{
	double
  **num,
  *value0,        // [n_period]: option value of no purchase option

  **X,        // [n_obs][2]: covariates
  **Z,        // [n_obs][n_inst]: instruments
  **Z_t,      // [n_inst][n_obs]: tranposed Z
  **invPHI,   // [n_inst][n_inst]: inverse of PHI = Z_t' Z
  **XZPHIZXZ,
  *Xtheta1,
  *e,         // [n_obs]: residual of least squares regression of delta wrt X
  *Ze,
  *PHIZe,
  
  *weight, // [n_grid]: transition probability weights for value function computation
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

struct PRED{
  double
  **ms,           // [n_period][n_product]: predicted market share
  *survivalRate;  // [n_period]: participation rate of consumers
};

struct DIAG{
  unsigned long long int
  lipschitzIter,
  blpIter,
  bellmanIter;
  
  int
  exit,
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
  int
  mcID,
  seed;
  
  double
  tic;
  
  struct DAT dat;     // store dat
  struct FXP fxp;     // store fixed points
  struct PAR par;     // store accepted parameters
  struct RND rnd;     // store random numbers
  struct VAR var;     // store frequently used variables
  struct DIAG diag;   // store diagnostics including Lipschitz constat, BLP/Bellman iterations
  struct INCV incv;   // store inclusive value related variables
  struct PRED pred;
};