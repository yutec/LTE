struct DATA{ // market share, prices, unobserved shocks
  double
  **ms,           // [n_period][n_product+1]: market share
  **xi,           // [n_period][n_product]: unobserved demand shock
  *delta;         // [n_obs]: mean utility vector
};

struct PAR{ // demand and price-equation parameters
  double
  *theta1,		// [n_theta1]: linear parameters
  *theta2;    // [n_theta2]: nonlinear parameters
};

struct HIST{
  int
  reset,
  *accept;       // history of acceptance
  
  double
  *paramIn,          // input parameter vector
  ***paramSeed,
  **param,
  **expDelta,
  ***valuefun;
};

struct INCV{
  double
  varPsi,   // variance of psi (error in autoregression of inclusive values)
  *lambda,  // [2]: coefficient vector of autoregression of inclusive values
  *wGrid,   // [n_grid]: discretized grid points
  *w,       // [n_obs]: inclusive values
  
  *y, *wy, *wlambda,  // supplemental variables to compute least squares estiamtes
  **w_lag, **w_lag_t, **ww, **ww_inv; // of lambda & varPsi

};

struct RND{ // random simulation of consumer draws
	double
  **nu;  // [n_person]: simulation draws for unobserved consumer heterogeneity
};

struct VAR{ // collection of variables commonly used

	double
  tol_blp,
  tol_bellman,
  *util,
  **num,
  ***expMu0,  // [n_person][n_period][n_product]: theta2*nu*X
  ***expMu1,  // [n_person][n_period][n_product]: theta2*nu*X
  ***expMuIn, // [n_person][n_period][n_product]: theta2*nu*X
  ***expMuOut,
  *mu_mkl,
  *expMu_mkl,

  *value0,    // [n_period]: option value of no purchase option

  **X,        // [n_obs][2]: covariates
  **Z,        // [n_obs][n_inst]: instruments
  **Z_t,      // [n_inst][n_obs]: tranposed Z
  **invPHI,   // [n_inst][n_inst]: inverse of PHI = Z_t' Z
  **XZPHIZXZ,
  *Xtheta1,
  *delta,
  *e,         // [n_obs]: residual of least squares regression of delta wrt X
  *Ze,
  *PHIZe,
  
  *weight,    // [n_grid]: transition probability for value function computation
  **weightMat;
};

struct FXP{
  double
  *expDeltaIn,
  *expDeltaOut,
  **valuefunIn,
  **valuefunOut;
};

struct PRED{
  double
  **ms,           // [n_period][n_product]: predicted market share
  *survivalRate;  // [n_period]: participation rate of consumers
};

struct POST{  // posterior means & std. deviations
  double
  *mean, // posterior means of theta1 & theta2
  *stderr,     // posterior standard deviations
  **theta;
};

struct DIAG{  // diagnostic variables
  int
  exit,
  Jmax,
  nIterPrev,
  n_burnin,
  maxReject,
  numReject,
  numRejectTotal,
  cumulAccept;
  
  unsigned long long int
  lipschitzIter,
  blpIter,
  bellmanIter;

  double
  tic, toc,
  lipschitz;
};



