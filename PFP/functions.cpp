#include <stdio.h>
#include "tools.h"
#include "global.h"
#include "extern.h"
#include "functions.h"
//#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
//#include <mkl.h>

void inclusiveValue
(
 struct INCV &incv,
 double *deltaIn,
 struct PAR &par,
 struct VAR &var,
 double *nu
//=======================================================
)
{
  int j, k, n, s, t;
  double
  det,
  expInclusiveValue,            // exponentiated logit inclusive value
  max_w = -DBL_MAX,
  min_w = DBL_MAX;
  
  // Compute logit inclusive values w[t] for consumer i
  for (t=1; t<=n_period; t++){
    expInclusiveValue = 0;
    for (j=1; j<=n_product; j++){
      var.util[j] = deltaIn[idx[t][j]];
      for (k=1; k<=n_theta2; k++){
        var.util[j] += par.theta2[k] * nu[k] * var.X[idx[t][j]][rcID[k]];
      }
      var.num[t][j] = exp(var.util[j]);
      expInclusiveValue += var.num[t][j];
    }
    incv.w[t] = log(expInclusiveValue);
    max_w = fmax(max_w, incv.w[t]);
    min_w = fmin(min_w, incv.w[t]);
  }
  
  for (t=1; t<=n_period-1; t++){
    incv.y[t] = incv.w[t+1];
    incv.w_lag[t][1] = 1;
    incv.w_lag[t][2] = incv.w[t];
  }
  
  // Estimate lambda in transition process of inclusive values
  transpose(incv.w_lag, incv.w_lag_t, n_period-1, 2);
  matmul(incv.w_lag_t, incv.w_lag, incv.ww, 2, n_period-1, 2);
  matvecmul(incv.w_lag_t, incv.y, incv.wy, 2, n_period-1);
  
  det = fabs(incv.ww[1][1]*incv.ww[2][2] - incv.ww[1][2]*incv.ww[2][1]);
  incv.ww_inv[1][1] = incv.ww[2][2] / det;
  incv.ww_inv[2][2] = incv.ww[1][1] / det;
  incv.ww_inv[1][2] = incv.ww_inv[2][1] = -incv.ww[1][2] / det;
  matvecmul(incv.ww_inv, incv.wy, incv.lambda, 2, 2); // make sure lambda works
  matvecmul(incv.w_lag, incv.lambda, incv.wlambda, n_period-1, 2);
  
  // Estimate variance of psi
  incv.varPsi = 0;
  for (n=1; n<=n_period-1; n++){
    incv.varPsi += (incv.y[n] - incv.wlambda[n]) * (incv.y[n] - incv.wlambda[n]);
  }
  incv.varPsi /= n_period - 2;
  
  // Set up grid space
  max_w = max_w + 2.575829*sqrt(incv.varPsi);
  min_w = min_w - 2.575829*sqrt(incv.varPsi);
  for (s=1; s<=n_grid; s++){
    incv.wGrid[s] = min_w + (s-1) * (max_w - min_w)/(n_grid-1);
  }
}

void inclusiveValue
(
 //=======================================================
 struct INCV &incv,
 double *expDeltaIn,
 double **expMu,
 struct VAR &var
//=======================================================
)
{
  int j, n, s, t;
  double
  det,
  expInclusiveValue,            // exponentiated logit inclusive value
  max_w = -DBL_MAX,
  min_w = DBL_MAX;
  
  // Compute logit inclusive values w[t] for consumer i
  for (t=1; t<=n_period; t++){
    expInclusiveValue = 0;
    for (j=1; j<=n_product; j++){
      var.num[t][j] = expDeltaIn[idx[t][j]] * expMu[t][j];
      expInclusiveValue += var.num[t][j];
    }
    incv.w[t] = log(expInclusiveValue);
    max_w = fmax(max_w, incv.w[t]);
    min_w = fmin(min_w, incv.w[t]);
  }
  
  for (t=1; t<=n_period-1; t++){
    incv.y[t] = incv.w[t+1];
    incv.w_lag[t][1] = 1;
    incv.w_lag[t][2] = incv.w[t];
  }
  
  // Estimate lambda in transition process of inclusive values
  transpose(incv.w_lag, incv.w_lag_t, n_period-1, 2);
  matmul(incv.w_lag_t, incv.w_lag, incv.ww, 2, n_period-1, 2);
  matvecmul(incv.w_lag_t, incv.y, incv.wy, 2, n_period-1);
  
  det = fabs(incv.ww[1][1]*incv.ww[2][2] - incv.ww[1][2]*incv.ww[2][1]);
  incv.ww_inv[1][1] = incv.ww[2][2] / det;
  incv.ww_inv[2][2] = incv.ww[1][1] / det;
  incv.ww_inv[1][2] = incv.ww_inv[2][1] = -incv.ww[1][2] / det;
  matvecmul(incv.ww_inv, incv.wy, incv.lambda, 2, 2); // make sure lambda works
  matvecmul(incv.w_lag, incv.lambda, incv.wlambda, n_period-1, 2);
  
  // Estimate variance of psi
  incv.varPsi = 0;
  for (n=1; n<=n_period-1; n++){
    incv.varPsi += (incv.y[n] - incv.wlambda[n]) * (incv.y[n] - incv.wlambda[n]);
  }
  incv.varPsi /= n_period - 2;
  
  // Set up grid space
  max_w = max_w + 2.575829*sqrt(incv.varPsi);
  min_w = min_w - 2.575829*sqrt(incv.varPsi);
  for (s=1; s<=n_grid; s++){
    incv.wGrid[s] = min_w + (s-1) * (max_w - min_w)/(n_grid-1);
  }
  
}

void choicePath
(
 //=======================================================
 struct PRED &predOut,
 struct VAR  &var,
 struct DIAG &diag
//=======================================================
)
{
  int j, t;
  double den;
  
  // Initialize survivalRateard (var.survivalRate) at t=1 for all i
  predOut.survivalRate[1] = 1;
  
  // Compute predicted market share given deltaIn
  for (t=1; t<=n_period; t++){
    var.num[t][0] = exp(beta * var.value0[t]);
    den = 0;
    for (j=0; j<=n_product; j++){
      den += var.num[t][j];
    }
    
    // Add choice prob * (1-survivalRateard) to predicted market share
    if (isnan(den) || isinf(den)){
      //printf("NaN den\n");
      //printf("Num %d %.8f %.8f\n", 0, var.num[t][0], var.value0[t]);
      //printf("beta: %.8f expValfun: %.8f\n", beta, var.value0[t]);
      diag.exit = -1;
      return;
    }
    
    for (j=1; j<=n_product; j++){
      predOut.ms[t][j] += var.num[t][j] / den * predOut.survivalRate[t];
    }
    
    // Update the next-period survivalRateard
    if (t < n_period){
      predOut.survivalRate[t+1] = var.num[t][0] / den * predOut.survivalRate[t];
    }
    
  } // t loop
}

void initializeMu
(
 //=======================================================
 double ***expMuOut,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var
//=======================================================
)
{
  int i, j, k, t, id;
  double mu;
  
  for (i=1; i<=n_person; i++){
    for (t=1; t<=n_period; t++){
      for (j=1; j<=n_product; j++){
        id = idx[t][j];
        mu = 0;
        for (k=1; k<=n_theta2; k++)
          mu += par.theta2[k] * rnd.nu[i][k] * var.X[id][rcID[k]];
        expMuOut[i][t][j] = exp(mu);
      }
    }
  }
}

void initializeMui
(
 //=======================================================
 double **expMu,
 double *nu,
 struct PAR &par,
 struct VAR &var
//=======================================================
)
{
  int j, k, t, id;
  double mu;
  
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      id = idx[t][j];
      mu = 0;
      for (k=1; k<=n_theta2; k++)
        mu += par.theta2[k] * nu[k] * var.X[id][rcID[k]];
      expMu[t][j] = exp(mu);
    }
  }
}

void transProb
(
 //=======================================================
 struct VAR &var,
 struct INCV &incv
//=======================================================
)
{
  int s0, s1;
  double condMean, sumWeight;
  
  for (s0=1; s0<=n_grid; s0++){
    condMean = incv.lambda[1] + incv.lambda[2] * incv.wGrid[s0];
    sumWeight = 0;
    for (s1=1; s1<=n_grid; s1++){
      var.weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += var.weight[s1];
    }
    for (s1=1; s1<=n_grid; s1++){
      var.weightMat[s0][s1] = var.weight[s1]/sumWeight;
    }
  }
  
}

void expectedValuefun
(
 //=======================================================
 double *expectValuefun,
 double *valuefunIn,
 double *weight,
 struct INCV &incv
//=======================================================
)
{
  int s1, t;
  double
  condMean,     // conditional mean of next-period state given observed states
  value0,       // unweighted expected value function given observed states
  sumWeight;    // sum of transition prob weight
  
  for (t=1; t<=n_period; t++){
    // Compute conditional prob of state wGrid[s1] given state w_i[m][t]
    condMean = incv.lambda[1] + incv.lambda[2] * incv.w[t];
    sumWeight = 0;
    for (s1=1; s1<=n_grid; s1++){
      weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += weight[s1];
    }
    
    // Compute expected value function
    value0 = 0;
    for (s1=1; s1<=n_grid; s1++){
      value0 += weight[s1] * valuefunIn[s1];
    }
    expectValuefun[t] = value0 / sumWeight;
    
  } // t loop
}

void bellmanOperator
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR  &var,
 struct INCV &incv
//=======================================================
)
{
  int s0, s1;
  double value0;
  
  for (s0=1; s0<=n_grid; s0++){
    // value0: expected continuation value
    value0 = 0;
    for (s1=1; s1<=n_grid; s1++){
      value0 += valuefunIn[s1] * var.weightMat[s0][s1];
    }
    value0 *= beta;
    valuefunOut[s0] = log(1 + exp(incv.wGrid[s0] - value0)) + value0;
  } // end of s0 loop
  
}

void nfpBellman
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR  &var,
 struct DIAG &diag,
 struct INCV &incv
//=======================================================
)
{
  int s0, iter = 0;
  double dist = DBL_MAX;
  
  // Solve for exact value function
  while (dist > var.tol_bellman && iter <= 5000){
    bellmanOperator(valuefunOut, valuefunIn, var, incv);
    dist = 0;
    for (s0=1; s0<=n_grid; s0++){
      dist = fmax(dist, fabs(valuefunIn[s0] - valuefunOut[s0]));
      valuefunIn[s0] = valuefunOut[s0];
    }
    iter++;
    diag.bellmanIter++;
  }
  
}

int findNeighbor
(
 //=======================================================
 int    &n_neighbor,
 double **paramHistIn,
 struct PAR &par
//=======================================================
)
{
  int it, k, whichMin;
  double diffAbs, dist, min_dist;
  
  min_dist = DBL_MAX;
  whichMin = 1;
  
  // Find the nearest neighbor of par.theta2
  for (it=1; it<=n_neighbor; it++){
    dist = 0;
    for (k=1; k<=n_theta2; k++){
      diffAbs = par.theta2[k] - paramHistIn[it][k];
      dist += diffAbs*diffAbs;
    }
    if (dist < min_dist){
      whichMin = it;
      min_dist = dist;
    }
  } // end of (it) loop
  
  return whichMin;
}

void moveHistory
(
 //=======================================================
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &n_adjust
 //=======================================================
)
{
  int upperBound = 9 + imin(floor(pow(imh-n_adjust, power)), n_history-9);
  vmh++;
  if (vmh>n_neighbor){
    (vmh<=upperBound) ? (n_neighbor=vmh):(vmh=1);
  }
}

void moveHistory2
(
 int &vmh,
 int &n_neighbor
 )
{
  vmh++;
  if (vmh > n_neighbor){
    n_neighbor = imin(n_neighbor+1, n_history);
    vmh = 1;
  }
}

void computeMatrices
(
 //=======================================================
 double **invPHI,
 double **XZPHIZXZ,
 double **X,
 double **Z,
 double **Z_t
//=======================================================
)
{
  int i, j;
  double
  **X_t,
  **PHI,
  **XZ,
  **ZX,
  **XZPHI,
  **XZPHIZ,
  **XZPHIZX,
  **XZPHIZXinv;
  
  X_t         = dmatrix(1, n_theta1, 1, n_obs);
  PHI 				= dmatrix(1, n_inst, 1, n_inst);
  XZ 					= dmatrix(1, n_theta1, 1, n_inst);
  ZX 					= dmatrix(1, n_inst, 1, n_theta1);
  XZPHI 			= dmatrix(1, n_theta1, 1, n_inst);
  XZPHIZ 			= dmatrix(1, n_theta1, 1, n_obs);
  XZPHIZX 		= dmatrix(1, n_theta1, 1, n_theta1);
  XZPHIZXinv  = dmatrix(1, n_theta1, 1, n_theta1);
  
  gsl_matrix
  *gsl_invPHI      = gsl_matrix_alloc(n_inst, n_inst),
  *gsl_XZPHIZXinv  = gsl_matrix_alloc(n_theta1, n_theta1);
  
  // Compute 1st-stage GMM weight matrix
  transpose(Z, Z_t, n_obs, n_inst);
  matmul(Z_t, Z, PHI, n_inst, n_obs, n_inst);     				// PHI = Z'Z
  for (i=1; i<=n_inst; i++){
    for (j=1; j<=n_inst; j++){
      gsl_matrix_set (gsl_invPHI, i-1, j-1, PHI[i][j]);
    }
  }
  gsl_linalg_cholesky_decomp(gsl_invPHI);
  gsl_linalg_cholesky_invert(gsl_invPHI);
  for (i=1; i<=n_inst; i++){
    for (j=1; j<=n_inst; j++){
      invPHI[i][j] = gsl_matrix_get (gsl_invPHI, i-1, j-1);
    }
  }
  transpose(X, X_t, n_obs, n_theta1);
  matmul(X_t, Z, XZ, n_theta1, n_obs, n_inst);    				// XZ = X_t*Z = X'Z
  transpose(XZ, ZX, n_theta1, n_inst);                      	// ZX = (XZ)' = Z'X
  matmul(XZ, invPHI, XZPHI, n_theta1, n_inst, n_inst);   // XZPHI = (X'Z)inv(Z'Z)
  matmul(XZPHI, ZX, XZPHIZX, n_theta1, n_inst, n_theta1); // XZPHIZX = (X'Z)inv(Z'Z)(Z'X)
  for (i=1; i<=n_theta1; i++){
    for (j=1; j<=n_theta1; j++){
      gsl_matrix_set (gsl_XZPHIZXinv, i-1, j-1, XZPHIZX[i][j]);
    }
  }
  gsl_linalg_cholesky_decomp(gsl_XZPHIZXinv);
  gsl_linalg_cholesky_invert(gsl_XZPHIZXinv);
  for (i=1; i<=n_theta1; i++){
    for (j=1; j<=n_theta1; j++){
      XZPHIZXinv[i][j] = gsl_matrix_get (gsl_XZPHIZXinv, i-1, j-1);
    }
  }
  matmul(XZPHI, Z_t, XZPHIZ, n_theta1, n_inst, n_obs);    // XZPHIZ = (X'Z)inv(Z'Z)Z'
  // XZPHIZXZ = inv((X'Z)inv(Z'Z)(Z'X))(X'Z)inv(Z'Z)Z'
  matmul(XZPHIZXinv, XZPHIZ, XZPHIZXZ, n_theta1, n_theta1, n_obs);
  
  gsl_matrix_free (gsl_invPHI);
  gsl_matrix_free (gsl_XZPHIZXinv);
  
  free_dmatrix(X_t       	, 1, n_theta1, 1, n_obs);
  free_dmatrix(PHI 				, 1, n_inst, 1, n_inst);
  free_dmatrix(XZ 				, 1, n_theta1, 1, n_inst);
  free_dmatrix(ZX 				, 1, n_inst, 1, n_theta1);
  free_dmatrix(XZPHI 			, 1, n_theta1, 1, n_inst);
  free_dmatrix(XZPHIZ 		, 1, n_theta1, 1, n_obs);
  free_dmatrix(XZPHIZX 		, 1, n_theta1, 1, n_theta1);
  free_dmatrix(XZPHIZXinv , 1, n_theta1, 1, n_theta1);
  
}

void printMCMC
(
 //=======================================================
 int    &mc_id,
 int    &vmh,
 int    &imh,
 int    &seed,
 int    &accept,
 int    &neighborNew,
 int    &n_neighbor,
 double &llh,
 double &pllh,
 struct PAR &par_a,
 struct DIAG &diag,
 FILE *fileout
 //=======================================================
)
{
  int k;
  double llhAccept, llhReject;
  
  llhAccept = (accept==1) ? llh:pllh;
  llhReject = (accept==1) ? pllh:llh;
  diag.toc = gettimeofday_sec();
  fprintf(fileout, "%4d %7d %5d ", mc_id, imh+diag.nIterPrev,seed);
  for (k=1; k<=n_theta2; k++)
    fprintf(fileout, "%20.14f ", par_a.theta2[k]);
  for (k=1; k<=n_theta1; k++)
    fprintf(fileout, "%20.14f ", par_a.theta1[k]);
  
  fprintf(fileout,"%16.8f %16.8f %4d %4d %4d %4d %8.3f ",
          llhAccept,llhReject,diag.numReject,
          vmh,neighborNew,n_neighbor,float(diag.cumulAccept)/imh);
  fprintf(fileout,"%14.4f\n", diag.toc-diag.tic);
}

void mcmcDiag
(
 char   *filenameEstim,
 int    &imh,
 int    &mcID,
 int    &seed,
 struct DIAG &diag,
 struct POST &post,
 RInside &R
)
{
  int Jmax;
  
  printf("\nMCID: %2d, Seed: %2d. MCMC diagnostics at mcmc: %6d\n",mcID,seed,imh);
  R["mcmcDraws"] = createMatrix(post.theta, imh, n_theta);
  R["thin"] = n_thin;
  Rcpp::NumericVector M((SEXP) R.parseEval("h<- heidelberger(mcmcDraws,thin)"));
  
  if (M[1]>=0){ // Zero indexing for matrix M!!!
    diag.exit = 0;
    diag.n_burnin = (M[0]-1)*n_thin;
  }
  else if (M[1]==-300){
    diag.exit = -1;
    printf("Exiting MCMC simulation due to no movement.");
  }
  else{
    Jmax = fmin(heidel*diag.Jmax, n_mh);
    diag.Jmax = (Jmax%10>0) ? imin(Jmax+(10-Jmax%10),n_mh) : Jmax;
  }
}

void printOutput
(
 char   *filenameEstim,
 int    &imh,
 int    &mcID,
 int    &seed,
 struct DIAG &diag,
 struct POST &post
)
{
  int k;
  FILE * outp_stat;
  
  posteriorMean(imh, diag, post);
  outp_stat = fopen(filenameEstim, "a");
  fprintf(outp_stat, "%4d %4d %4d %6d ", mcID,seed,diag.exit,imh+diag.nIterPrev);
  for (k=1; k<=n_theta; k++)
    fprintf(outp_stat, "%12.6f %12.6f ", post.mean[k], post.stderr[k]);
  fprintf(outp_stat, "%12lld %14.4f\n", diag.bellmanIter, diag.toc-diag.tic);
  fclose(outp_stat);
  
}

Rcpp::NumericMatrix createMatrix
(
 double **param,
 const int m,
 const int n
)
{
  Rcpp::NumericMatrix M(m,n);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      M(i,j) = param[i+1][j+1];
    }
  }
  return(M);
}

Rcpp::NumericMatrix createMatrix
(
 double **param,
 const int m1,
 const int m2,
 const int n
)
{
  Rcpp::NumericMatrix M(m2-m1,n);
  for (int i=0; i<m2-m1; i++) {
    for (int j=0; j<n; j++) {
      M(i,j) = param[i+1][j+1];
    }
  }
  return(M);
}

void posteriorMean
(
 int    &imh,
 struct DIAG &diag,
 struct POST &post
)
{
  int i, k;
  
  if (diag.exit>=0){
    for (k=1; k<=n_theta; k++){
      post.mean[k] = 0;
      i = diag.n_burnin + 1;
      //for (i=diag.n_burnin+1; i<=imh; i++){
      while (i<=imh){
        post.mean[k] += post.theta[i][k];
        i += n_thin;
      }
      post.mean[k] /= (imh-diag.n_burnin)/n_thin;
      //post.mean[k] /= (double) (imh-diag.n_burnin);
      
      post.stderr[k] = 0;
      i = diag.n_burnin + 1;
      while (i<=imh){
        post.stderr[k] += (post.theta[i][k]-post.mean[k]) * (post.theta[i][k]-post.mean[k]);
        i += n_thin;
      }
      post.stderr[k] = sqrt(post.stderr[k] / ((imh-diag.n_burnin)/n_thin-1));
      //post.stderr[k] = sqrt(post.stderr[k] / (imh - diag.n_burnin - 1));
    }
  }
  else if (diag.exit<0){
    for (k=1; k<=n_theta; k++){
      post.mean[k] = 0;
      post.stderr[k] = 0;
    }
  }
}

void readData
(
 int   &mcID,
 FILE* inpData,
 struct VAR  &var,
 struct DATA &data
)
{
  int j, k, t, id, imc, itime, iproduct, f=0;
  
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      id = idx[t][j];
      f = fscanf(inpData, "%d %d %d", &imc, &itime, &iproduct);
      f = fscanf(inpData, "%lf %lf %lf ", &data.ms[t][j], &data.ms[t][0], &data.xi[t][j]);
      for (k=1; k<=n_theta1; k++)
        f = fscanf(inpData, "%lf ", &var.X[id][k]);
      for (k=1; k<=n_inst; k++)
        f = fscanf(inpData, "%lf ", &var.Z[id][k]);
      if (imc > mc_last || itime > n_period || iproduct > n_product || f==EOF){
        error("Dimensions of data and estimation model do not match.");
      }
    }
  }
}

void readRand
(
 int   &mcID,
 FILE* inpRand,
 struct RND rnd
)
{
  int i, k, f=0;
  
  for (i=1; i<=n_person; i++){
    for (k=1; k<=n_theta2; k++)
      f = fscanf(inpRand, "%lf ", &rnd.nu[i][k]);
  }
  if (f==EOF)
    error("EOF error in readRand.");
}

void readSeed
(
 double ***paramSeed,
 char* filenameSeed
 )
{
  int i, k, s, f=0;
  FILE* inpSeed;
  
  inpSeed = fopen(filenameSeed,"r");
  for (i=1; i<=n_mc; i++){
    for (s=1; s<=n_seed; s++){
      for (k=1; k<=n_theta2; k++){
        f = fscanf(inpSeed, "%lf ", &paramSeed[i][s][k]);
      }
    }
  }
  if (f==EOF)
    error("EOF error in readSeed.");
  fclose(inpSeed);
}

void initialize
(
 int &mcID,
 int &seed,
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &neighborPast,
 int &ptAccpt,
 struct FXP &fxp,
 struct PAR &par_a,
 struct PAR &par_c,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct HIST &hist,
 gsl_rng *rng
)
{
  int i, j, k, t, s0;
  
  computeMatrices(var.invPHI, var.XZPHIZXZ, var.X, var.Z, var.Z_t);
  
  for (k=1; k<=n_theta2; k++){
    hist.param[1][k] = par_c.theta2[k] = par_a.theta2[k] = hist.paramSeed[mcID][seed][k];
  }
  imh = 1;
  vmh = 1;
  n_neighbor = 10;     // max until mcmc iteration reaches n_history
  neighborPast = 1;
  
  // Reset history
  resetHistory(2, hist);
  
  // Initialize BLP & Bellman fixed points
  fxp.expDeltaIn = hist.expDelta[1]; // Assign pointers
  fxp.valuefunIn = hist.valuefun[1]; // Assign pointers
  for (i=1; i<=n_person; i++){
    for (s0=1; s0<=n_grid; s0++){
      fxp.valuefunIn[i][s0] = 0;
    }
  }
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      // fxp.expDeltaIn[idx[t][j]] = exp(data.delta[idx[t][j]]);
      fxp.expDeltaIn[idx[t][j]] = data.ms[t][j] / data.ms[t][0];
    }
  }
  
  diag.lipschitz = 0;
  diag.lipschitzIter = 1;
  diag.blpIter = 0;
  diag.bellmanIter = 0;
  diag.cumulAccept = 0;
  hist.reset = 0;
  hist.accept[vmh] = 1;
  diag.maxReject = 10;
  diag.numReject = 0;
  diag.numRejectTotal = 0;
  ptAccpt = 0;
  diag.exit = 1;
  diag.Jmax = 0.1*n_mh;
  diag.nIterPrev = 0;
  
  gsl_rng_set(rng, n_obs*1000+n_person*n_grid*100+mcID*10+seed);
}

void reinitialize
(
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &neighborPast,
 int &ptAccpt,
 struct FXP &fxp,
 struct PAR &par_a,
 struct DIAG &diag,
 struct HIST &hist
)
{
  int k;
  
  for (k=1; k<=n_theta2; k++){
    hist.param[1][k] = par_a.theta2[k];
  }
  copyVector(hist.expDelta[neighborPast], hist.expDelta[1], n_obs);
  copyMatrix(hist.valuefun[neighborPast], hist.valuefun[1], n_person,n_grid);
  fxp.expDeltaIn = hist.expDelta[1]; // Assign pointers
  fxp.valuefunIn = hist.valuefun[1]; // Assign pointers

  diag.nIterPrev = imh;
  imh = 0;
  vmh = 1;
  n_neighbor = 10;     // max until mcmc iteration reaches n_history
  neighborPast = 1;
  
  // Reset history
  resetHistory(2, hist);
  
  diag.lipschitz = 0;
  diag.lipschitzIter = 1;
  diag.blpIter = 0;
  diag.bellmanIter = 0;
  diag.cumulAccept = 0;
  hist.reset = 0;
  hist.accept[vmh] = 1;
  diag.maxReject = 10;
  diag.numReject = 0;
  diag.numRejectTotal = 0;
  ptAccpt = 0;
  diag.exit = 1;
  diag.Jmax = 0.1*n_mh;
  
}

void resetHistory
(
 int indx,
 struct HIST &hist
)
{
  int i, k, n, it, s0;
  for (it=indx; it<=n_history; it++){
    for (k=1; k<=n_theta2; k++)
      hist.param[it][k] = DBL_MAX;
    for (n=1; n<=n_obs; n++)
      hist.expDelta[it][n] = 0;
    for (i=1; i<=n_person; i++){
      for (s0=1; s0<=n_grid; s0++)
        hist.valuefun[it][i][s0] = 0;
    }
  }
}

void resetHistParam
(
 int idFirst,
 int idLast,
 struct HIST &hist
)
{
  int k, it;
  for (it=idFirst; it<=idLast; it++){
    for (k=1; k<=n_theta2; k++)
      hist.param[it][k] = DBL_MAX;
  }
}

void memoryAllocate
(
 double **jump,
 struct FXP &fxp,
 struct PAR &par_c,
 struct PAR &par_a,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct HIST &hist,
 struct INCV &incv,
 struct POST &post,
 struct PRED &pred
)
{
  idx            = imatrix(1, n_period, 1, n_product);
  
  *jump          = dvector(1, n_theta2);
  hist.accept    = ivector(1, n_history);
  hist.paramIn   = dvector(1, n_param);
  hist.paramSeed = darray3d(1, n_mc, 1, n_seed, 1, n_theta2);
  hist.param     = dmatrix(1, n_history, 1, n_theta2);
  hist.expDelta  = dmatrix(1, n_history, 1, n_obs);
  hist.valuefun  = darray3d(1, n_history, 1, n_person, 1, n_grid);
  
  fxp.expDeltaOut = dvector(1, n_obs);
  fxp.valuefunOut = dmatrix(1, n_person, 1, n_grid);
  
  par_c.theta1  = dvector(1, n_theta1);
  par_c.theta2  = dvector(1, n_theta2);
  par_a.theta1  = dvector(1, n_theta1);
  par_a.theta2  = dvector(1, n_theta2);
  
  rnd.nu        = dmatrix(1, n_person, 1, n_theta2);
  
  var.util      = dvector(1, n_product);
  var.num       = dmatrix(1, n_period, 0, n_product);
  var.expMu0    = darray3d(1, n_person, 1, n_period, 1, n_product);
  var.expMu1    = darray3d(1, n_person, 1, n_period, 1, n_product);
  var.value0    = dvector (1, n_period);
  var.X         = dmatrix(1, n_obs, 1, n_theta1);
  var.Z         = dmatrix(1, n_obs, 1, n_inst);
  var.Z_t       = dmatrix(1, n_inst, 1, n_obs);
  var.invPHI    = dmatrix(1, n_inst, 1, n_inst);
  var.XZPHIZXZ  = dmatrix(1, n_theta1, 1, n_obs);
  var.Xtheta1   = dvector(1, n_obs);
  var.delta     = dvector(1, n_obs);
  var.e         = dvector(1, n_obs);
  var.Ze        = dvector(1, n_inst);
  var.PHIZe     = dvector(1, n_inst);
  var.weight    = dvector(1, n_grid);
  var.weightMat = dmatrix(1, n_grid, 1, n_grid);
  
  data.ms       = dmatrix(1, n_period, 0, n_product);
  data.xi       = dmatrix(1, n_period, 1, n_product);
  data.delta    = dvector(1, n_obs);
  
  incv.w        = dvector(1, n_period),
  incv.wGrid    = dvector(1, n_grid),
  incv.lambda   = dvector(1, 2);
  incv.y        = dvector(1, n_period - 1);
  incv.wy       = dvector(1, 2);
  incv.wlambda  = dvector(1, n_period - 1);
  incv.w_lag    = dmatrix(1, n_period - 1, 1, 2);
  incv.w_lag_t  = dmatrix(1, 2, 1, n_period - 1);
  incv.ww       = dmatrix(1, 2, 1, 2);
  incv.ww_inv   = dmatrix(1, 2, 1, 2);
  
  post.mean     = dvector(1, n_theta);
  post.stderr   = dvector(1, n_theta);
  post.theta    = dmatrix(1, n_mh, 1, n_theta);
  
  pred.ms            = dmatrix(1, n_period, 1, n_product);
  pred.survivalRate  = dvector(1, n_period);
}

void memoryDeallocate
(
 double **jump,
 struct FXP &fxp,
 struct PAR &par_c,
 struct PAR &par_a,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct HIST &hist,
 struct INCV &incv,
 struct POST &post,
 struct PRED &pred
)
{
  free_imatrix(idx,           1, n_period, 1, n_product);
  free_dvector(*jump,         1, n_theta2);
  
  free_dvector(fxp.expDeltaOut, 1, n_obs);
  free_dmatrix(fxp.valuefunOut, 1, n_person, 1, n_grid);
  
  free_dvector(par_c.theta1,  1, n_theta1);
  free_dvector(par_c.theta2,  1, n_theta2);
  free_dvector(par_a.theta1,  1, n_theta1);
  free_dvector(par_a.theta2,  1, n_theta2);
  
  free_dmatrix(rnd.nu,  1, n_person, 1, n_theta2);
  
  free_dvector(var.util,      1, n_product);
  free_dmatrix(var.num,       1, n_period, 0, n_product);
  free_darray3d(var.expMu0,   1, n_person, 1, n_period, 1, n_product);
  free_darray3d(var.expMu1,   1, n_person, 1, n_period, 1, n_product);
  free_dvector(var.value0,    1, n_period);
  free_dmatrix(var.X,         1, n_obs, 1, n_theta1);
  free_dmatrix(var.Z,         1, n_obs, 1, n_inst);
  free_dmatrix(var.Z_t,       1, n_inst, 1, n_obs);
  free_dmatrix(var.invPHI,    1, n_inst, 1, n_inst);
  free_dmatrix(var.XZPHIZXZ,  1, n_theta1, 1, n_obs);
  free_dvector(var.Xtheta1,   1, n_obs);
  free_dvector(var.delta,     1, n_obs);
  free_dvector(var.e,         1, n_obs);
  free_dvector(var.Ze,        1, n_inst);
  free_dvector(var.PHIZe,     1, n_inst);
  free_dvector(var.weight,    1, n_grid);
  free_dmatrix(var.weightMat, 1, n_grid, 1, n_grid);
  
  free_dmatrix(data.ms,       1, n_period, 0, n_product);
  free_dmatrix(data.xi,       1, n_period, 1, n_product);
  free_dvector(data.delta,    1, n_obs);
  
  free_ivector(hist.accept,     1, n_history);
  free_dvector(hist.paramIn,    1, n_param);
  free_darray3d(hist.paramSeed, 1, n_mc, 1, n_seed, 1, n_theta2);
  free_dmatrix(hist.param,      1, n_history, 1, n_theta2);
  free_dmatrix(hist.expDelta,   1, n_history, 1, n_obs);
  free_darray3d(hist.valuefun,  1, n_history, 1, n_person, 1, n_grid);
  
  free_dvector(incv.w,        1, n_period);
  free_dvector(incv.wGrid,    1, n_grid);
  free_dvector(incv.lambda,   1, 2);
  free_dvector(incv.y,        1, n_period - 1);
  free_dvector(incv.wy,       1, 2);
  free_dvector(incv.wlambda,  1, n_period - 1);
  free_dmatrix(incv.w_lag,    1, n_period - 1, 1, 2);
  free_dmatrix(incv.w_lag_t,  1, 2, 1, n_period - 1);
  free_dmatrix(incv.ww,       1, 2, 1, 2);
  free_dmatrix(incv.ww_inv,   1, 2, 1, 2);
  
  free_dvector(post.mean,    1, n_theta);
  free_dvector(post.stderr,  1, n_theta);
  free_dmatrix(post.theta,   1, n_mh, 1, n_theta);
  
  free_dmatrix(pred.ms,           1, n_period, 1, n_product);
  free_dvector(pred.survivalRate, 1, n_period);
  
}
