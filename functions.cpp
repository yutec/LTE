#include <stdio.h>
#include "tools.h"
#include "global.h"
#include "extern.h"
#include "functions.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

void printOutput
(
 char   *filenameEstim,
 int    &mcID,
 int    &seed,
 int    &nStatus,
 double &obj,
 double &toc,
 struct USERSTRUCT &us
)
{
  int k;
  double modulus = us.diag.lipschitz/us.diag.lipschitzIter;
  FILE * outp_stat;
  
  outp_stat = fopen(filenameEstim, "a");
  fprintf(outp_stat, "%6d %6d %6d ", mcID,seed,nStatus);
  for (k=1; k<=n_theta2; k++)
    fprintf(outp_stat, "%12.6f %12.6f ", us.par.theta2[k], us.par.sd2[k]);
  for (k=1; k<=n_theta1; k++)
    fprintf(outp_stat, "%12.6f %12.6f ", us.par.theta1[k], us.par.sd1[k]);
  fprintf(outp_stat, "%14.5f %8d %8d %8lld %12lld %8.4f %14.4f\n", obj,us.diag.iters,
          us.diag.FC_evals,us.diag.blpIter,us.diag.bellmanIter,modulus,toc-us.tic);
  fclose(outp_stat);
  printf("MC experiment %d seed %d complete. Total time elapsed: %f\n",
         mcID, seed, toc-us.tic);

}

void initializeMu
(
 //=======================================================
 double ***expMu,
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
        expMu[i][t][j] = exp(mu);
      }
    }
  }
}

void initialize
(
 int mcID,
 int seed,
 double *xInitial,
 double *paramSeed,
 struct USERSTRUCT &us
)
{
  int i, j, k, t, s0;
  
  us.mcID = mcID;
  us.seed = seed;
  
  for (k=1; k<=n_theta2; k++){
    xInitial[k-1] = paramSeed[k];
  }

  computeMatrices(us.var.invPHI,us.var.XZPHIZXZ,us.var.X,us.var.Z,us.var.Z_t);
  
  // Initialize BLP & Bellman fixed points
  for (i=1; i<=n_person; i++){
    for (s0=1; s0<=n_grid; s0++){
      us.fxp.valuefunIn[i][s0] = 0;
    }
  }
  
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      us.fxp.expDeltaIn[idx[t][j]] = us.dat.ms[t][j] / us.dat.ms[t][0];
    }
  }
  us.diag.exit = 0;
  us.diag.lipschitz = 0;
  us.diag.lipschitzIter = 1;
  us.diag.blpIter = 0;
  us.diag.bellmanIter = 0;
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
  free_dmatrix(XZPHIZXinv  , 1, n_theta1, 1, n_theta1);
  
}

void readData
(
 int mcID,
 FILE* inpData,
 FILE* inpRand,
 struct DAT &dat,
 struct RND &rnd,
 struct VAR &var
)
{
  int i, j, k, t, id, imc, itime, iproduct;
  
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      id = idx[t][j];
      fscanf(inpData, "%d %d %d", &imc, &itime, &iproduct);
      fscanf(inpData, "%lf %lf %lf ", &dat.ms[t][j], &dat.ms[t][0], &dat.xi[t][j]);
      for (k=1; k<=n_theta1; k++)
        fscanf(inpData, "%lf ", &var.X[id][k]);
      for (k=1; k<=n_inst; k++)
        fscanf(inpData, "%lf ", &var.Z[id][k]);
      if (imc >= n_mc+mc_init || itime > n_period || iproduct > n_product){
        error("Dimensions of data and estimation model do not match.");
      }
    }
  }
  
  for (i=1; i<=n_person; i++){
    for (k=1; k<=n_theta2; k++)
      fscanf(inpRand, "%lf ", &rnd.nu[i][k]);
  }
}

void readSeed
(
 double ***paramSeed,
 char* filenameSeed
 )
{
  int i, k, s;
  FILE* inpSeed;
  
  inpSeed = fopen(filenameSeed,"r");
  for (i=1; i<=n_mc; i++){
    for (s=1; s<=n_seed; s++){
      for (k=1; k<=n_theta2; k++){
        fscanf(inpSeed, "%lf ", &paramSeed[i][s][k]);
      }
    }
  }
  fclose(inpSeed);
}

void allocateMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
)
{
  par.theta1  = dvector(1, n_theta1);
  par.theta2  = dvector(1, n_theta2);
  par.sd1     = dvector(1, n_theta1);
  par.sd2     = dvector(1, n_theta2);
  
  dat.ms     = dmatrix(1, n_period, 0, n_product);
  dat.xi     = dmatrix(1, n_period, 1, n_product);
  dat.delta  = dvector(1, n_obs);

  fxp.delta         = dvector(1, n_obs);
  fxp.expDeltaIn    = dvector(1, n_obs);
  fxp.expDeltaOut   = dvector(1, n_obs);
  fxp.valuefunIn    = dmatrix(1, n_person, 1, n_grid);
  fxp.valuefunOut   = dmatrix(1, n_person, 1, n_grid);
  fxp.expMu         = darray3d(1, n_person, 1, n_period, 1, n_product);

  rnd.nu = dmatrix(1, n_person, 1, n_theta2);
  
  idx           = imatrix(1, n_period, 1, n_product);
  var.num       = dmatrix(1, n_period, 0, n_product);
  var.value0    = dvector (1, n_period);
  
  var.X         = dmatrix(1, n_obs, 1, n_theta1);
  var.Z         = dmatrix(1, n_obs, 1, n_inst);
  var.Z_t       = dmatrix(1, n_inst, 1, n_obs);
  var.invPHI	  = dmatrix(1, n_inst, 1, n_inst);
  var.XZPHIZXZ  = dmatrix(1, n_theta1, 1, n_obs);
  var.Xtheta1 	= dvector(1, n_obs);
  var.e 				= dvector(1, n_obs);
  var.Ze 				= dvector(1, n_inst);
  var.PHIZe 		= dvector(1, n_inst);
  var.weight    = dvector(1, n_grid);
  var.weightMat = dmatrix(1, n_grid, 1, n_grid);
  
  incv.w      = dvector(1, n_period);
  incv.wGrid  = dvector(1, n_grid);
  incv.lambda = dvector(1, 2);
  
  incv.y       = dvector(1, n_period - 1);
  incv.wy      = dvector(1, 2);
  incv.wlambda = dvector(1, n_period - 1);
  incv.w_lag   = dmatrix(1, n_period - 1, 1, 2);
  incv.w_lag_t = dmatrix(1, 2, 1, n_period - 1);
  incv.ww      = dmatrix(1, 2, 1, 2);
  incv.ww_inv  = dmatrix(1, 2, 1, 2);

  pred.ms            = dmatrix(1, n_period, 1, n_product);
  pred.survivalRate  = dvector(1, n_period);

}

void releaseMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
)
{
  free_dvector(par.theta1,  1, n_theta1);
  free_dvector(par.theta2,  1, n_theta2);
  free_dvector(par.sd1,     1, n_theta1);
  free_dvector(par.sd2,     1, n_theta2);
  
  free_dmatrix(dat.ms,       1, n_period, 0, n_product);
  free_dmatrix(dat.xi,       1, n_period, 1, n_product);
  free_dvector(dat.delta,    1, n_obs);

  free_dvector(fxp.delta,           1, n_obs);
  free_dvector(fxp.expDeltaIn,      1, n_obs);
  free_dvector(fxp.expDeltaOut,     1, n_obs);
  free_dmatrix(fxp.valuefunIn,      1, n_person, 1, n_grid);
  free_dmatrix(fxp.valuefunOut,     1, n_person, 1, n_grid);
  free_darray3d(fxp.expMu,    1, n_person, 1, n_period, 1, n_product);

  free_dmatrix(rnd.nu,        1, n_person, 1, n_theta2);
  
  free_imatrix(idx,           1, n_period, 1, n_product);
  free_dmatrix(var.num,       1, n_period, 0, n_product);
  free_dvector(var.value0,    1, n_period);
  
  free_dmatrix(var.X,         1, n_obs, 1, n_theta1);
  free_dmatrix(var.Z,         1, n_obs, 1, n_inst);
  free_dmatrix(var.Z_t,       1, n_inst, 1, n_obs);
  free_dmatrix(var.invPHI,    1, n_inst, 1, n_inst);
  free_dmatrix(var.XZPHIZXZ,  1, n_theta1, 1, n_obs);
  free_dvector(var.Xtheta1,   1, n_obs);
  free_dvector(var.e,         1, n_obs);
  free_dvector(var.Ze,        1, n_inst);
  free_dvector(var.PHIZe,     1, n_inst);
  free_dvector(var.weight,    1, n_grid);
  free_dmatrix(var.weightMat, 1, n_grid, 1, n_grid);
  
  free_dvector(incv.w,      1, n_period);
  free_dvector(incv.wGrid,  1, n_grid);
  free_dvector(incv.lambda, 1, 2);
  
  free_dvector(incv.y,           1, n_period - 1);
  free_dvector(incv.wy,          1, 2);
  free_dvector(incv.wlambda,     1, n_period - 1);
  free_dmatrix(incv.w_lag,       1, n_period - 1, 1, 2);
  free_dmatrix(incv.w_lag_t,     1, 2, 1, n_period - 1);
  free_dmatrix(incv.ww,          1, 2, 1, 2);
  free_dmatrix(incv.ww_inv,      1, 2, 1, 2);
  
  free_dmatrix(pred.ms,           1, n_period, 1, n_product);
  free_dvector(pred.survivalRate, 1, n_period);

}