#include <stdio.h>
#include "tools.h"
//#include "global.h"
#include "extern.h"
#include "functions.h"
#include "nfpIteration.h"
#include "mpec.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
//#include "knitro.h"

void printOutputMPEC
(
 const int &nnzJ,
 const double * const x,
 double * const obj,
 double * const objGrad,
 struct USERSTRUCT *us
)
{
  double toc;
    
  FILE *outp;
  toc = gettimeofday_sec();
  outp = fopen("outputMPEC.txt","a");
  fprintf(outp, "%6d %6d ",us->mcID, us->seed);
  for (int k=0; k<n_theta; k++)
    fprintf(outp, "%20.14f ",x[k]);
  fprintf(outp, "%16.5f %14.2f\n",*obj,toc-us->tic);
  fclose(outp);
}

void printOutput
(
 char   *filenameEstim,
 int    &mcID,
 int    &seed,
 int    &nStatus,
 double &obj,
 double &toc,
 double *x,
 struct USERSTRUCT &us
)
{
  int k;
  FILE * outp_stat;
  
  outp_stat = fopen(filenameEstim, "a");
  fprintf(outp_stat, "%6d %6d %6d ", mcID,seed,nStatus);
  for (k=0; k<n_theta2; k++)
    fprintf(outp_stat, "%12.6f %12.6f ", x[k], us.par.sd2[k]);
  for (k=0; k<n_theta1; k++)
    fprintf(outp_stat, "%12.6f %12.6f ", x[n_theta2+k], us.par.sd1[k]);
  fprintf(outp_stat, "%14.5f %8d %8d %8d %8d %12lld %14.4f\n", obj,us.diag.iters,
          us.diag.FC_evals,us.diag.GA_evals,us.diag.H_evals,
          us.diag.bellmanIter,toc-us.tic);
  fclose(outp_stat);
  printf("MC experiment %d seed %d complete. Total time elapsed: %f\n",
         mcID, seed, toc-us.tic);

}

void initializeMu
(
 //=======================================================
 adouble **expMu,
 double *nu,
 struct PAR_AD &par_ad,
 struct VAR &var
//=======================================================
)
{
  int j, k, t, id;
  
  adouble mu;
  
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      id = idx[t][j];
      mu = 0;
      for (k=0; k<n_theta2; k++)
        mu += par_ad.theta2[k] * nu[k] * var.X[id][rcID[k]];
      expMu[t][j] = exp(mu);
    }
  }
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
  
  for (i=0; i<n_person; i++){
    for (t=0; t<n_period; t++){
      for (j=0; j<n_product; j++){
        id = idx[t][j];
        mu = 0;
        for (k=0; k<n_theta2; k++)
          mu += par.theta2[k] * rnd.nu[i][k] * var.X[id][rcID[k]];
        expMu[i][t][j] = exp(mu);
      }
    }
  }
}

void initialize
(
 int    &mcID,
 int    &seed,
 double *xInitial,
 double *paramSeed,
 struct USERSTRUCT &us
)
{
  int i, j, k, n, t, s0;
  
  us.mcID = mcID;
  us.seed = seed;

  computeMatrices(us.var.invPHI, us.var.XZPHIZXZ, us.var.X, us.var.Z, us.var.Z_t);

  for (k=0; k<n_theta2; k++){
    xInitial[k] = paramSeed[k];
    us.par.theta2[k] = paramSeed[k];
  }
  
  /*
  par.theta2[0] = 0.5;
  par.theta2[1] = 0.5;
  par.theta2[2] = 0.25;
  par.theta1[0] = 6.0;
  par.theta1[1] = 1.0;
  par.theta1[2] = 1.0;
  par.theta1[3] = 0.5;
  par.theta1[4] = 2.0;

  xInitial[0] = 0.5;
  xInitial[1] = 0.5;
  xInitial[2] = 0.25;
  xInitial[3] = 6.0;
  xInitial[4] = 1.0;
  xInitial[5] = 1.0;
  xInitial[6] = 0.5;
  xInitial[7] = 2.0;
  */
  matvecmul(us.var.X, us.par.theta1, us.var.Xtheta1, n_obs, n_theta1);  // X*theta1
  
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      us.fxp.expDeltaIn[idx[t][j]] = us.dat.ms[t][j] / us.dat.ms0[t];
      //fxp.expDeltaOut[idx[t][j]] = exp(var.Xtheta1[idx[t][j]] + dat.xi[t][j]);
    }
  }
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++){
      us.fxp.valuefunIn[i][s0] = 0;
      //fxp.valuefunOut[i][s0] = dat.valfun[i][s0];
    }
  }

  us.var.tol_blp = 1e-3;
  us.var.tol_bellman = 1e-5;
  us.diag.lipschitz = 0;
  us.diag.lipschitzIter = 1;
  us.diag.blpIter = 0;
  us.diag.bellmanIter = 0;

  nfpBLP(us.fxp,us.dat,us.par,us.rnd,us.var,us.diag,us.incv,us.pred);
  logVector(us.fxp.expDeltaOut, us.fxp.delta, n_obs);
  matvecmul(us.var.XZPHIZXZ, us.fxp.delta, us.par.theta1, n_theta1, n_obs);
  matvecmul(us.var.X, us.par.theta1, us.var.Xtheta1, n_obs, n_theta1);

  k = n_theta2;
  for (i=0; i<n_theta1; i++){
    xInitial[k] = us.par.theta1[i];
    k++;
  }

  /*
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      par.xi[idx[t][j]] = dat.xi[t][j];
    }
  } */
  
  for (n=0; n<n_obs; n++){
    us.par.xi[n] = us.fxp.delta[n] - us.var.Xtheta1[n];
    xInitial[k] = us.par.xi[n];
    k++;
  }
  
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++){
      xInitial[k] = us.fxp.valuefunOut[i][s0];
      k++;
    }
  }

  matvecmul(us.var.Z_t, us.par.xi, us.par.g, n_inst, n_obs);
  for (i=0; i<n_inst; i++){
    xInitial[k] = us.par.g[i];
    k++;
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
  
  X_t         = dmatrix(n_theta1, n_obs);
  PHI 				= dmatrix(n_inst, n_inst);
  XZ 					= dmatrix(n_theta1, n_inst);
  ZX 					= dmatrix(n_inst, n_theta1);
  XZPHI 			= dmatrix(n_theta1, n_inst);
  XZPHIZ 			= dmatrix(n_theta1, n_obs);
  XZPHIZX 		= dmatrix(n_theta1, n_theta1);
  XZPHIZXinv  = dmatrix(n_theta1, n_theta1);
  
  gsl_matrix
  *gsl_invPHI      = gsl_matrix_alloc(n_inst, n_inst),
  *gsl_XZPHIZXinv  = gsl_matrix_alloc(n_theta1, n_theta1);
  
  // Compute 1st-stage GMM weight matrix
  transpose(Z, Z_t, n_obs, n_inst);
  matmul(Z_t, Z, PHI, n_inst, n_obs, n_inst);     				// PHI = Z'Z
  for (i=0; i<n_inst; i++){
    for (j=0; j<n_inst; j++){
      gsl_matrix_set (gsl_invPHI, i, j, PHI[i][j]);
    }
  }
  gsl_linalg_cholesky_decomp(gsl_invPHI);
  gsl_linalg_cholesky_invert(gsl_invPHI);
  for (i=0; i<n_inst; i++){
    for (j=0; j<n_inst; j++){
      invPHI[i][j] = gsl_matrix_get (gsl_invPHI, i, j);
    }
  }
  transpose(X, X_t, n_obs, n_theta1);
  matmul(X_t, Z, XZ, n_theta1, n_obs, n_inst);    				// XZ = X_t*Z = X'Z
  transpose(XZ, ZX, n_theta1, n_inst);                      	// ZX = (XZ)' = Z'X
  matmul(XZ, invPHI, XZPHI, n_theta1, n_inst, n_inst);   // XZPHI = (X'Z)inv(Z'Z)
  matmul(XZPHI, ZX, XZPHIZX, n_theta1, n_inst, n_theta1); // XZPHIZX = (X'Z)inv(Z'Z)(Z'X)
  for (i=0; i<n_theta1; i++){
    for (j=0; j<n_theta1; j++){
      gsl_matrix_set (gsl_XZPHIZXinv, i, j, XZPHIZX[i][j]);
    }
  }
  gsl_linalg_cholesky_decomp(gsl_XZPHIZXinv);
  gsl_linalg_cholesky_invert(gsl_XZPHIZXinv);
  for (i=0; i<n_theta1; i++){
    for (j=0; j<n_theta1; j++){
      XZPHIZXinv[i][j] = gsl_matrix_get (gsl_XZPHIZXinv, i, j);
    }
  }
  matmul(XZPHI, Z_t, XZPHIZ, n_theta1, n_inst, n_obs);    // XZPHIZ = (X'Z)inv(Z'Z)Z'
  // XZPHIZXZ = inv((X'Z)inv(Z'Z)(Z'X))(X'Z)inv(Z'Z)Z'
  matmul(XZPHIZXinv, XZPHIZ, XZPHIZXZ, n_theta1, n_theta1, n_obs);
  
  gsl_matrix_free (gsl_invPHI);
  gsl_matrix_free (gsl_XZPHIZXinv);
  
  free_dmatrix(X_t);
  free_dmatrix(PHI);
  free_dmatrix(XZ);
  free_dmatrix(ZX);
  free_dmatrix(XZPHI);
  free_dmatrix(XZPHIZ);
  free_dmatrix(XZPHIZX);
  free_dmatrix(XZPHIZXinv);
  
}

void readData
(
 int    &mcID,
 FILE*  inpData,
 FILE*  inpRand,
 char*  filenameValfun,
 struct DAT &dat,
 struct RND &rnd,
 struct VAR &var
 )
{
  int i, j, k, t, s0, id, imc, itime, iproduct, status;
  
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      id = idx[t][j];
      status = fscanf(inpData, "%d %d %d ", &imc, &itime, &iproduct);
      status = fscanf(inpData, "%lf %lf %lf ", &dat.ms[t][j], &dat.ms0[t], &dat.xi[t][j]);
      for (k=0; k<n_theta1; k++)
        status = fscanf(inpData, "%lf ", &var.X[id][k]);
      for (k=0; k<n_inst; k++)
        status = fscanf(inpData, "%lf ", &var.Z[id][k]);

      if (imc>n_mc || itime>n_period || iproduct>n_product || imc!=mcID+1){
        error("Dimensions of data and estimation model do not match.");
      }
    }
  }
  
  FILE* inpVfun = fopen(filenameValfun,"r");
  for (i=0; i<n_person; i++){
    for (k=0; k<n_theta2; k++)
      status = fscanf(inpRand, "%lf ", &rnd.nu[i][k]);
    for (s0=0; s0<n_grid; s0++)
      status = fscanf(inpVfun, "%lf ", &dat.valfun[i][s0]);
  }
  fclose(inpVfun);
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
  for (i=0; i<n_mc; i++){
    for (s=0; s<n_seed; s++){
      for (k=0; k<n_theta2; k++){
        fscanf(inpSeed, "%lf ", &paramSeed[i][s][k]);
      }
    }
  }
  fclose(inpSeed);
}

void initializeKnitro
(
 int    &m,
 int    &n,
 int    &nnzJ,
 int    &nnzH,
 int    *cType,
 double *xInitial,
 double *xLoBnds,
 double *xUpBnds,
 double *cLoBnds,
 double *cUpBnds,
 double *paramBound,
 struct USERSTRUCT &us
)
{
  int i,j,k;
  
  for (i=0; i<n_theta2; i++){
    xLoBnds[i] = 0.0;
    xUpBnds[i] = paramBound[i];
  }
  for (i=n_theta2; i<n; i++){
    xLoBnds[i] = -KTR_INFBOUND;
    xUpBnds[i] = KTR_INFBOUND;
  }
  for (j=0; j<m; j++){
    cType[j] = KTR_CONTYPE_GENERAL;
    cLoBnds[j] = 0.0;
    cUpBnds[j] = 0.0;
  }
  for (j=n_obs+dimValfun; j<m; j++){
    cType[j] = KTR_CONTYPE_LINEAR;
  }
  
  us.tagJacobian = 0;     // index of Adolc tape for Jacobian
  us.tagHessian  = 1;     // index of Adolc tape for Hessian
  us.optionsJac[0] = 0;
  us.optionsJac[1] = 0;
  us.optionsJac[2] = 0;
  us.optionsJac[3] = 0;
  
  nnzJ = NULL;
  us.nnzHconst = NULL;
  
  double *constraint  = dvector(dimConst);
  us.tic = gettimeofday_sec();
  // Generate permanent tapes
  tapeJacobian(constraint, xInitial, us);
  free_dvector(constraint);
  
  // Detect sparsity patterns
  sparse_jac(us.tagJacobian, dimConst, dimOptim,0, xInitial, &nnzJ,
             &us.jacIndexCons, &us.jacIndexVars, &us.jac, us.optionsJac);
  //printf("Constraint Jacobian sparsity pattern detected.\n");
  
  us.jacIndexVarsInt  = ivector(nnzJ);
  us.jacIndexConsInt  = ivector(nnzJ);
  
  // Copy Adolc jacobian index to KNITRO
  for (k=0; k<nnzJ; k++){
    us.jacIndexConsInt[k] = us.jacIndexCons[k];
    us.jacIndexVarsInt[k] = us.jacIndexVars[k];
  }
  
  // UPPER TRIANGULAR Hessian sparsity
  k = 0;
  for (i=0; i<n_theta+n_obs+dimValfun; i++){
    for (j=i; j<n_theta+n_obs+dimValfun; j++){
      us.hessIndexRows[k] = i;
      us.hessIndexCols[k] = j;
      k++;
    }
  }
  us.nnzHconst = k;
  for (i=n_theta+n_obs+dimValfun; i<dimOptim; i++){
    for (j=i; j<dimOptim; j++){
      us.hessIndexRows[k] = i;
      us.hessIndexCols[k] = j;
      k++;
    }
  }
  nnzH = k;
}

void allocateKnitro
(
 int    &m,
 int    &n,
 //int    &nnzJ,
 int    &nnzH,
 int    **cType,
 double **xLoBnds,
 double **xUpBnds,
 double **xInitial,
 double **x,
 double **cLoBnds,
 double **cUpBnds,
 double **lambda,
 struct USERSTRUCT &us
)
{
  *xLoBnds  = dvector(n);
  *xUpBnds  = dvector(n);
  *xInitial = dvector(n);
  *x        = dvector(n);
  *cType    = ivector(m);
  *cLoBnds  = dvector(m);
  *cUpBnds  = dvector(m);
  *lambda   = dvector(m+n);

  allocateActiveMemory(us.par_ad, us.fxp_ad, us.var_ad, us.incv_ad, us.pred_ad);
  
  us.c_ad  = advector(dimConst);
  us.xcopy = dvector(dimOptim);

  us.jac = NULL;
  us.hessConstraint = dmatrix(dimOptim, dimOptim);
  us.jacIndexCons = NULL;
  us.jacIndexVars = NULL;
  us.hessIndexRows = ivector(nnzH);
  us.hessIndexCols = ivector(nnzH);
  
}

void freeKnitro
(
 int    **cType,
 double **xLoBnds,
 double **xUpBnds,
 double **xInitial,
 double **x,
 double **cLoBnds,
 double **cUpBnds,
 double **lambda,
 struct USERSTRUCT &us
)
{
  free_dvector(*xLoBnds);
  free_dvector(*xUpBnds);
  free_dvector(*xInitial);
  free_dvector(*x);
  free_ivector(*cType);
  free_dvector(*cLoBnds);
  free_dvector(*cUpBnds);
  free_dvector(*lambda);
 
  releaseActiveMemory(us.par_ad, us.fxp_ad, us.var_ad, us.incv_ad, us.pred_ad);
  
  free_advector(us.c_ad);
  free_dvector(us.xcopy);

  free(us.jac);
  free_dmatrix(us.hessConstraint);
  free(us.jacIndexVars);
  free(us.jacIndexCons);
  free_ivector(us.hessIndexRows);
  free_ivector(us.hessIndexCols);

  free_ivector(us.jacIndexVarsInt);
  free_ivector(us.jacIndexConsInt);

}

void allocateMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct INCV &incv,
 struct PRED &pred
)
{
  par.theta1  = dvector(n_theta1);
  par.theta2  = dvector(n_theta2);
  par.sd1     = dvector(n_theta1);
  par.sd2     = dvector(n_theta2);
  par.xi      = dvector(n_obs);
  par.valfun  = dmatrix(n_person, n_grid);
  par.g       = dvector(n_inst);

  dat.ms     = dmatrix(n_period, n_product);
  dat.ms0    = dvector(n_period);
  dat.xi     = dmatrix(n_period, n_product);
  dat.delta  = dvector(n_obs);
  dat.valfun = dmatrix(n_person, n_grid);
  
  fxp.delta         = dvector(n_obs);
  fxp.expDeltaIn    = dvector(n_obs);
  fxp.expDeltaOut   = dvector(n_obs);
  fxp.valuefunIn    = dmatrix(n_person, n_grid);
  fxp.valuefunOut   = dmatrix(n_person, n_grid);
  fxp.expMu         = darray3d(n_person, n_period, n_product);
  
  rnd.nu = dmatrix(n_person, n_theta2);
  
  idx           = imatrix(n_period, n_product);
  
  var.X         = dmatrix(n_obs, n_theta1);
  var.Z         = dmatrix(n_obs, n_inst);
  var.Z_t       = dmatrix(n_inst, n_obs);
  var.invPHI	  = dmatrix(n_inst, n_inst);
  var.XZPHIZXZ  = dmatrix(n_theta1, n_obs);
  
  var.num       = dmatrix(n_period, n_product+1);
  var.value0    = dvector (n_period);
  var.Xtheta1   = dvector(n_obs);
  var.e 			  = dvector(n_obs);
  var.Ze 			  = dvector(n_inst);
  var.PHIZe 	  = dvector(n_inst);
  var.weightMat = dmatrix(n_grid, n_grid);

  incv.w       = dvector(n_period);
  incv.wGrid   = dvector(n_grid);
  incv.lambda  = dvector(2);
  incv.y       = dvector(n_period-1);
  incv.wy      = dvector(2);
  incv.wlambda = dvector(n_period-1);
  incv.w_lag   = dmatrix(n_period-1, 2);
  incv.w_lag_t = dmatrix(2, n_period-1);
  incv.ww      = dmatrix(2, 2);
  incv.ww_inv  = dmatrix(2, 2);
  
  pred.ms           = dmatrix(n_period, n_product);
  pred.survivalRate = dvector(n_period);

}

void allocateActiveMemory
(
 struct PAR_AD &par_ad,
 struct FXP_AD &fxp_ad,
 struct VAR_AD &var_ad,
 struct INCV_AD &incv_ad,
 struct PRED_AD &pred_ad
)
{
  par_ad.theta1  = advector(n_theta1);
  par_ad.theta2  = advector(n_theta2);
  par_ad.xi      = advector(n_obs);
  par_ad.valfun  = admatrix(n_person, n_grid);
  par_ad.g       = advector(n_inst);
  
  fxp_ad.expDeltaIn  = advector(n_obs);
  fxp_ad.expMu       = admatrix(n_period, n_product);
  fxp_ad.valuefunOut = admatrix(n_person, n_grid);

  var_ad.num        = admatrix(n_period, n_product+1);
  var_ad.value0     = advector (n_period);
  var_ad.Xtheta1    = advector(n_obs);
  var_ad.e 				  = advector(n_obs);
  var_ad.Ze 				= advector(n_inst);
  var_ad.PHIZe 		  = advector(n_inst);
  var_ad.weightMat  = admatrix(n_grid, n_grid);
  
  incv_ad.w       = advector(n_period);
  incv_ad.wGrid   = advector(n_grid);
  incv_ad.lambda  = advector(2);
  incv_ad.y       = advector(n_period-1);
  incv_ad.wy      = advector(2);
  incv_ad.wlambda = advector(n_period-1);
  incv_ad.w_lag   = admatrix(n_period-1, 2);
  incv_ad.w_lag_t = admatrix(2, n_period-1);
  incv_ad.ww      = admatrix(2, 2);
  incv_ad.ww_inv  = admatrix(2, 2);

  pred_ad.ms           = admatrix(n_period, n_product);
  pred_ad.survivalRate = advector(n_period);

}

void releaseMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct INCV &incv,
 struct PRED &pred
)
{
  free_dvector(par.theta1);
  free_dvector(par.theta2);
  free_dvector(par.sd1);
  free_dvector(par.sd2);
  free_dvector(par.xi);
  free_dmatrix(par.valfun);
  free_dvector(par.g);

  free_dmatrix(dat.ms);
  free_dvector(dat.ms0);
  free_dmatrix(dat.xi);
  free_dvector(dat.delta);
  free_dmatrix(dat.valfun);
  
  free_dvector(fxp.delta);
  free_dvector(fxp.expDeltaIn);
  free_dvector(fxp.expDeltaOut);
  free_dmatrix(fxp.valuefunIn);
  free_dmatrix(fxp.valuefunOut);
  free_darray3d(fxp.expMu);

  free_dmatrix(rnd.nu);
  
  free_imatrix(idx);

  free_dmatrix(var.X);
  free_dmatrix(var.Z);
  free_dmatrix(var.Z_t);
  free_dmatrix(var.invPHI);
  free_dmatrix(var.XZPHIZXZ);
  
  free_dmatrix(var.num);
  free_dvector(var.value0);
  free_dvector(var.Xtheta1);
  free_dvector(var.e);
  free_dvector(var.Ze);
  free_dvector(var.PHIZe);
  free_dmatrix(var.weightMat);
  
  free_dvector(incv.w);
  free_dvector(incv.wGrid);
  free_dvector(incv.lambda);
  free_dvector(incv.y);
  free_dvector(incv.wy);
  free_dvector(incv.wlambda);
  free_dmatrix(incv.w_lag);
  free_dmatrix(incv.w_lag_t);
  free_dmatrix(incv.ww);
  free_dmatrix(incv.ww_inv);
  
  free_dmatrix(pred.ms);
  free_dvector(pred.survivalRate);
  
}

void releaseActiveMemory
(
 struct PAR_AD &par_ad,
 struct FXP_AD &fxp_ad,
 struct VAR_AD &var_ad,
 struct INCV_AD &incv_ad,
 struct PRED_AD &pred_ad
)
{
  free_advector(par_ad.theta1);
  free_advector(par_ad.theta2);
  free_advector(par_ad.xi);
  free_admatrix(par_ad.valfun);
  free_advector(par_ad.g);
  
  free_advector(fxp_ad.expDeltaIn);
  free_admatrix(fxp_ad.expMu);
  free_admatrix(fxp_ad.valuefunOut);
  
  free_admatrix(var_ad.num);
  free_advector(var_ad.value0);
  free_advector(var_ad.Xtheta1);
  free_advector(var_ad.e);
  free_advector(var_ad.Ze);
  free_advector(var_ad.PHIZe);
  free_admatrix(var_ad.weightMat);
  
  free_advector(incv_ad.w);
  free_advector(incv_ad.wGrid);
  free_advector(incv_ad.lambda);
  free_advector(incv_ad.y);
  free_advector(incv_ad.wy);
  free_advector(incv_ad.wlambda);
  free_admatrix(incv_ad.w_lag);
  free_admatrix(incv_ad.w_lag_t);
  free_admatrix(incv_ad.ww);
  free_admatrix(incv_ad.ww_inv);
  
  free_admatrix(pred_ad.ms);
  free_advector(pred_ad.survivalRate);
  
}