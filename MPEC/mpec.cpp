#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "extern.h"
#include "tools.h"
#include "functions.h"
#include "mpec.h"
#include "nfpIteration.h"
//#include "knitro.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void constraintOnly
(
 double * const c,
 const double * const x,
 struct USERSTRUCT *us
 )
{
  int i, j, k, t, s0;
  
  for (i=0; i<n_theta2; i++)
    us->par.theta2[i] = x[i];
  for (i=0; i<n_theta1; i++)
    us->par.theta1[i] = x[n_theta2+i];
  for (i=0; i<n_obs; i++)
    us->par.xi[i] = x[n_theta+i];
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++)
      us->par.valfun[i][s0] = x[n_theta+n_obs+n_grid*i+s0];
  }
  for (i=0; i<n_inst; i++)
    us->par.g[i] = x[n_theta+n_obs+n_grid*n_person+i];
  
  matvecmul(us->var.X, us->par.theta1, us->var.Xtheta1, n_obs, n_theta1);
  
  for (i=0; i<n_obs; i++)
    us->fxp.expDeltaIn[i] = exp(us->var.Xtheta1[i] + us->par.xi[i]);
  
  initializeMu(us->fxp.expMu, us->par, us->rnd, us->var);
  
  for (i=0; i<n_person; i++){
    computeInclusiveValue(us->incv, us->fxp.expDeltaIn, us->fxp.expMu[i], us->var);
    transProb(us->var, us->incv);
    bellmanOperator(us->fxp.valuefunOut[i], us->par.valfun[i], us->var, us->incv);
    computeExpectedValuefun(us->var.value0, us->fxp.valuefunOut[i], us->incv);
    computeChoicePath(us->pred, us->var);
  }
  k = 0;
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      us->pred.ms[t][j] /= n_person;
      c[k] = log(us->dat.ms[t][j]/us->pred.ms[t][j]);
      k++;
    }
  }
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++){
      c[k] = us->fxp.valuefunOut[i][s0] - us->par.valfun[i][s0];
      k++;
    }
  }
  matvecmul(us->var.Z_t, us->par.xi, us->var.Ze, n_inst, n_obs);
  
  for (i=0; i<n_inst; i++){
    c[k] = us->par.g[i] - us->var.Ze[i];
    k++;
  }
}

void tapeJacobian
(
 double * const c,
 const double * const x,
 struct USERSTRUCT &us
)
{
  int i, j, k, t, s0;
  //size_t tapeStats[STAT_SIZE];
  
  trace_on(us.tagJacobian);

  for (i=0; i<n_theta2; i++)
    us.par_ad.theta2[i] <<= x[i];
  for (i=0; i<n_theta1; i++)
    us.par_ad.theta1[i] <<= x[n_theta2+i];
  for (i=0; i<n_obs; i++)
    us.par_ad.xi[i] <<= x[n_theta+i];
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++)
      us.par_ad.valfun[i][s0] <<= x[n_theta+n_obs+n_grid*i+s0];
  }
  for (i=0; i<n_inst; i++)
    us.par_ad.g[i] <<= x[n_theta+n_obs+dimValfun+i];
  
  matvecmul(us.var.X, us.par_ad.theta1, us.var_ad.Xtheta1, n_obs, n_theta1);
  
  for (i=0; i<n_obs; i++)
    us.fxp_ad.expDeltaIn[i] = exp(us.var_ad.Xtheta1[i] + us.par_ad.xi[i]);
  
  for (i=0; i<n_person; i++){
    initializeMu(us.fxp_ad.expMu, us.rnd.nu[i], us.par_ad, us.var);
    computeInclusiveValue(us.incv_ad, us.fxp_ad.expDeltaIn, us.fxp_ad.expMu, us.var_ad);
    transProb_AD(us.var_ad, us.incv_ad);
    bellmanOperator(us.fxp_ad.valuefunOut[i], us.par_ad.valfun[i], us.var_ad, us.incv_ad);
    computeExpectedValuefun(us.var_ad.value0, us.fxp_ad.valuefunOut[i], us.incv_ad);
    computeChoicePath_AD(us.pred_ad, us.var_ad);
  }
  k = 0;
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      us.pred_ad.ms[t][j] /= n_person;
      us.c_ad[k] = log(us.dat.ms[t][j]/us.pred_ad.ms[t][j]);
      //us.c_ad[k] = us.dat.ms[t][j] - us.pred_ad.ms[t][j];
      k++;
    }
  }
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++){
      us.c_ad[k] = us.fxp_ad.valuefunOut[i][s0] - us.par_ad.valfun[i][s0];
      k++;
    }
  }
  matvecmul(us.var.Z_t, us.par_ad.xi, us.var_ad.Ze, n_inst, n_obs);
  
  for (i=0; i<n_inst; i++){
    us.c_ad[k] = us.par_ad.g[i] - us.var_ad.Ze[i];
    k++;
  }
  
  for (i=0; i<dimConst; i++){
    us.c_ad[i] >>= c[i];
  }
  trace_off();
  /*
  tapestats(tag,tapeStats);
  printf("N of indep: %lu\n", tapeStats[0]);
  printf("N of dep: %lu\n", tapeStats[1]);
  printf("N of maxlive: %lu\n", tapeStats[2]);
  printf("Tape size: %lu\n", tapeStats[3]);
  */
  //jacobian(us.tag, dimConst, dimOptim, x, us.jac);
  
}

void constraintJacobian
(
 const int &nnzJ,
 double * const c,
 double * jac,
 const double * const x,
 struct USERSTRUCT &us
)
{
  //tapeJacobian(nnzJ, c, x, us);
  /*
  jacobian(us.tagJacobian, dimConst, dimOptim, x, us.jac);
  for (int k=0; k<nnzJ; k++){
    jac[k] = us.jac[us.jacIndexCons[k]][us.jacIndexVars[k]];
  }
  */
  int nnz = nnzJ;
  sparse_jac(us.tagJacobian,dimConst,dimOptim,1,x,&nnz,&(us.jacIndexCons),&(us.jacIndexVars),&jac,us.optionsJac);
  
}

void tapeHessian
(
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
)
{
  int i, j, k, t, s0;
  double constraint;
  
  // Constraint Hessian
  trace_on(us->tagHessian);
  
  for (i=0; i<n_theta2; i++)
    us->par_ad.theta2[i] <<= x[i];
  for (i=0; i<n_theta1; i++)
    us->par_ad.theta1[i] <<= x[n_theta2+i];
  for (i=0; i<n_obs; i++)
    us->par_ad.xi[i] <<= x[n_theta+i];
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++)
      us->par_ad.valfun[i][s0] <<= x[n_theta+n_obs+n_grid*i+s0];
  }
  for (i=0; i<n_inst; i++)
    us->par_ad.g[i] <<= x[n_theta+n_obs+dimValfun+i];

  matvecmul(us->var.X, us->par_ad.theta1, us->var_ad.Xtheta1, n_obs, n_theta1);
  
  for (i=0; i<n_obs; i++)
    us->fxp_ad.expDeltaIn[i] = exp(us->var_ad.Xtheta1[i] + us->par_ad.xi[i]);
  
  for (i=0; i<n_person; i++){
    initializeMu(us->fxp_ad.expMu, us->rnd.nu[i], us->par_ad, us->var);
    computeInclusiveValue(us->incv_ad, us->fxp_ad.expDeltaIn, us->fxp_ad.expMu, us->var_ad);
    transProb_AD(us->var_ad, us->incv_ad);
    bellmanOperator(us->fxp_ad.valuefunOut[i], us->par_ad.valfun[i], us->var_ad, us->incv_ad);
    computeExpectedValuefun(us->var_ad.value0, us->fxp_ad.valuefunOut[i], us->incv_ad);
    computeChoicePath_AD(us->pred_ad, us->var_ad);
  }
  k = 0;
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      us->pred_ad.ms[t][j] /= n_person;
      us->c_ad[k] = log(us->dat.ms[t][j]/us->pred_ad.ms[t][j]);
      //us->c_ad[k] = us->dat.ms[t][j] - us->pred_ad.ms[t][j];
      k++;
    }
  }
  for (i=0; i<n_person; i++){
    for (s0=0; s0<n_grid; s0++){
      us->c_ad[k] = us->fxp_ad.valuefunOut[i][s0] - us->par_ad.valfun[i][s0];
      k++;
    }
  }
  matvecmul(us->var.Z_t, us->par_ad.xi, us->var_ad.Ze, n_inst, n_obs);
  
  for (i=0; i<n_inst; i++){
    us->c_ad[k] = us->par_ad.g[i] - us->var_ad.Ze[i];
    k++;
  }
  
  us->constraint = 0;
  for (i=0; i<dimConst; i++){
    //us->constraint += us->par_ad.lambda[i] * us->c_ad[i];
    us->constraint += lambda[i] * us->c_ad[i];
  }
  us->constraint >>= constraint;
  trace_off();

}

void lagrangianHessian
(
 const int &nnzH,
 double * const hess,
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
)
{
  int i, j, k;
  lagrangianHessianNoF(nnzH, hess, x, lambda, us);
  
  k = us->nnzHconst;
  for (i=0; i<n_inst; i++){
    for (j=i; j<n_inst; j++){
      hess[k] = 2*us->var.invPHI[i][j];
      k++;
    }
  }
}

void lagrangianHessianNoF
(
 const int &nnzH,
 double * const hess,
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
)
{
  int i, k;
  
  // Generate temporary tape for Hessian
  tapeHessian(x, lambda, us);

  // Compute LOWER TRIANGULAR Hessian for the constraint
  for (i=0; i<dimOptim; i++)
    us->xcopy[i] = x[i];

  hessian(us->tagHessian, dimOptim, us->xcopy, us->hessConstraint);
  
  // Copy to UPPER TRIANGULAR Hessian for KNITRO
  for (k=0; k<us->nnzHconst; k++)
    hess[k] = us->hessConstraint[us->hessIndexCols[k]][us->hessIndexRows[k]];
  for (k=us->nnzHconst; k<nnzH; k++)
    hess[k] = 0;
}

void gmmObj_and_gradient
(
 //=======================================================
 double *obj,
 double *grad,
 const double * const x,
 struct USERSTRUCT *us
//=======================================================
)
{
  int i;
  for (i=0; i<n_inst; i++)
    us->par.g[i] = x[n_theta+n_obs+dimValfun+i];
  
  matvecmul(us->var.invPHI, us->par.g, us->var.PHIZe, n_inst, n_inst);   // PHIZe = inv(Z'Z)(Z'e)
  *obj = 0.0;
  for (i=0; i<n_theta+n_obs+dimValfun; i++)
    grad[i] = 0.0;
  
  for (i=0; i<n_inst; i++){
    *obj += us->var.PHIZe[i] * us->par.g[i];
    grad[n_theta+n_obs+dimValfun+i] = 2*us->var.PHIZe[i];
  }
  
}

int callbackEvalFC
(
 const int evalRequestCode,
 const int n,
 const int m,
 const int nnzJ,
 const int nnzH,
 const double * const x,
 double * const lambda,
 double * const obj,
 double * const c,
 double * const objGrad,
 double * const jac,
 double * const hessian,
 double * const hessVector,
 void * userParams
)
{
  struct USERSTRUCT *us = (struct USERSTRUCT*) userParams;
  
  if (evalRequestCode == KTR_RC_EVALFC){
    gmmObj_and_gradient(obj, objGrad, x, us);
    //constraintOnly(c,x,us);

    for (int i=0; i<n; i++)
      us->xcopy[i] = x[i];
    function(us->tagJacobian,m,n,us->xcopy,c);
    //printOutputMPEC(x, obj, us);
    return 0;
  }
  else {
    printf("callbackEvalFC incorrectly called with eval code %d\n",evalRequestCode);
    return -1;
  }
}

int callbackEvalGA
(
 const int evalRequestCode,
 const int n,
 const int m,
 const int nnzJ,
 const int nnzH,
 const double * const x,
 double * const lambda,
 double * const obj,
 double * const c,
 double * const objGrad,
 double * const jac,
 double * const hessian,
 double * const hessVector,
 void * userParams
 )
{
  struct USERSTRUCT *us = (struct USERSTRUCT*) userParams;
  
  if (evalRequestCode == KTR_RC_EVALGA){
    gmmObj_and_gradient(obj, objGrad, x, us);
    constraintJacobian(nnzJ,c,jac,x,*us);
    printOutputMPEC(nnzJ,x, obj, objGrad,us);
    return 0;
  }
  else {
    printf("callbackEvalGA incorrectly called with eval code %d\n",evalRequestCode);
    return -1;
  }
}

int callbackEvalH
(
 const int evalRequestCode,
 const int n,
 const int m,
 const int nnzJ,
 const int nnzH,
 const double * const x,
 double * const lambda,
 double * const obj,
 double * const c,
 double * const objGrad,
 double * const jac,
 double * const hessian,
 double * const hessVector,
 void * userParams
 )
{
  struct USERSTRUCT *us = (struct USERSTRUCT*) userParams;
  
  if (evalRequestCode == KTR_RC_EVALH){
    lagrangianHessian(nnzH,hessian,x,lambda,us);
    return 0;
  }

  else if (evalRequestCode == KTR_RC_EVALH_NO_F){
    lagrangianHessianNoF(nnzH,hessian,x,lambda,us);
    return 0;
  }
  else {
    printf("callbackEvalH incorrectly called with eval code %d\n",evalRequestCode);
    return -1;
  }
}
