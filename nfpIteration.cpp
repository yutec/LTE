#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "extern.h"
#include "global.h"
#include "tools.h"
#include "functions.h"
#include "nfpIteration.h"
//#include "knitro.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

double computeLogQuasiposterior
(//=======================================================
 struct PAR &par,
 struct FXP &fxp,
 struct RND &rnd,
 struct VAR &var,
 struct DAT &dat,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  int i, n;
  double obj;
  
  nfpBLP(fxp, dat,par,rnd,var,diag,incv,pred);
  if (diag.exit<0){
    copyVector(fxp.expDeltaIn, fxp.expDeltaOut, n_obs);
    copyMatrix(fxp.valuefunIn, fxp.valuefunOut, n_person, n_grid);
    return 1e12;
  }
  
  logVector(fxp.expDeltaOut, fxp.delta, n_obs);
  matvecmul(var.XZPHIZXZ, fxp.delta, par.theta1, n_theta1, n_obs);
  matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
  for (n=1; n<=n_obs; n++){
    var.e[n] = fxp.delta[n] - var.Xtheta1[n];
  }
  
  matvecmul(var.Z_t, var.e, var.Ze, n_inst, n_obs);            // Z'e
  matvecmul(var.invPHI, var.Ze, var.PHIZe, n_inst, n_inst);   // PHIZe = inv(Z'Z)(Z'e)
  
  obj = 0;
  for (i=1; i<=n_inst; i++)
    obj += var.Ze[i] * var.PHIZe[i];
  
  return(obj);
}

void computeExpectedValuefun
(//=======================================================
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
  sumWeight;   // sum of transition prob weight
  
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

void nfpBellman
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv
//=======================================================
)
{
  int s0, iter = 0;
  double dist = DBL_MAX;
  
  // Solve for exact value function
  while (dist > tol_blm && iter <= 3000){
    bellmanOperator(valuefunOut, valuefunIn, var, incv);
    dist = 0;
    for (s0=1; s0<=n_grid; s0++){
      dist = max(dist, fabs(valuefunIn[s0] - valuefunOut[s0]));
      valuefunIn[s0] = valuefunOut[s0];
    }
    iter++;
    diag.bellmanIter++;
  }
  
}

void nfpBLP
(//=======================================================
 struct FXP &fxp,
 struct DAT &dat,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  int n, it;
  double dist, distPrev=1;
  
  initializeMu(fxp.expMu, par, rnd, var);
  
  it = 0;
  dist = DBL_MAX;
  diag.exit = 0;
  while (dist > tol_blp && it <= 2000 && diag.exit>=0){
    blpMap(fxp.expDeltaOut, fxp.valuefunOut, fxp.expDeltaIn, fxp.valuefunIn, fxp.expMu,
           dat, par, var, diag, incv, pred);
    
    dist = 0;
    for (n=1; n<=n_obs; n++){
      dist = fmax(dist, fabs(log(fxp.expDeltaIn[n]/fxp.expDeltaOut[n])));
      fxp.expDeltaIn[n] = fxp.expDeltaOut[n];
    }
    it++;
    
    if (diagPrint == 1){
      printf("Iter: %5d  Dist: %18.15f  BLP mod: %f\n",it, dist, dist/distPrev);
    }
    
    if (it > 2){
      diag.lipschitz += dist / distPrev;
      diag.lipschitzIter++;
    }
    diag.blpIter++;
    distPrev = dist;
    
  }
  if (diagPrint == 1){
    printf("--------------------------------\n");
  }
  
}


void bellmanOperator
(//=======================================================
  double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct INCV &incv
//=======================================================
)
{
  int s0, s1;
  double value0;
  
  for (s0=1; s0<=n_grid; s0++){
    // value0: expected value function of delaying purchase
    value0 = 0;
    for (s1=1; s1<=n_grid; s1++){
      value0 += valuefunIn[s1] * var.weightMat[s0][s1];
    }
    value0 *= beta;
    valuefunOut[s0] = log(1 + exp(incv.wGrid[s0] - value0)) + value0;
  } // end of s0 loop
  
}

void blpMap
(
 //=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct DAT &dat,
 struct PAR &par,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  int i, j, t;
  
  for (t=1; t<=n_period; t++) {
    for (j=1; j<=n_product; j++) {
      pred.ms[t][j] = 0;
    }
  }
  
  for (i=1; i<=n_person; i++){
    computeInclusiveValue(incv, expDeltaIn, expMuIn[i], var);
    transProb(var, incv);
    nfpBellman(valuefunOut[i], valuefunIn[i], var, diag, incv);
    computeExpectedValuefun(var.value0, valuefunOut[i], var.weight, incv);
    computeChoicePath(pred, var);
  }
  diag.exit = deltaUpdate(expDeltaOut, expDeltaIn, dat, pred);
}

void computeInclusiveValue
(//=======================================================
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

void computeChoicePath
(
 //=======================================================
 struct PRED &pred,
 struct VAR &var
//=======================================================
)
{
  int j, t;
  double den;
  
  // Initialize survivalRateard (var.survivalRate) at t=1 for all i
  pred.survivalRate[1] = 1;
  
  // Compute predicted market share given deltaIn
  for (t=1; t<=n_period; t++){
    var.num[t][0] = exp(beta * var.value0[t]);
    den = 0;
    for (j=0; j<=n_product; j++){
      den += var.num[t][j];
    }
    
    // Add choice prob * (1-survivalRateard) to predicted market share
    for (j=1; j<=n_product; j++){
      pred.ms[t][j] += var.num[t][j] / den * pred.survivalRate[t];
    }
    
    // Update the next-period survivalRateard
    if (t < n_period){
      pred.survivalRate[t+1] = var.num[t][0] / den * pred.survivalRate[t];
    }
    
  } // t loop
}

int callback
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
  int k;
  double toc;
  FILE *outp;
  struct USERSTRUCT *us = (struct USERSTRUCT*) userParams;
  
  if (evalRequestCode == KTR_RC_EVALFC){
    for (k=0; k<n_theta2; k++)
      us->par.theta2[k+1] = x[k];
    
    *obj = computeLogQuasiposterior(us->par, us->fxp, us->rnd, us->var, us->dat,
                                    us->diag, us->incv, us->pred);
    toc = gettimeofday_sec();
    outp = fopen("outputNFP.txt","a");
    fprintf(outp, "%6d %6d ",us->mcID, us->seed);
    for (k=1; k<=n_theta2; k++)
      fprintf(outp, "%20.14f ",us->par.theta2[k]);
    for (k=1; k<=n_theta1; k++)
      fprintf(outp, "%20.14f ", us->par.theta1[k]);
    fprintf(outp, "%21.14f %14lld %8.4f %14.2f\n", *obj, us->diag.bellmanIter,
            us->diag.lipschitz/us->diag.lipschitzIter, toc-us->tic);
    fclose(outp);
    return 0;
  }
  else{
    printf ("Wrong evalRequestCode in callback function.\n");
    return -1;
  }
}

void stdError
(
 int    &nStatus,
 double *x,
 double &toc,
 struct FXP &fxp,
 struct DAT &dat,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred,
 KTR_context *kc
 )
{
  int i, j, k, l, n;
  double
  epsilon = 1e-6,
  *xi0,
  *xiP,
  **dXi,
  **dZe,
  **dZeTrans,
  **PhiZe,
  **ZePhiZe;
  
  gsl_matrix *gsl_ZePhiZeInv;
  
  struct PAR parP;
  
  if (diag.exit!=0)
    nStatus = -1;
  
  if (nStatus == 0){
    gsl_ZePhiZeInv = gsl_matrix_alloc(n_theta, n_theta);
    xi0      = dvector(1, n_obs);
    xiP      = dvector(1, n_obs);
    dXi      = dmatrix(1, n_obs, 1, n_theta);
    dZe      = dmatrix(1, n_inst, 1, n_theta);
    dZeTrans = dmatrix(1, n_theta, 1, n_inst);
    PhiZe    = dmatrix(1, n_inst, 1, n_theta);
    ZePhiZe  = dmatrix(1, n_theta, 1, n_theta);
    
    parP.theta1 = dvector(1, n_theta1);
    parP.theta2 = dvector(1, n_theta2);
    
    for (k=0; k<n_theta2; k++)
      par.theta2[k+1] = x[k];
    
    nfpBLP(fxp,dat,par,rnd,var,diag,incv,pred);
    logVector(fxp.expDeltaOut, fxp.delta, n_obs);
    matvecmul(var.XZPHIZXZ, fxp.delta, par.theta1, n_theta1, n_obs);
    matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
    for (n=1; n<=n_obs; n++){
      xi0[n] = fxp.delta[n] - var.Xtheta1[n];
    }
    
    for (k=1; k<=n_theta2; k++){
      for (l=1; l<=n_theta2; l++)
        parP.theta2[l] = par.theta2[l];
      parP.theta2[k] = par.theta2[k] + epsilon;
      
      nfpBLP(fxp,dat,parP,rnd,var,diag,incv,pred);
      logVector(fxp.expDeltaOut, fxp.delta, n_obs);
      //matvecmul(var.XZPHIZXZ, fxp.delta, parP.theta1, n_theta1, n_obs);
      //matvecmul(var.X, parP.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
      for (n=1; n<=n_obs; n++){
        xiP[n] = fxp.delta[n] - var.Xtheta1[n];
        dXi[n][k] = (xi0[n] - xiP[n]) / epsilon;
      }
    }
    for (n=1; n<=n_obs; n++){
      for (k=1; k<=n_theta1; k++)
        dXi[n][n_theta2+k] = -var.X[n][k];
    }
    matmul(var.Z_t, dXi, dZe, n_inst, n_obs, n_theta);
    transpose(dZe, dZeTrans, n_inst, n_theta);
    matmul(var.invPHI, dZe, PhiZe, n_inst, n_inst, n_theta);
    matmul(dZeTrans, PhiZe, ZePhiZe, n_theta, n_inst, n_theta);
    for (i=1; i<=n_theta; i++){
      for (j=1; j<=n_theta; j++){
        gsl_matrix_set(gsl_ZePhiZeInv, i-1, j-1, ZePhiZe[i][j]);
      }
    }
    gsl_linalg_cholesky_decomp(gsl_ZePhiZeInv);
    gsl_linalg_cholesky_invert(gsl_ZePhiZeInv);
    
    for (k=1; k<=n_theta2; k++){
      par.sd2[k] = sqrt(gsl_matrix_get(gsl_ZePhiZeInv, k-1, k-1));
      //printf("SE for theta2: %21.14f\n",par.sd2[k]);
    }
    for (k=1; k<=n_theta1; k++)
      par.sd1[k] = sqrt(gsl_matrix_get(gsl_ZePhiZeInv, n_theta2+k-1,n_theta2+k-1));
    
    gsl_matrix_free(gsl_ZePhiZeInv);
    
    free_dvector(xi0, 1, n_obs);
    free_dvector(xiP, 1, n_obs);
    free_dmatrix(dXi, 1, n_obs, 1, n_theta);
    free_dmatrix(dZe, 1, n_inst, 1, n_theta);
    free_dmatrix(dZeTrans, 1, n_theta, 1, n_inst);
    free_dmatrix(PhiZe, 1, n_inst, 1, n_theta);
    free_dmatrix(ZePhiZe, 1, n_theta, 1, n_theta);
    
    free_dvector(parP.theta1, 1, n_theta1);
    free_dvector(parP.theta2, 1, n_theta2);
  }
  else {
    for (k=0; k<n_theta2; k++)
      par.sd2[k] = 0;
    for (k=0; k<n_theta1; k++)
      par.sd1[k] = 0;
    
  }
  toc = gettimeofday_sec();
  
  diag.FC_evals    = KTR_get_number_FC_evals(kc);
  diag.GA_evals    = KTR_get_number_GA_evals(kc);
  diag.H_evals     = KTR_get_number_H_evals(kc);
  diag.HV_evals    = KTR_get_number_HV_evals(kc);
  diag.iters       = KTR_get_number_iters(kc);
  diag.major_iters = KTR_get_number_major_iters(kc);
  diag.minor_iters = KTR_get_number_minor_iters(kc);
  
}

int deltaUpdate
(
 double *expDeltaOut,
 double *expDeltaIn,
 struct DAT  &dat,
 struct PRED &pred
)
{
  int j, t;
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      if (!isinf(pred.ms[t][j]) && !isnan(pred.ms[t][j])){
        pred.ms[t][j] /= n_person;
        expDeltaOut[idx[t][j]] = expDeltaIn[idx[t][j]] * (dat.ms[t][j] / pred.ms[t][j]);
      }
      else{
        printf("t: %d j: %d obs MS: %.14f pred MS: %.14f\n", t, j, dat.ms[t][j], pred.ms[t][j]);
        printf("NaN delta. Aborting blp mapping.\n");
        return -1;
      }
    }
  }
  return 0;
}
