#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "extern.h"
#include "tools.h"
#include "functions.h"
#include "nfpIteration.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

void computeExpectedValuefun
(
 //=======================================================
 // Computes i-specific expected value E[V_i(I_{t+1})|I_t] of
 // no purchase option given observed states I_t at time t.
 double *expectValuefun,
 double *valuefunIn,
 struct INCV &incv
//=======================================================
)
{
  int s1, t;
  double
  condMean,     // conditional mean of next-period state given observed states
  value0,       // unweighted expected value function given observed states
  sumWeight,   // sum of transition prob weight
  weight[n_grid];
  
  for (t=0; t<n_period; t++){
    // Compute conditional prob of state wGrid[s1] given state w_i[m][t]
    condMean = incv.lambda[0] + incv.lambda[1] * incv.w[t];
    sumWeight = 0;
    for (s1=0; s1<n_grid; s1++){
      weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += weight[s1];
    }
    
    // Compute expected value function
    value0 = 0;
    for (s1=0; s1<n_grid; s1++){
      value0 += weight[s1] * valuefunIn[s1];
    }
    expectValuefun[t] = value0 / sumWeight;
  } // t loop
}

void computeExpectedValuefun
(
 //=======================================================
 adouble *expectValuefun,
 adouble *valuefunIn,
 struct INCV_AD &incv
//=======================================================
)
{
  int s1, t;
  adouble
  condMean,     // conditional mean of next-period state given observed states
  value0,       // unweighted expected value function given observed states
  sumWeight,   // sum of transition prob weight
  weight[n_grid];
  
  for (t=0; t<n_period; t++){
    // Compute conditional prob of state wGrid[s1] given state w_i[m][t]
    condMean = incv.lambda[0] + incv.lambda[1] * incv.w[t];
    sumWeight = 0;
    for (s1=0; s1<n_grid; s1++){
      weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += weight[s1];
    }
    
    // Compute expected value function
    value0 = 0;
    for (s1=0; s1<n_grid; s1++){
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
  while (dist > var.tol_bellman && iter <= 5000){
    bellmanOperator(valuefunOut, valuefunIn, var, incv);
    dist = 0;
    for (s0=0; s0<n_grid; s0++){
      dist = fmax(dist, fabs(valuefunIn[s0] - valuefunOut[s0]));
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
  while (dist > var.tol_blp && it <= 5000){
    blpMap(fxp.expDeltaOut, fxp.valuefunOut, fxp.expDeltaIn, fxp.valuefunIn, fxp.expMu,
           dat, par, var, diag, incv, pred);
    
    dist = 0;
    for (n=0; n<n_obs; n++){
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
 // Returns consumer i-specific value function (valuefunOut)
 // after applying the Bellman operator once to value function
 // (valuefunIn), given states (wGrid).
 double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct INCV &incv
//=======================================================
)
{
  int s0, s1;
  double value0;
  
  for (s0=0; s0<n_grid; s0++){
    // value0: expected value function of delaying purchase
    value0 = 0;
    for (s1=0; s1<n_grid; s1++){
      value0 += valuefunIn[s1] * var.weightMat[s0][s1];
    }
    value0 *= beta;
    valuefunOut[s0] = log(1 + exp(incv.wGrid[s0] - value0)) + value0;
  } // end of s0 loop
  
}

void bellmanOperator
(//=======================================================
 adouble *valuefunOut,
 adouble *valuefunIn,
 struct VAR_AD &var,
 struct INCV_AD &incv
//=======================================================
)
{
  int s0, s1;
  adouble value0;
  
  for (s0=0; s0<n_grid; s0++){
    // value0: expected value function of delaying purchase
    value0 = 0;
    for (s1=0; s1<n_grid; s1++){
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
  
  for (t=0; t<n_period; t++) {
    for (j=0; j<n_product; j++) {
      pred.ms[t][j] = 0;
    }
  }
  
  for (i=0; i<n_person; i++){
    computeInclusiveValue(incv, expDeltaIn, expMuIn[i], var);
    transProb(var, incv);
    nfpBellman(valuefunOut[i], valuefunIn[i], var, diag, incv);
    computeExpectedValuefun(var.value0, valuefunOut[i], incv);
    computeChoicePath(pred, var);
  }
  
  for (t=0; t<n_period; t++){
    for (j=0; j<n_product; j++){
      pred.ms[t][j] /= n_person;
      expDeltaOut[idx[t][j]] = expDeltaIn[idx[t][j]] * (dat.ms[t][j] / pred.ms[t][j]);
      if (isinf(expDeltaOut[idx[t][j]]) || isnan(expDeltaOut[idx[t][j]])){
        printf("NaN delta in blpMap routine.\n");
        for (int k=0; k<n_theta2; k++)
          printf("Theta[%d]: %21.16f\n",k,par.theta2[k]);
      }
    }
  }
}

void computeInclusiveValue
(//=======================================================
 struct INCV_AD &incv,
 adouble *expDeltaIn,
 adouble **expMu,
 struct VAR_AD &var
//=======================================================
)
{
  //--------------------------------------------------------------------
  int j, n, s, t;
  adouble
  max_w = -DBL_MAX,
  min_w = DBL_MAX,
  det,
  expInclusiveValue;
  
  // Compute logit inclusive values w[t] for consumer i
  for (t=0; t<n_period; t++){
    expInclusiveValue = 0;
    for (j=0; j<n_product; j++){
      var.num[t][j] = expDeltaIn[idx[t][j]] * expMu[t][j];
      expInclusiveValue += var.num[t][j];
    }
    incv.w[t] = log(expInclusiveValue);
    condassign(max_w, incv.w[t]-max_w, incv.w[t], max_w);
    condassign(min_w, min_w-incv.w[t], incv.w[t], min_w);
    
    //max_w = (max_w < incv.w[t]) ? incv.w[t] : max_w;
    //min_w = (min_w > incv.w[t]) ? incv.w[t] : min_w;
  }
  
  for (t=0; t<n_period-1; t++){
    incv.y[t] = incv.w[t+1];
    //incv.w_lag[t][0].set_value(1.0);
    incv.w_lag[t][0] = 1;
    incv.w_lag[t][1] = incv.w[t];
  }
  
  // Estimate lambda in transition process of inclusive values
  transpose(incv.w_lag, incv.w_lag_t, n_period-1, 2);
  matmul(incv.w_lag_t, incv.w_lag, incv.ww, 2, n_period-1, 2);
  matvecmul(incv.w_lag_t, incv.y, incv.wy, 2, n_period-1);
  
  det = fabs(incv.ww[0][0]*incv.ww[1][1] - incv.ww[0][1]*incv.ww[1][0]);
  incv.ww_inv[0][0] = incv.ww[1][1] / det;
  incv.ww_inv[1][1] = incv.ww[0][0] / det;
  incv.ww_inv[0][1] = incv.ww_inv[1][0] = -incv.ww[0][1] / det;
  matvecmul(incv.ww_inv, incv.wy, incv.lambda, 2, 2); // make sure lambda works
  matvecmul(incv.w_lag, incv.lambda, incv.wlambda, n_period-1, 2);
  
  // Estimate variance of psi
  incv.varPsi = 0;
  for (n=0; n<n_period-1; n++){
    incv.varPsi += (incv.y[n] - incv.wlambda[n]) * (incv.y[n] - incv.wlambda[n]);
  }
  incv.varPsi /= n_period - 2;
  
  // Set up grid space
  max_w = max_w + 2.575829*sqrt(incv.varPsi);
  min_w = min_w - 2.575829*sqrt(incv.varPsi);
  for (s=0; s<n_grid; s++){
    incv.wGrid[s] = min_w + s * (max_w - min_w)/(n_grid-1);
  }
  
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
  //--------------------------------------------------------------------
  int j, n, s, t;
  double
  det,
  expInclusiveValue,            // exponentiated logit inclusive value
  max_w = -DBL_MAX,
  min_w = DBL_MAX;
  
  // Compute logit inclusive values w[t] for consumer i
  for (t=0; t<n_period; t++){
    expInclusiveValue = 0;
    for (j=0; j<n_product; j++){
      var.num[t][j] = expDeltaIn[idx[t][j]] * expMu[t][j];
      expInclusiveValue += var.num[t][j];
    }
    incv.w[t] = log(expInclusiveValue);
    max_w = (max_w < incv.w[t]) ? incv.w[t] : max_w;
    min_w = (min_w > incv.w[t]) ? incv.w[t] : min_w;
  }
  
  for (t=0; t<n_period-1; t++){
    incv.y[t] = incv.w[t+1];
    incv.w_lag[t][0] = 1;
    incv.w_lag[t][1] = incv.w[t];
  }
  
  // Estimate lambda in transition process of inclusive values
  transpose(incv.w_lag, incv.w_lag_t, n_period-1, 2);
  matmul(incv.w_lag_t, incv.w_lag, incv.ww, 2, n_period-1, 2);
  matvecmul(incv.w_lag_t, incv.y, incv.wy, 2, n_period-1);
  
  det = fabs(incv.ww[0][0]*incv.ww[1][1] - incv.ww[0][1]*incv.ww[1][0]);
  incv.ww_inv[0][0] = incv.ww[1][1] / det;
  incv.ww_inv[1][1] = incv.ww[0][0] / det;
  incv.ww_inv[0][1] = incv.ww_inv[1][0] = -incv.ww[0][1] / det;
  matvecmul(incv.ww_inv, incv.wy, incv.lambda, 2, 2); // make sure lambda works
  matvecmul(incv.w_lag, incv.lambda, incv.wlambda, n_period-1, 2);
  
  // Estimate variance of psi
  incv.varPsi = 0;
  for (n=0; n<n_period-1; n++){
    incv.varPsi += (incv.y[n] - incv.wlambda[n]) * (incv.y[n] - incv.wlambda[n]);
  }
  incv.varPsi /= n_period - 2;
  
  // Set up grid space
  max_w = max_w + 2.575829*sqrt(incv.varPsi);
  min_w = min_w - 2.575829*sqrt(incv.varPsi);
  for (s=0; s<n_grid; s++){
    incv.wGrid[s] = min_w + s * (max_w - min_w)/(n_grid-1);
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
  double condMean, sumWeight, weight[n_grid];
  
  for (s0=0; s0<n_grid; s0++){
    condMean = incv.lambda[0] + incv.lambda[1] * incv.wGrid[s0];
    sumWeight = 0;
    for (s1=0; s1<n_grid; s1++){
      weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += weight[s1];
    }
    for (s1=0; s1<n_grid; s1++){
      var.weightMat[s0][s1] = weight[s1]/sumWeight;
    }
  }
}

void transProb_AD
(
 //=======================================================
 struct VAR_AD &var,
 struct INCV_AD &incv
//=======================================================
)
{
  int s0, s1;
  adouble condMean, sumWeight, weight[n_grid];
  
  for (s0=0; s0<n_grid; s0++){
    condMean = incv.lambda[0] + incv.lambda[1] * incv.wGrid[s0];
    sumWeight = 0;
    for (s1=0; s1<n_grid; s1++){
      weight[s1] = exp(-0.5*(incv.wGrid[s1]-condMean)*(incv.wGrid[s1]-condMean)/incv.varPsi);
      sumWeight += weight[s1];
    }
    for (s1=0; s1<n_grid; s1++){
      var.weightMat[s0][s1] = weight[s1]/sumWeight;
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
  
  // Initialize survivalRateard (var.survivalRate) at t=0 for all i
  pred.survivalRate[0] = 1;
  
  // Compute predicted market share given deltaIn
  for (t=0; t<n_period; t++){
    var.num[t][n_product] = exp(beta * var.value0[t]);
    den = 0;
    for (j=0; j<n_product+1; j++){
      den += var.num[t][j];
    }
    
    // Add choice prob * (1-survivalRateard) to predicted market share
    for (j=0; j<n_product; j++){
      pred.ms[t][j] += var.num[t][j] / den * pred.survivalRate[t];
    }
    
    // Update the next-period survivalRateard
    if (t < n_period-1){
      pred.survivalRate[t+1] = var.num[t][n_product] / den * pred.survivalRate[t];
    }
    
  } // t loop
}

void computeChoicePath_AD
(
 //=======================================================
 struct PRED_AD &pred,
 struct VAR_AD &var
//=======================================================
)
{
  int j, t;
  adouble den;
  
  // Initialize survivalRateard (var.survivalRate) at t=0 for all i
  pred.survivalRate[0] = 1;
  
  // Compute predicted market share given deltaIn
  for (t=0; t<n_period; t++){
    var.num[t][n_product] = exp(beta * var.value0[t]);
    den = 0;
    for (j=0; j<n_product+1; j++){
      den += var.num[t][j];
    }
    
    // Add choice prob * (1-survivalRateard) to predicted market share
    /*
    if (isnan(den) || isinf(den)){
      printf("NaN den\n");
      printf("Numerator[n_product] %.8f %.8f\n",
             var.num[t][n_product].value(), var.value0[t].value());
      printf("Beta: %.8f expValfun: %.8f\n", beta, var.value0[t].value());
      abort();
    }
    */
    
    for (j=0; j<n_product; j++){
      pred.ms[t][j] += var.num[t][j] / den * pred.survivalRate[t];
    }
    
    // Update the next-period survivalRateard
    if (t < n_period-1){
      pred.survivalRate[t+1] = var.num[t][n_product] / den * pred.survivalRate[t];
    }
  } // t loop
}

void stdError
(
 int    &nStatus,
 double *x,
 double &toc,
 struct DAT &dat,
 struct FXP &fxp,
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
  
  struct PAR  parP;
  
  if (nStatus > -200){
    gsl_ZePhiZeInv = gsl_matrix_alloc(n_theta, n_theta);
    
    xi0 = dvector(n_obs);
    xiP = dvector(n_obs);
    dXi = dmatrix(n_obs, n_theta);
    dZe = dmatrix(n_inst, n_theta);
    dZeTrans   = dmatrix(n_theta, n_inst);
    PhiZe      = dmatrix(n_inst, n_theta);
    ZePhiZe    = dmatrix(n_theta, n_theta);
    
    parP.theta1 = dvector(n_theta1);
    parP.theta2 = dvector(n_theta2);
    
    var.tol_blp = tol_blp;
    var.tol_bellman = tol_blm;
    for (k=0; k<n_theta2; k++)
      par.theta2[k] = x[k];
    for (k=0; k<n_theta1; k++)
      par.theta1[k] = x[n_theta2+k];
    for (n=0; n<n_obs; n++)
      xi0[n] = x[n_theta+n];
    
    // To minimize time for NFP
    matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);
    for (n=0; n<n_obs; n++){
      fxp.expDeltaIn[n] = exp(var.Xtheta1[n] + xi0[n]);
    }
    for (i=0; i<n_person; i++){
      for (int s0=0; s0<n_grid; s0++){
        fxp.valuefunIn[i][s0] = x[n_theta+n_obs+n_grid*i+s0];
      }
    }
    
    nfpBLP(fxp,dat,par,rnd,var,diag,incv,pred);
    logVector(fxp.expDeltaOut, fxp.delta, n_obs);
    matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
    for (n=0; n<n_obs; n++){
      xi0[n] = fxp.delta[n] - var.Xtheta1[n];
      //if (x[n_theta+n] != xi0[n])
      //  printf("NFP xi not match MPEC xi: %18.14f %18.14f %5d\n",xi0[n],x[n_theta+n],n);
    }
    
    for (k=0; k<n_theta2; k++){
      for (l=0; l<n_theta2; l++)
        parP.theta2[l] = par.theta2[l];
      parP.theta2[k] = par.theta2[k] + epsilon;
      
      nfpBLP(fxp,dat,parP,rnd,var,diag,incv,pred);
      logVector(fxp.expDeltaOut, fxp.delta, n_obs);
      for (n=0; n<n_obs; n++){
        xiP[n] = fxp.delta[n] - var.Xtheta1[n];
        dXi[n][k] = (xi0[n] - xiP[n]) / epsilon;
      }
    }
    for (n=0; n<n_obs; n++){
      for (k=0; k<n_theta1; k++)
        dXi[n][n_theta2+k] = -var.X[n][k];
    }
    matmul(var.Z_t, dXi, dZe, n_inst, n_obs, n_theta);
    transpose(dZe, dZeTrans, n_inst, n_theta);
    matmul(var.invPHI, dZe, PhiZe, n_inst, n_inst, n_theta);
    matmul(dZeTrans, PhiZe, ZePhiZe, n_theta, n_inst, n_theta);
    for (i=0; i<n_theta; i++){
      for (j=0; j<n_theta; j++){
        gsl_matrix_set(gsl_ZePhiZeInv, i, j, ZePhiZe[i][j]);
      }
    }
    gsl_linalg_cholesky_decomp(gsl_ZePhiZeInv);
    gsl_linalg_cholesky_invert(gsl_ZePhiZeInv);
    
    for (k=0; k<n_theta2; k++)
      par.sd2[k] = sqrt(gsl_matrix_get(gsl_ZePhiZeInv, k, k));
    for (k=0; k<n_theta1; k++)
      par.sd1[k] = sqrt(gsl_matrix_get(gsl_ZePhiZeInv, n_theta2+k,n_theta2+k));
    
    gsl_matrix_free(gsl_ZePhiZeInv);
    
    free_dvector(xi0);
    free_dvector(xiP);
    free_dmatrix(dXi);
    free_dmatrix(dZe);
    free_dmatrix(dZeTrans);
    free_dmatrix(PhiZe);
    free_dmatrix(ZePhiZe);
    
    free_dvector(parP.theta1);
    free_dvector(parP.theta2);
    
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