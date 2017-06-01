#include <stdio.h>
#include <stdlib.h>
#include "fxpMap.h"
#include "global.h"
#include "extern.h"
#include "functions.h"
#include "tools.h"
#include <math.h>
#include <float.h>

void nfpBLPcore
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR  &par,
 struct VAR  &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  int n, it = 0;
  
  double
  dist = DBL_MAX,
  distPrev = 1;
  
  diag.exit = 1;
  while (dist>var.tol_blp && it<=1000 && diag.exit>0){
    blpMapMu(expDeltaOut, valuefunOut, expDeltaIn, valuefunIn, expMuIn,
             par, var, data, diag, incv, pred);
    
    dist = 0;
    for (n=1; n<=n_obs; n++){
      dist = fmax(dist, fabs(log(expDeltaIn[n]/expDeltaOut[n])));
      expDeltaIn[n] = expDeltaOut[n];
    }
    it++;
    
    if (diagPrint==1){
      printf("Iter: %5d  Dist: %18.15f  BLP mod: %f\n",it, dist, dist/distPrev);
    }
    
    if (it>2){
      diag.lipschitz += dist / distPrev;
      diag.lipschitzIter++;
    }
    distPrev = dist;
  }
  diag.blpIter += it;
  
  if (diagPrint==1){
    printf("--------------------------------\n");
  }
}

void nfpBLP
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  nfpBLPcore(expDeltaOut, valuefunOut,
             expDeltaIn, valuefunIn, expMuIn,
             par, var, data, diag, incv, pred);
}

void nfpBLP0
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double ***expMuOut,
 double *expDeltaIn,
 double **valuefunIn,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
)
{
  initializeMu(expMuOut, par, rnd, var);
  nfpBLPcore(expDeltaOut, valuefunOut,
             expDeltaIn, valuefunIn, expMuOut,
             par, var, data, diag, incv, pred);
  
}

void blpMapMu
(
 //=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
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
    inclusiveValue(incv, expDeltaIn, expMuIn[i], var);
    transProb(var, incv);
    nfpBellman(valuefunOut[i], valuefunIn[i], var, diag, incv);
    expectedValuefun(var.value0, valuefunOut[i], var.weight, incv);
    choicePath(pred, var, diag);
  }
  deltaUpdate(expDeltaOut, expDeltaIn, data, diag, pred);
}

void blpMap
(
 //=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double ***expMu,
 double *expDeltaIn,
 double **valuefunIn,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
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
  
  //initializeMu_mkl(expMu, par, rnd, var);
  for (i=1; i<=n_person; i++){
    initializeMui(expMu[i], rnd.nu[i], par, var);
    inclusiveValue(incv, expDeltaIn, expMu[i], var);
    transProb(var, incv);
    nfpBellman(valuefunOut[i], valuefunIn[i], var, diag, incv);
    expectedValuefun(var.value0, valuefunOut[i], var.weight, incv);
    choicePath(pred, var, diag);
  }
  deltaUpdate(expDeltaOut, expDeltaIn, data, diag, pred);
}

void deltaUpdate
(
 double *expDeltaOut,
 double *expDeltaIn,
 struct DATA &data,
 struct DIAG &diag,
 struct PRED &pred
 )
{
  int j, t;
  
  for (t=1; t<=n_period; t++){
    for (j=1; j<=n_product; j++){
      pred.ms[t][j] /= n_person;
      expDeltaOut[idx[t][j]] = expDeltaIn[idx[t][j]] * (data.ms[t][j] / pred.ms[t][j]);
    }
  }
  
}
