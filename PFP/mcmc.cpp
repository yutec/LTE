//
//  mcmc.cpp
//  dynblp4LTE
//
//  Created by Yutec Sun on 7/24/15.
//  Copyright (c) 2015 Yutec Sun. All rights reserved.
//
#include <math.h>
#include <float.h>
#include "mcmc.h"
#include "tools.h"
#include "global.h"
#include "extern.h"
#include "fxpMap.h"
#include "functions.h"

void mcmc
(
 int &mcID,
 int &seed,
 int &imh,
 int &vmh,
 int &neighborPast,
 int &n_neighbor,
 int &ptAccpt,
 double &pllh,
 double *jump,
 struct FXP &fxp,
 struct PAR &par_a,
 struct PAR &par_c,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct HIST &hist,
 struct INCV &incv,
 struct POST &post,
 struct PRED &pred,
 gsl_rng *rng,
 FILE *outp
)
{
  int
  k,
  neighborNew;
  
  double
  u,
  llh,
  pllhbackup;
  
  for (k=1; k<=n_theta2; k++){
    par_c.theta2[k] = par_a.theta2[k] + jump[k] * gsl_ran_gaussian(rng, 1);
    while (par_c.theta2[k]<0)
      par_c.theta2[k] = par_a.theta2[k] + jump[k] * gsl_ran_gaussian(rng, 1);
  }

  u = gsl_rng_uniform_pos(rng);
  var.tol_bellman = fmin(1e-4, 1.0/imh);
  var.expMuIn = (ptAccpt==0) ? var.expMu0 : var.expMu1;

  if (diag.numRejectTotal/diag.maxReject>=5){
    diag.exit = -1;
    return;
  }
  else if (diag.numReject>diag.maxReject && n_neighbor<n_history){
    rollbackTest(imh,vmh,n_neighbor,neighborPast,pllh,pllhbackup,
                 fxp,par_a,var,data,diag,hist,incv,pred);
  }
  else{
    blpMapMu(fxp.expDeltaOut, fxp.valuefunOut,
             hist.expDelta[neighborPast], hist.valuefun[neighborPast], var.expMuIn,
             par_a, var, data, diag, incv, pred);
    logQuasiposterior(pllh, fxp.expDeltaOut, par_a, var);
    copyVector(fxp.expDeltaOut, hist.expDelta[neighborPast], n_obs);
    copyMatrix(fxp.valuefunOut, hist.valuefun[neighborPast], n_person, n_grid);
  }
  
  neighborNew = findNeighbor(n_neighbor, hist.param, par_c);
  moveHistory(imh, vmh, n_neighbor, hist.reset);
  var.expMuOut = (ptAccpt==0) ? var.expMu1 : var.expMu0;
  blpMap(hist.expDelta[vmh], hist.valuefun[vmh], var.expMuOut,
         hist.expDelta[neighborNew], hist.valuefun[neighborNew],
         par_c, rnd, var, data, diag, incv, pred);
  logQuasiposterior(llh, hist.expDelta[vmh], par_c, var);
  
  // Accept candidate parameter with prob min(llh-pllh,0)
  pllhbackup = pllh;
  if (log(u) <= fmin(llh-pllh, 0)){
    neighborPast = vmh;
    hist.accept[vmh] = 1;
    copyVector(par_c.theta2, par_a.theta2, n_theta2);
    copyVector(par_c.theta1, par_a.theta1, n_theta1);
    ptAccpt = 1 - ptAccpt;
    pllh = llh;
    diag.numReject = 0;
    diag.numRejectTotal = 0;
    diag.cumulAccept++;
  }
  else{
    hist.accept[vmh] = 0;
    diag.numReject++;
    diag.numRejectTotal++;
  }

  //-------------------------------------------------------------------------
  // Store Markov chain
  //-------------------------------------------------------------------------
  // Store previously accepted & candidate parameters
  copyVector(par_c.theta2, hist.param[vmh], n_theta2);
  copyVector(par_a.theta2, post.theta[imh], n_theta2);
  for (k=1; k<=n_theta1; k++)
    post.theta[imh][n_theta2+k] = par_a.theta1[k];
  
  printMCMC(mcID,vmh,imh,seed,hist.accept[vmh],neighborNew,n_neighbor,
            llh,pllhbackup,par_a,diag,outp);

}

void logQuasiposterior
(
 //=======================================================
 double &llh,
 struct FXP &fxp,
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
  //--------------------------------------------------------------------
  // Uses nested fixed point algorithm
  //--------------------------------------------------------------------
  
  int i, n;
  double obj;
  
  nfpBLP0(fxp.expDeltaOut, fxp.valuefunOut, var.expMu0,
          fxp.expDeltaIn, fxp.valuefunIn,
          par, rnd, var, data, diag, incv, pred);
  
  logVector(fxp.expDeltaOut, var.delta, n_obs);
  matvecmul(var.XZPHIZXZ, var.delta, par.theta1, n_theta1, n_obs);
  matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
  for (n=1; n<=n_obs; n++){
    var.e[n] = var.delta[n] - var.Xtheta1[n];
  }
  
  matvecmul(var.Z_t, var.e, var.Ze, n_inst, n_obs);            // Z'e
  matvecmul(var.invPHI, var.Ze, var.PHIZe, n_inst, n_inst);   // PHIZe = inv(Z'Z)(Z'e)
  
  obj = 0;
  for (i=1; i<=n_inst; i++)
    obj -= var.Ze[i] * var.PHIZe[i];
  
  llh = obj * 0.5;
  
}

void logQuasiposterior
(
 //=======================================================
 double &llh,
 double *expDeltaIn,
 struct PAR &par,
 struct VAR &var
//=======================================================
)
{
  //--------------------------------------------------------------------
  // Uses pseudo fixed point algorithm
  //--------------------------------------------------------------------
  
  int i, n;
  double obj;
  
  logVector(expDeltaIn, var.delta, n_obs);
  //logVector_mkl(expDeltaIn, var.delta, n_obs);
  matvecmul(var.XZPHIZXZ, var.delta, par.theta1, n_theta1, n_obs);
  matvecmul(var.X, par.theta1, var.Xtheta1, n_obs, n_theta1);  // X*theta1
  for (n=1; n<=n_obs; n++){
    var.e[n] = var.delta[n] - var.Xtheta1[n];
  }
  
  matvecmul(var.Z_t, var.e, var.Ze, n_inst, n_obs);           // Z'e
  matvecmul(var.invPHI, var.Ze, var.PHIZe, n_inst, n_inst);   // PHIZe = inv(Z'Z)(Z'e)
  
  obj = 0;
  for (i=1; i<=n_inst; i++)
    obj -= var.Ze[i] * var.PHIZe[i];
  
  llh = obj * 0.5;
}

void rollbackTest
(
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &neighborPast,
 double &pllh,
 double &pllhbackup,
 struct FXP &fxp,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct HIST &hist,
 struct INCV &incv,
 struct PRED &pred
)
{
  int tolBellman;
  
  tolBellman = var.tol_bellman;
  var.tol_bellman = tol_blm_pfp;
  nfpBLP(fxp.expDeltaOut, fxp.valuefunOut,
         hist.expDelta[neighborPast], hist.valuefun[neighborPast], var.expMuIn,
         par, var, data, diag, incv, pred);
  logQuasiposterior(pllh, fxp.expDeltaOut, par, var);
  copyVector(fxp.expDeltaOut, hist.expDelta[neighborPast], n_obs);
  copyMatrix(fxp.valuefunOut, hist.valuefun[neighborPast], n_person, n_grid);
  
  if (fabs(pllh-pllhbackup)>0.1){
    //hist.reset = imh;
    if (vmh>neighborPast){
      resetHistParam(neighborPast+1,vmh,hist);
    }
    else{
      resetHistParam(1,vmh,hist);
      resetHistParam(neighborPast,n_neighbor,hist);
    }
    //n_neighbor = neighborPast;
    vmh = neighborPast;
  }
  diag.maxReject++;
  diag.numReject = 0;
  var.tol_bellman = tolBellman;
  //imh -= diag.numRejectTotal;
}


