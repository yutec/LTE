/*******************************************************************************
 
 This program estimates a simple dynamic adoption model with infinite time horizon.
 
 Estimation procedure options:
 BLP fixed point procedure
 1. Pseudo fixed point
 2. Nested fixed point
 
 Bellman fixed point procedure
 1. Pseudo fixed point
 2. Nested fixed point
 
 First Version: May 4, 2012
 Updated: September 24, 2012
 
 Created by Yutec Sun on 2012-08-20.
 
 *******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
//#include <stddef.h>
#include <string.h>
//#include <unistd.h>
//#include <math.h>
//#include <time.h>
//#include <float.h>
//#include <sys/time.h>
#include "extern.h"
#include "global.h"
#include "tools.h"
#include "functions.h"
#include "nfpIteration.h"
#include "knitro.h"

int main(void)
{
  //======================================================================
  // Declare utility variables
  //======================================================================
  int
  i, j, k, t,
  seed,
  mcID, // index for Monte Carlo experiments
  nStatus;
  
  double
  *paramTrue,         // input parameter vector
  *paramBound,
  ***paramSeed,
  toc;
  
  struct USERSTRUCT us;
  
  char
  directoryInput[1024],
  filenameData[1024],
  filenameRand[1024],
  filenameParam[1024],
  filenameSeed[1024],
  filenameEstim[1024],
  paramNameDummy[512],
  parName[n_theta1+n_theta2][512];
  
  FILE *inpPar, *inpData, *inpRand, *outp, *outp_stat;
  
  //======================================================================
  // Memory allocation
  //======================================================================
  paramBound = dvector(1, n_theta2);
  paramTrue  = dvector(1, n_theta2);
  paramSeed  = darray3d(1, n_mc, 1, n_seed, 1, n_theta2);
  
  allocateMemory(us.dat,us.fxp,us.par,us.rnd,us.var,us.diag,us.incv,us.pred);
  
  //============================================================================
  // Read in the starting parameter values
  //============================================================================
  strcpy(directoryInput,"./data/");
  strcpy(filenameData, directoryInput);
  strcpy(filenameRand, directoryInput);
  strcpy(filenameParam, directoryInput);
  strcpy(filenameSeed, directoryInput);
  strcpy(filenameEstim, "estimResult.txt");
  strcat(filenameData, "dat_dynblp4.txt");
  strcat(filenameRand, "rand.txt");
  strcat(filenameParam, "paramInput.txt");
  strcat(filenameSeed, "seed.txt");
  
  inpPar = fopen(filenameParam,"r");
  for (i=1; i<=n_theta2; i++){
    fscanf(inpPar, "%lf %s", &paramTrue[i], paramNameDummy);
    //printf("%12s = %6.3f\n", paramNameDummy, paramTrue[i]);
  }
  fclose(inpPar);
  
  strcpy(parName[0], "sdBeta1");
  strcpy(parName[1], "sdBeta2");
  strcpy(parName[2], "sdAlpha");
  strcpy(parName[3], "beta0");
  strcpy(parName[4], "beta1");
  strcpy(parName[5], "beta2");
  strcpy(parName[6], "beta3");
  strcpy(parName[7], "alpha");
  
  //printf("Enter seed number for optimization: ");
  //scanf("%d", &seed);
  
  //============================================================================
  // Open file streams
  //============================================================================
  // Erase pre-existing output files
  outp = fopen("outputNFP.txt", "w");
  fprintf(outp, "%6s %6s ","mcID","Seed");
  for (k=1; k<=n_theta1+n_theta2; k++)
    fprintf(outp, "%20s ", parName[k-1]);
  fprintf(outp, "%21s %14s %8s %14s\n","LogQpost","BellmanIter","Modulus","CPUtime");
  fclose(outp);
  
  outp_stat = fopen(filenameEstim, "w");
  fprintf(outp_stat,"%6s %6s %6s ","MC","Seed","Exit");
  for (k=1; k<=n_theta2; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  for (k=n_theta2+1; k<=n_theta; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  fprintf(outp_stat,"%14s %8s %8s %8s %12s %8s %14s\n",
          "Obj","Iter","F_eval","BLPiter","BellmanIter","Modulus","CPUtime");
  fclose(outp_stat);
  
  // Open input file streams
  inpData = fopen(filenameData, "r");
  
  // Read consumer simulation draw
  inpRand = fopen(filenameRand, "r");
  
  // Index for obs
  for (t=1; t<=n_period; t++)
    for (j=1; j<=n_product; j++)
      idx[t][j] = n_product*(t-1) + j;
  
  paramBound[1] = 5;
  paramBound[2] = 5;
  paramBound[3] = 2.5;
  
  readSeed(paramSeed, filenameSeed);
  
  //======================================================================
  // KNITRO setup
  //======================================================================
  KTR_context *kc;
  int n, m, nnzJ, nnzH, objGoal, objType;
  int *cType;
  int *jacIndexVars, *jacIndexCons;
  double obj, *x, *lambda;
  double *xLoBnds, *xUpBnds, *xInitial, *cLoBnds, *cUpBnds;
  
  n = n_theta2;
  m = 0;
  nnzJ = n*m;
  nnzH = 0;
  x = (double *) malloc(n*sizeof(double));
  lambda = (double *) malloc((m+n)*sizeof(double));
  
  xLoBnds = (double *) malloc(n*sizeof(double));
  xUpBnds = (double *) malloc(n*sizeof(double));
  xInitial = (double *) malloc(n*sizeof(double));
  cType = (int *) malloc(m*sizeof(int));
  cLoBnds = (double *) malloc(m*sizeof(double));
  cUpBnds = (double *) malloc(m*sizeof(double));
  jacIndexVars = (int *) malloc(nnzJ*sizeof(int));
  jacIndexCons = (int *) malloc(nnzJ*sizeof(int));
  
  objType = KTR_OBJTYPE_GENERAL;
  objGoal = KTR_OBJGOAL_MINIMIZE;
  
  for (i=0; i<n; i++){
    xLoBnds[i] = 0.0;
    xUpBnds[i] = paramBound[i+1];
  }
  for (j=0; j<m; j++){
    cType[j] = KTR_CONTYPE_GENERAL;
    cLoBnds[j] = 0.0;
    cUpBnds[j] = KTR_INFBOUND;
  }
  
  k = 0;
  for (i=0; i<n; i++){
    for (j=0; j<m; j++){
      jacIndexCons[k] = j;
      jacIndexVars[k] = i;
      k++;
    }
  }
  
  //======================================================================
  // Monte Carlo loop
  //======================================================================
  // Skip dat prior to mcID
  for (mcID=1; mcID<=mc_init-1; mcID++){
    readData(mcID, inpData, inpRand, us.dat, us.rnd, us.var);
  }
  
  for (mcID=mc_init; mcID<=mc_last; mcID++){
    readData(mcID, inpData, inpRand, us.dat, us.rnd, us.var);
    
    for (seed=1; seed<=n_seed; seed++){
      
      kc = KTR_new();
      nStatus = KTR_load_param_file(kc, "knitroopt.txt");
      if (KTR_set_func_callback(kc, (KTR_callback * const) callback) != 0)
        exit(-1);
      initialize(mcID, seed, xInitial, paramSeed[mcID][seed], us);
      nStatus = KTR_init_problem(kc, n, objGoal, objType,
                                 xLoBnds, xUpBnds,
                                 NULL, NULL, NULL, NULL,
                                 nnzJ, jacIndexVars, jacIndexCons,
                                 nnzH, NULL, NULL, xInitial, NULL);
      
      us.tic = gettimeofday_sec();
      nStatus = KTR_solve(kc, x, lambda, 0, &obj,
                          NULL, NULL, NULL, NULL, NULL, &us);
      stdError(nStatus,x,toc,us.fxp,us.dat,us.par,us.rnd,us.var,us.diag,us.incv,us.pred,kc);
      printOutput(filenameEstim,mcID,seed,nStatus,obj,toc,us);
      KTR_free(&kc);
      
    }
    
  } // end of mcID loop
  
  free(x);
  free(lambda);
  free(xLoBnds);
  free(xUpBnds);
  free(xInitial);
  free(cType);
  free(cLoBnds);
  free(cUpBnds);
  free(jacIndexVars);
  free(jacIndexCons);
  
  fclose(inpData);
  fclose(inpRand);
  
  //======================================================================
  // Free up Memory
  //======================================================================
  free_dvector(paramBound,  1, n_theta2);
  free_dvector(paramTrue,   1, n_theta2);
  free_darray3d(paramSeed,  1, n_mc, 1, n_seed, 1, n_theta2);
  
  releaseMemory(us.dat,us.fxp,us.par,us.rnd,us.var,us.diag,us.incv,us.pred);
  
  printf("Monte Carlo experiment complete.\n");
  return 0;
}



