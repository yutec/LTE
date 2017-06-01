#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "extern.h"
#include "tools.h"
#include "functions.h"
#include "nfpIteration.h"
#include "mpec.h"
//#include "knitro.h"

int main(void) {
  //======================================================================
  // Declare utility variables
  //======================================================================
  int
  i, j, k, t,
  seed,
  mcID, // index for Monte Carlo experiments
  nStatus;
  
  double
  *paramIn,         // input parameter vector
  *paramBound,
  ***paramSeed,
  toc;
  
  struct USERSTRUCT us;
  
  char
  directoryInput[1024],
  filenameData[1024],
  filenamePar[1024],
  filenameRand[1024],
  filenameValfun[1024],
  filenameSeed[1024],
  filenameEstim[1024],
  paramNameDummy[512],
  parName[n_theta1+n_theta2][512];
  
  FILE *inpPar, *inpData, *inpRand, *outp, *outp_stat;
  
  //======================================================================
  // Memory allocation
  //======================================================================
  paramBound  = dvector(n_theta2);
  paramIn     = dvector(n_theta2);
  paramSeed   = darray3d(n_mc, n_seed, n_theta2);
  
  allocateMemory(us.dat,us.fxp,us.par,us.rnd,us.var,us.incv,us.pred);
  
  //============================================================================
  // Read in the starting parameter values
  //============================================================================
  strcpy(directoryInput,"./data/");
  strcpy(filenameData, directoryInput);
  strcpy(filenamePar, directoryInput);
  strcpy(filenameRand, directoryInput);
  strcpy(filenameSeed, directoryInput);
  strcpy(filenameValfun, directoryInput);
  strcpy(filenameEstim, "estimResult.txt");
  strcat(filenameData, "dat_dynblp4.txt");
  strcat(filenamePar, "paramInput.txt");
  strcat(filenameRand, "rand.txt");
  strcat(filenameSeed, "seed.txt");
  strcat(filenameValfun, "valfun.txt");
  
  inpPar = fopen(filenamePar,"r");
  for (i=0; i<n_theta2; i++){
    nStatus = fscanf(inpPar, "%lf %s", &paramIn[i], paramNameDummy);
    printf("%12s = %6.3f\n", paramNameDummy, paramIn[i]);
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
  
  //printf("Enter seed number for optimization (1-5): ");
  //nStatus = scanf("%d", &seed);
  
  //============================================================================
  // Open file streams
  //============================================================================
  // Erase pre-existing output files
  outp = fopen("outputMPEC.txt", "w");
  fprintf(outp, "%6s %6s ","mcID","Seed");
  for (k=0; k<n_theta1+n_theta2; k++)
    fprintf(outp, "%20s ", parName[k]);
  fprintf(outp, "%16s %14s\n","Obj","CPUtime");
  fclose(outp);
  
  outp_stat = fopen(filenameEstim, "w");
  fprintf(outp_stat,"%6s %6s %6s ","MC","Seed","Exit");
  for (k=1; k<=n_theta2; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  for (k=n_theta2+1; k<=n_theta; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  fprintf(outp_stat,"%14s %8s %8s %8s %8s %12s %14s\n",
          "Obj","Iter","F_eval","G_eval","H_eval","BellmanIter","CPUtime");
  fclose(outp_stat);
  
  // Open input file streams
  inpData = fopen(filenameData, "r");
  
  // Read consumer simulation draw
  inpRand = fopen(filenameRand, "r");
  
  // Index for obs
  for (t=0; t<n_period; t++)
    for (j=0; j<n_product; j++)
      idx[t][j] = n_product*t + j;
  
  paramBound[0] = 5;
  paramBound[1] = 5;
  paramBound[2] = 2.5;
  
  readSeed(paramSeed, filenameSeed);
  
  //======================================================================
  // KNITRO setup
  //======================================================================
  KTR_context *kc;
  int n, m, nnzJ, nnzH, objGoal, objType, *cType;
  double obj, *x, *lambda;
  double *xLoBnds, *xUpBnds, *xInitial, *cLoBnds, *cUpBnds;
  
  n = dimOptim;
  m = dimConst;
  
  // Number of nonzero elements in sparse constraint Jacobian & Hessian
  //nnzJ = n_obs*(n_theta+n_obs+dimValfun) + dimValfun*(n_theta+n_obs+n_grid) + n_inst*(n_obs+1);
  nnzH = (n_theta+n_obs+dimValfun)*(n_theta+n_obs+dimValfun+1)/2 + n_inst*(n_inst+1)/2;
  
  objType = KTR_OBJTYPE_GENERAL;
  objGoal = KTR_OBJGOAL_MINIMIZE;
  
  //======================================================================
  // Monte Carlo loop
  //======================================================================
  // Skip dat prior to mcID
  for (mcID=0; mcID<mc_init; mcID++){
    readData(mcID, inpData, inpRand, filenameValfun, us.dat, us.rnd, us.var);
  }
  
  for (mcID=mc_init; mcID<mc_last; mcID++){
    readData(mcID, inpData, inpRand, filenameValfun, us.dat, us.rnd, us.var);
    
    for (seed=0; seed<n_seed; seed++){
      allocateKnitro(m,n,nnzH,&cType,&xLoBnds,&xUpBnds,&xInitial,&x,&cLoBnds,&cUpBnds,&lambda,us);
      initialize(mcID,seed,xInitial,paramSeed[mcID][seed],us);
      initializeKnitro(m,n,nnzJ,nnzH,cType,xInitial,xLoBnds,xUpBnds,cLoBnds,cUpBnds,paramBound,us);

      kc = KTR_new();
      nStatus = KTR_load_param_file(kc, "knitroopt.txt");
      
      if (KTR_set_func_callback(kc, (KTR_callback * const) callbackEvalFC) != 0)
        exit(-1);
      if (KTR_set_grad_callback(kc, (KTR_callback * const) callbackEvalGA) != 0)
        exit(-1);
      if (KTR_set_hess_callback(kc, (KTR_callback * const) callbackEvalH) != 0)
        exit(-1);
      
      nStatus = KTR_init_problem(kc, n, objGoal, objType,
                                 xLoBnds, xUpBnds,
                                 m, cType, cLoBnds, cUpBnds,
                                 nnzJ, us.jacIndexVarsInt, us.jacIndexConsInt,
                                 nnzH, us.hessIndexRows, us.hessIndexCols, xInitial, NULL);
      
      nStatus = KTR_solve(kc, x, lambda, 0, &obj,
                          NULL, NULL, NULL, NULL, NULL, &us);
      stdError(nStatus,x,toc,us.dat,us.fxp,us.par,us.rnd,us.var,us.diag,us.incv,us.pred,kc);
      printOutput(filenameEstim,mcID,seed,nStatus,obj,toc,x,us);
      KTR_free(&kc);
      freeKnitro(&cType,&xLoBnds,&xUpBnds,&xInitial,&x,&cLoBnds,&cUpBnds,&lambda,us);
    } // end of seed loop
  } // end of mcID loop
  
  //======================================================================
  // Free up Memory
  //======================================================================
  fclose(inpData);
  fclose(inpRand);
  
  free_dvector(paramBound);
  free_dvector(paramIn);
  free_darray3d(paramSeed);
  
  releaseMemory(us.dat,us.fxp,us.par,us.rnd,us.var,us.incv,us.pred);
  
  printf("Monte Carlo experiment complete.\n");
  
  return 0;
}
