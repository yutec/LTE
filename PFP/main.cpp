/*******************************************************************************
 This program runs the second & third Monte Carlo experiments of dynamic BLP
 demand estimation using Pseudo fixed point algorithm (Sun & Ishihara, 2017).
 
 Compiling and running the program requires GSL (GNU Scientific Libraray) package,
 which can be obtained at "http://www.gnu.org/software/gsl/". Specify the path to
 the headers and the link to the GSL library in the Makefile.
 
 The Makefile compiles the codes into a command-line console execution file "a.out"
 in the Linux & Mac OS. The Makefile uses the GCC compiler while standard
 C/C++ compilers should work with no problem.
 
 First Version: May 4, 2012
 Updated: June 30, 2017
 *******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include "tools.h"
#include "global.h"
#include "extern.h"
#include "functions.h"
#include "fxpMap.h"
#include "mcmc.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <RInside.h>

int main(int argc, char *argv[])
{
  RInside R(argc, argv);
  std::string cmdSetup =         // load library & module
  "library(coda);"
  "source('heidelberger.R')";
  std::string reset =
  "rm(mcmcDraws);"
  "if (exists('h')){rm(h)};"
  "if (exists('thin')){rm(thin)};"
  "gc()";
  
  R.parseEvalQ(cmdSetup);        // eval command, no return
  
  //======================================================================
  // Declare utility variables
  //======================================================================
  int
  i,        // index for consumers
  j,        // index for products
  t,        // index for time periods
  k,        // index for parameters
  imh,      // index for MC iterations
  vmh,      // index of current iter in Pseudo FP history (1<=vmh<=psd.n_neighbor)
  mcID,     // index for Monte Carlo experiments
  seed,
  ptAccpt,
  neighborPast, // index of nearest neighbor of accepted MCMC draw
  n_neighbor;
  
  double
  pllh,             // past-iteration log quasi-posterior
  *jump;            // vector of jump sizes for M-H procedure
  
  struct FXP fxp;       // store solution of nested fixed points
  struct PAR par_a;     // store accepted parameters
  struct PAR par_c;     // store candidate parameters
  struct RND rnd;       // store random numbers
  struct VAR var;       // store frequently used variables
  struct DATA data;     // store data
  struct DIAG diag;     // store diagnostics as Lipschitz constant, fixed point iterations
  struct HIST hist;
  struct INCV incv;     // store state transition distribution
  struct POST post;     // store posterior outputs
  struct PRED pred;     // store predicted values
  
  char
  filenameData[1024], filenameRand[1024], filenameParam[1024], filenameSeed[1024],
  paramNameDummy[n_param][512],
  parName[8][512],
  filenameMCMC[1024],
  filenameEstim[1024];
  
  FILE *inp_par, *inpData, *inpRand, *outp, *outp_stat;
  
  //======================================================================
  // Setup input files & parameters
  //======================================================================
  strcpy(filenameData,  directoryInput);
  strcpy(filenameRand,  directoryInput);
  strcpy(filenameParam, directoryInput);
  strcpy(filenameSeed,  directoryInput);
  strcat(filenameData,  "dat_dynblp4.txt");
  strcat(filenameRand,  "rand.txt");
  //strcat(filenameParam, "paramInputLow.txt");
  strcat(filenameParam, "paramInput.txt");
  strcat(filenameSeed,  "seed.txt");
  strcpy(filenameEstim, "estimResult.txt");
  strcpy(filenameMCMC, "outMCMC.txt");
  printf("File name for posterior mean output: %s\n", filenameEstim);
  printf("File name for MCMC simulation output: %s\n", filenameMCMC);
  
  var.tol_blp     = tol_blp_pfp;
  var.tol_bellman = tol_blm_pfp;
  
  //======================================================================
  // Memory allocation
  //======================================================================
  gsl_rng *rng = gsl_rng_alloc (gsl_rng_mt19937);
  memoryAllocate(&jump,fxp,par_c,par_a,rnd,var,data,hist,incv,post,pred);
  
  //============================================================================
  // Read in the starting parameter values
  //============================================================================
  inp_par = fopen(filenameParam, "r");
  for (i=1; i<=n_param; i++){
    fscanf(inp_par, "%lf %s", &hist.paramIn[i], paramNameDummy[i-1]);
    printf("%12s = %6.3f\n", paramNameDummy[i-1], hist.paramIn[i]);
  }
  fclose(inp_par);
  
  strcpy(parName[0], "sdBeta1");
  strcpy(parName[1], "sdBeta2");
  strcpy(parName[2], "sdAlpha");
  strcpy(parName[3], "beta0");
  strcpy(parName[4], "beta1");
  strcpy(parName[5], "beta2");
  strcpy(parName[6], "beta3");
  strcpy(parName[7], "alpha");
  
  for (k=1; k<=n_theta2; k++){
    jump[k] = hist.paramIn[n_theta2+k];
    printf("Random-walk MH step jump size for %10s = %6.3f\n",parName[k-1],jump[k]);
  }
  
  //============================================================================
  // Open file streams
  //============================================================================
  // open output file streams
  outp = fopen(filenameMCMC, "w");
  fprintf(outp, "%4s %7s %5s ","mcID","Iter","Seed");
  for (k=1; k<=n_theta1+n_theta2; k++)
    fprintf(outp, "%20s ", parName[k-1]);
  fprintf(outp,"%16s %16s %4s %4s %4s %4s %8s %14s\n",
          "AcceptQpost","RejectQpost","Rej","vmh","nNew","nHist","AccProb","CPUtime");
  fclose(outp);
  
  outp_stat = fopen(filenameEstim, "w");
  fprintf(outp_stat,"%4s %4s %4s %6s ","MC","Seed","Exit","Iter");
  for (k=1; k<=n_theta2; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  for (k=n_theta2+1; k<=n_theta2+n_theta1; k++){
    fprintf(outp_stat,"%12s %12s ",parName[k-1],"SE");
  }
  fprintf(outp_stat,"%12s %14s\n","BellmanIter","CPUtime");
  fclose(outp_stat);
  
  // Read consumer simulation draw
  inpRand = fopen(filenameRand, "r");
  
  // Open input file streams
  inpData = fopen(filenameData, "r");
  
  // Index for obs
  for (t=1; t<=n_period; t++)
    for (j=1; j<=n_product; j++)
      idx[t][j] = n_product*(t-1) + j;
  
  readSeed(hist.paramSeed,filenameSeed);
  
  //======================================================================
  // Monte Carlo loop
  //======================================================================
  // Skip data prior to mcID
  for (mcID=1; mcID<=mc_init-1; mcID++){
    readData(mcID, inpData, var, data);
    readRand(mcID, inpRand, rnd);
  }
  
  for (mcID=mc_init; mcID<=mc_last; mcID++){
    readData(mcID, inpData, var, data);
    readRand(mcID, inpRand, rnd);
    
    for (seed=1; seed<=n_seed; seed++){
      //----------------------------------------------------------------------
      // Preparation & Initialization
      //----------------------------------------------------------------------
      initialize(mcID,seed,imh,vmh,n_neighbor,neighborPast,ptAccpt,
                 fxp,par_a,par_c,var,data,diag,hist,rng);
      diag.tic = gettimeofday_sec();
      nfpBLP0(fxp.expDeltaOut, fxp.valuefunOut, var.expMu0,
              fxp.expDeltaIn, fxp.valuefunIn,
              par_a, rnd, var, data, diag, incv, pred);
      logQuasiposterior(pllh, fxp.expDeltaOut, par_a, var);
      
      // Write out the initial parameter values
      outp = fopen(filenameMCMC, "a");
      printMCMC(mcID,vmh,imh,seed,hist.accept[vmh],neighborPast,n_neighbor,pllh,pllh,par_a,diag,outp);
      
      //======================================================================
      // MCMC Loop
      //======================================================================
      while (diag.exit>0 & imh<n_mh){
        imh++;
        mcmc(mcID,seed,imh,vmh,neighborPast,n_neighbor,ptAccpt,pllh,jump,
             fxp,par_a,par_c,rnd,var,data,diag,hist,incv,post,pred,rng,outp);
        if (imh==diag.Jmax && diag.exit>0){ // Convergence test
          mcmcDiag(filenameEstim, imh, mcID, seed, diag, post, R);
        }
      }
      if (diag.exit!=0){
        printf("Faied to converge with exit code: %d\n",diag.exit);
        printf("Restarting MCMC...\n");
        reinitialize(imh,vmh,n_neighbor,neighborPast,ptAccpt,fxp,par_a,diag,hist);
        while (diag.exit>0 & imh<n_mh){
          imh++;
          mcmc(mcID,seed,imh,vmh,neighborPast,n_neighbor,ptAccpt,pllh,jump,
               fxp,par_a,par_c,rnd,var,data,diag,hist,incv,post,pred,rng,outp);
          if (imh==diag.Jmax){ // Convergence test
            mcmcDiag(filenameEstim, imh, mcID, seed, diag, post, R);
          }
        }
      }
      fclose(outp);
      
      // Store MC results in output file
      printf("\nMCID %2d, Seed %2d complete with exit:%4d. Time elapsed: %.2f\n",
             mcID,seed,diag.exit,diag.toc-diag.tic);
      printf("--------------------------------------------------------------------\n");
      printOutput(filenameEstim, imh, mcID, seed, diag, post);
      R.parseEval(reset);
    }
  } // end of mcID loop
  
  fclose(inpData);
  fclose(inpRand);
  
  //======================================================================
  // Free up Memory
  //======================================================================
  gsl_rng_free (rng);
  
  memoryDeallocate(&jump,fxp, par_c, par_a, rnd, var, data, hist,incv, post, pred);
  
  printf("Monte Carlo experiment complete.\n");
  return 0;
}


