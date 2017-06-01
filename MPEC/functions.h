#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void printOutputMPEC
(
 const int &nnzJ,
 const double * const x,
 double * const obj,
 double * const objGrad,
 struct USERSTRUCT *us
);

void printOutput
(
 char *filenameEstim,
 int    &mcID,
 int    &seed,
 int    &nStatus,
 double &obj,
 double &toc,
 double *x,
 struct USERSTRUCT &us
);

void initializeMu
(
 adouble **expMu,
 double *nu,
 struct PAR_AD &par_ad,
 struct VAR &var
);

void initializeMu
(
 //=======================================================
 double ***expMu,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var
//=======================================================
);

void initialize
(
 int    &mcID,
 int    &seed,
 double *xInitial,
 double *paramSeed,
 struct USERSTRUCT &us
);

void computeMatrices
(double **invPHI,
 double **XZPHIZXZ,
 double **X,
 double **Z,
 double **Z_t
);

void readData
(
 int    &mcID,
 FILE   *inpData,
 FILE   *inpRand,
 char   *filenameValfun,
 struct DAT &dat,
 struct RND &rnd,
 struct VAR &var
);

void readSeed
(
 double ***paramSeed,
 char* filenameSeed
);

void initializeKnitro
(
 int &m,
 int &n,
 int &nnzJ,
 int &nnzH,
 int *cType,
 double *xInitial,
 double *xLoBnds,
 double *xUpBnds,
 double *cLoBnds,
 double *cUpBnds,
 double *paramBound,
 struct USERSTRUCT &us
);

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
);

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
);

void allocateMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct INCV &incv,
 struct PRED &pred
);

void releaseMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct INCV &incv,
 struct PRED &pred
);

void allocateActiveMemory
(
 struct PAR_AD &par_ad,
 struct FXP_AD &fxp_ad,
 struct VAR_AD &var_ad,
 struct INCV_AD &incv_ad,
 struct PRED_AD &pred_ad
);

void releaseActiveMemory
(
 struct PAR_AD &par_ad,
 struct FXP_AD &fxp_ad,
 struct VAR_AD &var_ad,
 struct INCV_AD &incv_ad,
 struct PRED_AD &pred_ad
);

#endif