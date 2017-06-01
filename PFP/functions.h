#include <RInside.h>
#include <gsl/gsl_rng.h>

void inclusiveValue
(
 struct INCV &incv,
 double *deltaIn,
 struct PAR &par,
 struct VAR &var,
 double *nu
//=======================================================
);

void inclusiveValue
(
 //=======================================================
 struct INCV &incv,
 double *expDeltaIn,
 double **expMu,
 struct VAR &var
//=======================================================
);

void choicePath
(
 //=======================================================
 struct PRED &predOut,
 struct VAR  &var,
 struct DIAG &diag
//=======================================================
);

void initializeMu
(
 //=======================================================
 double ***expMuOut,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var
//=======================================================
);

void initializeMui
(
 //=======================================================
 double **expMu,
 double *nu,
 struct PAR &par,
 struct VAR &var
//=======================================================
);

void transProb
(
 //=======================================================
 struct VAR &var,
 struct INCV &incv
//=======================================================
);

void expectedValuefun
(
 //=======================================================
 double *expectValuefun,
 double *valuefunIn,
 double *weight,
 struct INCV &incv
//=======================================================
);

void bellmanOperator
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct INCV &incv
//=======================================================
);

void nfpBellman
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct DIAG &diag,
 struct INCV &incv
//=======================================================
);

int findNeighbor
(
 //=======================================================
 int    &n_neighbor,
 double **paramHistIn,
 struct PAR &par
//=======================================================
);

void moveHistory
(
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &n_adjust
);

void moveHistory2
(
 int &vmh,
 int &n_neighbor
);

void computeMatrices
(double **invPHI,
 double **XZPHIZXZ,
 double **X,
 double **Z,
 double **Z_t
);

void printMCMC
(
 int    &mcID,
 int    &vmh,
 int    &imh,
 int    &seed,
 int    &accept,
 int    &neighborNew,
 int    &n_neighbor,
 double &llh,
 double &pllh,
 struct PAR &par_a,
 struct DIAG &diag,
 FILE *fileout
);

void mcmcDiag
(
 char   *filenameEstim,
 int    &imh,
 int    &mcID,
 int    &seed,
 struct DIAG &diag,
 struct POST &post,
 RInside &R
);

void mcmcRunLength
(
 char   *filenameEstim,
 int    &imh,
 int    &mcID,
 int    &seed,
 struct DIAG &diag,
 struct POST &post,
 RInside &R
);

void printOutput
(
 char   *filenameEstim,
 int    &imh,
 int    &mcID,
 int    &seed,
 struct DIAG &diag,
 struct POST &post
);

Rcpp::NumericMatrix createMatrix
(
 double **param,
 const int m,
 const int n
);

Rcpp::NumericMatrix createMatrix
(
 //=======================================================
 double **param,
 const int m1,
 const int m2,
 const int n
//=======================================================
);

void posteriorMean
(
 //=======================================================
 int    &imh,
 struct DIAG &diag,
 struct POST &post
//=======================================================
);

void readData
(
 int   &mcID,
 FILE* inpData,
 struct VAR &var,
 struct DATA &data
);

void readRand
(
 int   &mcID,
 FILE* inpRand,
 struct RND rnd
);

void readSeed
(
 double ***paramSeed,
 char* filenameSeed
);

void initialize
(
 int &mcID,
 int &seed,
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &neighborPast,
 int &ptAccpt,
 struct FXP &fxp,
 struct PAR &par_a,
 struct PAR &par_c,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct HIST &hist,
 gsl_rng *rng
);

void reinitialize
(
 int &imh,
 int &vmh,
 int &n_neighbor,
 int &neighborPast,
 int &ptAccpt,
 struct FXP &fxp,
 struct PAR &par_a,
 struct DIAG &diag,
 struct HIST &hist
);

void resetHistory
(
 int indx,
 struct HIST &hist
);

void resetHistParam
(
 int idFirst,
 int idLast,
 struct HIST &hist
);

void memoryAllocate
(
 double **jump,
 struct FXP &fxp,
 struct PAR &par_c,
 struct PAR &par_a,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct HIST &hist,
 struct INCV &incv,
 struct POST &post,
 struct PRED &pred
 );

void memoryDeallocate
(
 double **jump,
 struct FXP &fxp,
 struct PAR &par_c,
 struct PAR &par_a,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct HIST &hist,
 struct INCV &incv,
 struct POST &post,
 struct PRED &pred
 );

