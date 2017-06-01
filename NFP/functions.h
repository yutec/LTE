void printOutput
(
 char   *filenameEstim,
 int    &mcID,
 int    &seed,
 int    &nStatus,
 double &obj,
 double &toc,
 struct USERSTRUCT &us
);

void initializeMu
(
 double ***expMu,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var
);

void initialize
(
 int mcID,
 int seed,
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
 int mcID,
 FILE* inpData,
 FILE* inpRand,
 struct DAT &dat,
 struct RND &rnd,
 struct VAR &var
);

void readSeed
(
 double ***paramSeed,
 char* filenameSeed
);

void allocateMemory
(
 struct DAT &dat,
 struct FXP &fxp,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DIAG &diag,
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
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
);




