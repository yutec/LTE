#include <stdio.h>
#include "knitro.h"

double computeLogQuasiposterior
(
 //=======================================================
 struct PAR &par,
 struct FXP &fxp,
 struct RND &rnd,
 struct VAR &var,
 struct DAT &dat,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);

void computeExpectedValuefun
(
 //=======================================================
 double *expectValuefun,
 double *valuefunIn,
 double *weight,
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
);

void bellmanOperator
(//=======================================================
 double *valuefunOut,
 double *valuefunIn,
 struct VAR &var,
 struct INCV &incv
//=======================================================
);

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
);

void computeInclusiveValue
(
 //=======================================================
 struct INCV &incv,
 double *expDeltaIn,
 double **expMu,
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

void computeChoicePath
(
 //=======================================================
 struct PRED &pred,
 struct VAR &var
//=======================================================
);

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
 );

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
);

int deltaUpdate
(
 double *expDeltaOut,
 double *expDeltaIn,
 struct DAT  &dat,
 struct PRED &pred
);
