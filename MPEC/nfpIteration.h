#ifndef NFPITERATION_H
#define NFPITERATION_H

#include <stdio.h>
#include "knitro.h"

void computeExpectedValuefun
(
 //=======================================================
 double *expectValuefun,
 double *valuefunIn,
 struct INCV &incv
//=======================================================
);

void computeExpectedValuefun
(
 //=======================================================
 adouble *expectValuefun,
 adouble *valuefunIn,
 struct INCV_AD &incv
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

void bellmanOperator
(//=======================================================
 adouble *valuefunOut,
 adouble *valuefunIn,
 struct VAR_AD &var,
 struct INCV_AD &incv
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
 struct INCV_AD &incv,
 adouble *expDeltaIn,
 adouble **expMu,
 struct VAR_AD &var
//=======================================================
);

void computeInclusiveValue
(//=======================================================
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
void transProb_AD
(
 //=======================================================
 struct VAR_AD &var,
 struct INCV_AD &incv
//=======================================================
);

void computeChoicePath
(
 //=======================================================
 struct PRED &pred,
 struct VAR &var
//=======================================================
);

void computeChoicePath_AD
(
 //=======================================================
 struct PRED_AD &pred,
 struct VAR_AD &var
//=======================================================
);

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
);

#endif