#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
);

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
);

void logQuasiposterior
(
 //=======================================================
 double &llh,
 double *expDeltaIn,
 struct PAR &par,
 struct VAR &var
//=======================================================
);

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
);
