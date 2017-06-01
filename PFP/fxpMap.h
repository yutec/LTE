void nfpBLPcore
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);

void nfpBLP
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);

void nfpBLP0
(//=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double ***expMuOut,
 double *expDeltaIn,
 double **valuefunIn,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);

void blpMapMu
(
 //=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double *expDeltaIn,
 double **valuefunIn,
 double ***expMuIn,
 struct PAR &par,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);

void blpMap
(
 //=======================================================
 double *expDeltaOut,
 double **valuefunOut,
 double ***expMu,
 double *expDeltaIn,
 double **valuefunIn,
 struct PAR &par,
 struct RND &rnd,
 struct VAR &var,
 struct DATA &data,
 struct DIAG &diag,
 struct INCV &incv,
 struct PRED &pred
//=======================================================
);


void deltaUpdate
(
 double *expDeltaOut,
 double *expDeltaIn,
 struct DATA &data,
 struct DIAG &diag,
 struct PRED &pred
);
