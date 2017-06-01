#ifndef MPEC_H
#define MPEC_H

void constraintOnly
(
 double * const c,
 const double * const x,
 struct USERSTRUCT *us
 );

void tapeJacobian
(
 double * const c,
 const double * const x,
 struct USERSTRUCT &us
);

void constraintJacobian
(
 const int &nnzJ,
 double * const c,
 double * jac,
 const double * const x,
 struct USERSTRUCT &us
);

void tapeHessian
(
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
);

void lagrangianHessian
(
 const int &nnzH,
 double * const hess,
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
);

void lagrangianHessianNoF
(
 const int &nnzH,
 double * const hess,
 const double * const x,
 double * const lambda,
 struct USERSTRUCT *us
);

void gmmObj_and_gradient
(
 //=======================================================
 double *obj,
 double *grad,
 const double * const x,
 struct USERSTRUCT *us
//=======================================================
);

int callbackEvalFC
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

int callbackEvalGA
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

int callbackEvalH
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

#endif