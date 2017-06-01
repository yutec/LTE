#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <errno.h>
#include "tools.h"

/*========================================================================
  Get time
========================================================================*/
double gettimeofday_sec(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}


/*========================================================================
  error
========================================================================*/
void error(const char error_text[])
{
	printf("%s\n",error_text);
  exit(1);
}


//========================================================================
//	Linear Algebra
//========================================================================

void matvecmul(double **A, double *B, double *AB, long m, long n)
{
  int i,j;
  for (i=0; i<m; i++){
    AB[i] = 0.0;
    for (j=0; j<n; j++){
        AB[i] += A[i][j]*B[j];
    }
  }
}

void matvecmul(double **A, double *B, adouble *AB, long m, long n)
{
  int i,j;
  for (i=0; i<m; i++){
    AB[i] = 0.0;
    for (j=0; j<n; j++){
      AB[i] += A[i][j]*B[j];
    }
  }
}

void matvecmul(double **A, adouble *B, adouble *AB, long m, long n)
{
  int i,j;
  for (i=0; i<m; i++){
    AB[i] = 0.0;
    for (j=0; j<n; j++){
      AB[i] += A[i][j]*B[j];
    }
  }
}
void matvecmul(adouble **A, adouble *B, adouble *AB, long m, long n)
{
  int i,j;
  for (i=0; i<m; i++){
    AB[i] = 0.0;
    for (j=0; j<n; j++){
      AB[i] += A[i][j]*B[j];
    }
  }
}

void matmul(double **A, double **B, double **AB, long m, long n, long l)
{
  int i,j,k;
  for (i=0; i<m; i++){
    for (j=0; j<l; j++){
      AB[i][j] = 0.0;
      for (k=0; k<n; k++){
        AB[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}
void matmul(adouble **A, adouble **B, adouble **AB, long m, long n, long l)
{
  int i,j,k;
  for (i=0; i<m; i++){
    for (j=0; j<l; j++){
      AB[i][j] = 0.0;
      for (k=0; k<n; k++){
        AB[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}

void transpose(double **A, double **A_t, int m, int n)
{
  int i,j;
  double temp_copy;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      temp_copy = A[i][j];
      A_t[j][i] = temp_copy;
    }
  }
}

void transpose(adouble **A, adouble **A_t, int m, int n)
{
  int i,j;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      A_t[j][i] = A[i][j];
    }
  }
}

void copyVector(double *source, double *target, int n)
{
  int i;
  for (i=0; i<n; i++){
    target[i] = source[i];
  }
}

void copyVector(const double *source, double *target, int n)
{
  int i;
  for (i=0; i<n; i++){
    target[i] = source[i];
  }
}


void copyMatrix(double **source, double **target, int m, int n)
{
  int i, j;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      target[i][j] = source[i][j];
    }
  }
}

void copy3darray(double ***source, double ***target, int m, int n, int l)
{
  int i, j, k;
  for (i=0; i<m; i++){
    for (j=0; j<n; j++){
      for (k=0; k<l; k++){
        target[i][j][k] = source[i][j][k];
      }
    }
  }
}

void logVector(double *source, double *target, int size)
{
  int i;
  for (i=0; i<size; i++)
    target[i] = log(source[i]);
}

/*========================================================================
  Memory Allocation for Arrays
========================================================================*/

// vector ----------------------------------------------------------------

double *dvector(long n)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
  v = new double[n];
  return v;
}

int *ivector(long n)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
  v = new int[n];
  return v;
}

unsigned int *uivector(long n)
/* allocate an int vector with subscript range v[nl..nh] */
{
  unsigned int *v;
  v = new unsigned int[n];
  return v;
}

// matrix ----------------------------------------------------------------

double **dmatrix(long nrow, long ncol)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i;
	double **m = new double*[nrow];
  m[0] = new double[nrow*ncol];
  for (i=1; i<nrow; i++)
    m[i] = &m[0][i*ncol];
	return m;
}

int **imatrix(long nrow, long ncol)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i;
  int **m = new int*[nrow];
  m[0] = new int[nrow*ncol];
  for (i=1; i<nrow; i++)
    m[i] = &m[0][i*ncol];
  return m;
}

// array3d ---------------------------------------------------------------

double ***darray3d(long n1, long n2, long n3)
{
  long i, j;
  double ***m = new double**[n1];
  m[0] = new double*[n1*n2];
  m[0][0] = new double[n1*n2*n3];
  
  for (i=0; i<n1; i++){
    if (i>0){
      m[i] = m[i-1] + n2;
      m[i][0] = m[i-1][0] + (n2*n3);
    }
    for (j=1; j<n2; j++){
      m[i][j] = m[i][j-1] + n3;
    }
  }
  return m;
}

int ***iarray3d(long n1, long n2, long n3)
{
  long i, j;
  int ***m = new int**[n1];
  m[0] = new int*[n1*n2];
  m[0][0] = new int[n1*n2*n3];
  
  for (i=0; i<n1; i++){
    if (i>0){
      m[i] = m[i-1] + n2;
      m[i][0] = m[i-1][0] + (n2*n3);
    }
    for (j=1; j<n2; j++){
      m[i][j] = m[i][j-1] + n3;
    }
  }
  return m;
}



// free vector ----------------------------------------------------

void free_dvector(double *v)
/* free a double vector allocated with dvector() */
{
  delete [] v;
}

void free_ivector(int *v)
/* free an int vector allocated with ivector() */
{
  delete [] v;
}

void free_uivector(unsigned int *v)
/* free an int vector allocated with ivector() */
{
  delete [] v;
}

// free matrix ----------------------------------------------------

void free_dmatrix(double **m)
/* free a double matrix allocated by dmatrix() */
{
  delete [] m[0];
  delete [] m;
}

void free_imatrix(int **m)
/* free an int matrix allocated by imatrix() */
{
  delete [] m[0];
  delete [] m;
}

// free array3d ----------------------------------------------------------

void free_darray3d(double ***m)
/* free a double 3-dimensional array allocated by darray3d() */
{
  delete [] m[0][0];
  delete [] m[0];
  delete m;
}

void free_iarray3d(int ***m)
/* free a int 3-dimensional array allocated by darray3d() */
{
  delete [] m[0][0];
  delete [] m[0];
  delete m;
}

/*========================================================================
 Memory Allocation for adouble Arrays
 ========================================================================*/

adouble *advector(long n)
/* allocate a double vector with subscript range v[nl..nh] */
{
  adouble *v;
  v = new adouble[n];
  if (!v) error("allocation failure in dvector()");
  return v;
}

adouble **admatrix(long nrow, long ncol)
{
  long i;
  adouble **m = new adouble*[nrow];
  m[0] = new adouble[nrow*ncol];
  for (i=1; i<nrow; i++)
    m[i] = &m[0][i*ncol];
  
  /* return pointer to array of pointers to rows */
  return m;
}

adouble ***adarray3d(long n1, long n2, long n3)
{
  long i, j;
  adouble ***m = new adouble**[n1];
  m[0] = new adouble*[n1*n2];
  m[0][0] = new adouble[n1*n2*n3];
  
  for (i=0; i<n1; i++){
    if (i>0){
      m[i] = m[i-1] + n2;
      m[i][0] = m[i-1][0] + (n2*n3);
    }
    for (j=1; j<n2; j++){
      m[i][j] = m[i][j-1] + n3;
    }
  }
  return m;
}

void free_advector(adouble *v)
/* free a double vector allocated with dvector() */
{
  delete [] v;
}

void free_admatrix(adouble **m)
/* free a double matrix allocated by dmatrix() */
{
  delete [] m[0];
  delete [] m;
}

void free_adarray3d(adouble ***m)
/* free a double 3-dimensional array allocated by darray3d() */
{
  delete [] m[0][0];
  delete [] m[0];
  delete m;
}

