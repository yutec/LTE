/* This tool box contains useful numerical algorithms not 
  defined in ANSI C standard libraries. Some of the algorithms 
  are taken from Numerical Recipes in C.

  First Version: May 23, 2008
  Second Version: July 4, 2008
  Updated
*/

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
  for (i=1; i<=m; i++){
    AB[i] = 0.0;
    for (j=1; j<=n; j++){
        AB[i] += A[i][j]*B[j];
    }
  }
}


void matmul(double **A, double **B, double **AB, long m, long n, long l)
{
  int i,j,k;
  for (i=1; i<=m; i++){
    for (j=1; j<=l; j++){
      AB[i][j] = 0.0;
      for (k=1; k<=n; k++){
        AB[i][j] += A[i][k]*B[k][j];
      }
    }
  }
}

void transpose(double **A, double **A_t, int m, int n)
{
  int i,j;
  double temp_copy;
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      temp_copy = A[i][j];
      A_t[j][i] = temp_copy;
    }
  }
}

void copyVector(double *source, double *target, int n)
{
  int i;
  for (i=1; i<=n; i++){
    target[i] = source[i];
  }
}

void copyMatrix(double **source, double **target, int m, int n)
{
  int i, j;
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      target[i][j] = source[i][j];
    }
  }
}

void copy3darray(double ***source, double ***target, int m, int n, int l)
{
  int i, j, k;
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      for (k=1; k<=l; k++){
        target[i][j][k] = source[i][j][k];
      }
    }
  }
}

void logVector(double *source, double *target, int size)
{
  int i;
  for (i=1; i<=size; i++)
    target[i] = log(source[i]);
}

/*========================================================================
  Memory Allocation for Arrays
========================================================================*/

#define NR_END 1
#define FREE_ARG char*

// vector ----------------------------------------------------------------

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) error("allocation failure in dvector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) error("allocation failure in ivector()");
	return v-nl+NR_END;
}

// matrix ----------------------------------------------------------------

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

// array3d ---------------------------------------------------------------

double ***darray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h] */
{
  long i, j, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1;
  double ***m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double ***) malloc(((n1+NR_END)*sizeof(double**)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double **)malloc(((n1*n2+NR_END)*sizeof(double *)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double *)malloc(((n1*n2*n3+NR_END)*sizeof(double)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }

  /* return pointer to array of pointers to array of pointers to a row of
     memories */
  return m;

}

int ***iarray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h] */
{
  long i, j, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1;
  int ***m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int ***) malloc(((n1+NR_END)*sizeof(int**)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int **)malloc(((n1*n2+NR_END)*sizeof(int *)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int *)malloc(((n1*n2*n3+NR_END)*sizeof(int)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }

  /* return pointer to array of pointers to rows */
  return m;

}

// array4d ---------------------------------------------------------------

double ****darray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1;
  double ****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double ****) malloc(((n1+NR_END)*sizeof(double***)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double ***)malloc(((n1*n2+NR_END)*sizeof(double **)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double **)malloc(((n1*n2*n3+NR_END)*sizeof(double *)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(double *)malloc(((n1*n2*n3*n4+NR_END)*sizeof(double)));
  if (!m[n1l][n2l][n3l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][n2l][n3l] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }

  /* return pointer to array of pointers to rows */
  return m;

}


int ****iarray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h)
/* allocate a double 3-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1;
  int ****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int ****) malloc(((n1+NR_END)*sizeof(int***)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int ***)malloc(((n1*n2+NR_END)*sizeof(int **)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int **)malloc(((n1*n2*n3+NR_END)*sizeof(int *)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(int *)malloc(((n1*n2*n3*n4+NR_END)*sizeof(int)));
  if (!m[n1l][n2l][n3l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][n2l][n3l] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }

  /* return pointer to array of pointers to rows */
  return m;

}


// array5d ---------------------------------------------------------------

double *****darray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		     long n4l, long n4h, long n5l, long n5h)
/* allocate a double 5-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, l, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1,
    n5=n5h-n5l+1;
  double *****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(double *****) malloc(((n1+NR_END)*sizeof(double ****)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(double ****)malloc(((n1*n2+NR_END)*sizeof(double ***)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(double ***)malloc(((n1*n2*n3+NR_END)*sizeof(double **)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(double **)malloc(((n1*n2*n3*n4+NR_END)*sizeof(double *)));
  if (!m[n1l][n2l][n3l]) error("allocation failure 4 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][j][k] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }


  m[n1l][n2l][n3l][n4l]=(double *)malloc(((n1*n2*n3*n4*n5+NR_END)*sizeof(double)));
  if (!m[n1l][n2l][n3l][n4l]) error("allocation failure 5 in matrix()");
  m[n1l][n2l][n3l][n4l] += NR_END;
  m[n1l][n2l][n3l][n4l] -= n5l;


  for(l=n4l+1; l<=n4h; l++) m[n1l][n2l][n3l][l] = m[n1l][n2l][n3l][l-1] + n5;

  for(k=n3l+1; k<=n3h; k++)
    for(l=n4l; l<=n4h; l++) {

      if (l==n4l)
	m[n1l][n2l][k][l] = m[n1l][n2l][k-1][n4h] + n5;
      else
	m[n1l][n2l][k][l] = m[n1l][n2l][k][l-1] + n5;

    }


  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++)
      for(l=n4l; l<=n4h; l++) {

	if ((k==n3l) & (l==n4l))
	  m[n1l][j][k][l] = m[n1l][j-1][n3h][n4h] + n5;
	else {
	  if (l==n4l)
	    m[n1l][j][k][l] = m[n1l][j][k-1][n4h] + n5;
	  else
	    m[n1l][j][k][l] = m[n1l][j][k][l-1] + n5;

	}

      }



  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++)
	for(l=n4l; l<=n4h; l++) {

	  if (((j==n2l) & (k==n3l)) & (l==n4l))
	    m[i][j][k][l] = m[i-1][n2h][n3h][n4h] + n5;
	  else if ((k==n3l) & (l==n4l))
	    m[i][j][k][l] = m[i][j-1][n3h][n4h] + n5;
	  else if (l==n4l)
	    m[i][j][k][l] = m[i][j][k-1][n4h] + n5;
	  else
	    m[i][j][k][l] = m[i][j][k][l-1] + n5;

	}

  /* return pointer to array of pointers to rows */
  return m;

}


int *****iarray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		  long n4l, long n4h, long n5l, long n5h)
/* allocate a double 5-dimensional array with subscript
   range m[n1l..n1h][n2l..n2h][n3l..n3h][n4l..n4h][n5l..n5h] */
{
  long i, j, k, l, n1=n1h-n1l+1, n2=n2h-n2l+1, n3=n3h-n3l+1, n4=n4h-n4l+1,
    n5=n5h-n5l+1;
  int *****m;

  /* allocate pointer-to-pointers to the first dimensional of the array */
  m=(int *****) malloc(((n1+NR_END)*sizeof(int****)));
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= n1l;

  /* allocate pointers and set pointer-to-pointers to them */
  m[n1l]=(int ****)malloc(((n1*n2+NR_END)*sizeof(int ***)));
  if (!m[n1l]) error("allocation failure 2 in matrix()");
  m[n1l] += NR_END;
  m[n1l] -= n2l;

  for(i=n1l+1; i<=n1h; i++) m[i] = m[i-1] + n2;

  /* allocate a row of memories and set pointers to them */
  m[n1l][n2l]=(int ***)malloc(((n1*n2*n3+NR_END)*sizeof(int **)));
  if (!m[n1l][n2l]) error("allocation failure 3 in matrix()");
  m[n1l][n2l] += NR_END;
  m[n1l][n2l] -= n3l;

  for(j=n2l+1; j<=n2h; j++) m[n1l][j] = m[n1l][j-1] + n3;

  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++) {

      if (j==n2l)
	m[i][j] = m[i-1][n2h] + n3;
      else
	m[i][j] = m[i][j-1] + n3;

    }


  m[n1l][n2l][n3l]=(int **)malloc(((n1*n2*n3*n4+NR_END)*sizeof(int *)));
  if (!m[n1l][n2l][n3l]) error("allocation failure 4 in matrix()");
  m[n1l][n2l][n3l] += NR_END;
  m[n1l][n2l][n3l] -= n4l;

  for(k=n3l+1; k<=n3h; k++) m[n1l][n2l][k] = m[n1l][n2l][k-1] + n4;

  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++) {

      if (k==n3l)
	m[n1l][j][k] = m[n1l][j-1][n3h] + n4;
      else
	m[n1l][j][k] = m[n1l][j][k-1] + n4;

    }


  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++) {

	if ((j==n2l) & (k==n3l))
	  m[i][j][k] = m[i-1][n2h][n3h] + n4;
	else {
	  if (k==n3l)
	    m[i][j][k] = m[i][j-1][n3h] + n4;
	  else
	    m[i][j][k] = m[i][j][k-1] + n4;

	}

      }


  m[n1l][n2l][n3l][n4l]=(int *)malloc(((n1*n2*n3*n4*n5+NR_END)*sizeof(int)));
  if (!m[n1l][n2l][n3l][n4l]) error("allocation failure 5 in matrix()");
  m[n1l][n2l][n3l][n4l] += NR_END;
  m[n1l][n2l][n3l][n4l] -= n5l;


  for(l=n4l+1; l<=n4h; l++) m[n1l][n2l][n3l][l] = m[n1l][n2l][n3l][l-1] + n5;

  for(k=n3l+1; k<=n3h; k++)
    for(l=n4l; l<=n4h; l++) {

      if (l==n4l)
	m[n1l][n2l][k][l] = m[n1l][n2l][k-1][n4h] + n5;
      else
	m[n1l][n2l][k][l] = m[n1l][n2l][k][l-1] + n5;

    }


  for(j=n2l+1; j<=n2h; j++)
    for(k=n3l; k<=n3h; k++)
      for(l=n4l; l<=n4h; l++) {

	if ((k==n3l) & (l==n4l))
	  m[n1l][j][k][l] = m[n1l][j-1][n3h][n4h] + n5;
	else {
	  if (l==n4l)
	    m[n1l][j][k][l] = m[n1l][j][k-1][n4h] + n5;
	  else
	    m[n1l][j][k][l] = m[n1l][j][k][l-1] + n5;

	}

      }



  for(i=n1l+1; i<=n1h; i++)
    for(j=n2l; j<=n2h; j++)
      for(k=n3l; k<=n3h; k++)
	for(l=n4l; l<=n4h; l++) {

	  if (((j==n2l) & (k==n3l)) & (l==n4l))
	    m[i][j][k][l] = m[i-1][n2h][n3h][n4h] + n5;
	  else if ((k==n3l) & (l==n4l))
	    m[i][j][k][l] = m[i][j-1][n3h][n4h] + n5;
	  else if (l==n4l)
	    m[i][j][k][l] = m[i][j][k-1][n4h] + n5;
	  else
	    m[i][j][k][l] = m[i][j][k][l-1] + n5;

	}


  /* return pointer to array of pointers to rows */
  return m;

}



// free vector ----------------------------------------------------

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

// free matrix ----------------------------------------------------

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

// free array3d ----------------------------------------------------------

void free_darray3d(double ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h)
/* free a double 3-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray3d(int ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h)
/* free a int 3-dimensional array allocated by darray3d() */
{
        free((m[n1l][n2l]+n3l-NR_END));
	free((m[n1l]+n2l-NR_END));
	free((m+n1l-NR_END));
}


// free array4d ----------------------------------------------------------

void free_darray4d(double ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h)
/* free a double 4-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray4d(int ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h)
/* free a double 4-dimensional array allocated by darray3d() */
{
        free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
	free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
	free((FREE_ARG) (m[n1l]+n2l-NR_END));
	free((FREE_ARG) (m+n1l-NR_END));
}


// free array5d ----------------------------------------------------------

void free_darray5d(double *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h)
/* free a double 5-dimensional array allocated by darray3d() */
{
  free((FREE_ARG) (m[n1l][n2l][n3l][n4l]+n5l-NR_END));
  free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
  free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
  free((FREE_ARG) (m[n1l]+n2l-NR_END));
  free((FREE_ARG) (m+n1l-NR_END));
}

void free_iarray5d(int *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h)
/* free a double 5-dimensional array allocated by darray3d() */
{

        free((FREE_ARG) (m[n1l][n2l][n3l][n4l]+n5l-NR_END));
        free((FREE_ARG) (m[n1l][n2l][n3l]+n4l-NR_END));
	free((FREE_ARG) (m[n1l][n2l]+n3l-NR_END));
	free((FREE_ARG) (m[n1l]+n2l-NR_END));
	free((FREE_ARG) (m+n1l-NR_END));
}

