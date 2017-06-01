#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

double gettimeofday_sec(void);
void error(const char error_text[]);

//========================================================================
// Linear Algebra
//========================================================================
void matvecmul(double **A, double *B, double *AB, long m, long n);
void matmul(double **A, double **B, double **AB, long m, long n, long l);
void transpose(double **A, double **A_t, int m, int n);
void copyVector(double *source, double *target, int n);
void copyMatrix(double **source, double **target, int m, int n);
void copy3darray(double ***source, double ***target, int m, int n, int l);
void logVector(double *source, double *target, int size);

/*========================================================================
  Memory Allocation for Arrays
========================================================================*/
double *dvector(long nl, long nh);
int *ivector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double ***darray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h);
int ***iarray3d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h);
double ****darray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
                		long n4l, long n4h);
int ****iarray4d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		long n4l, long n4h);
double *****darray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
            		     long n4l, long n4h, long n5l, long n5h);
int *****iarray5d(long n1l, long n1h, long n2l, long n2h, long n3l, long n3h,
		  long n4l, long n4h, long n5l, long n5h);
void free_dvector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_darray3d(double ***m, long n1l, long n1h, long n2l, long n2h,
             		   long n3l, long n3h);
void free_iarray3d(int ***m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h);
void free_darray4d(double ****m, long n1l, long n1h, long n2l, long n2h,
            		   long n3l, long n3h, long n4l, long n4h);
void free_iarray4d(int ****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h);
void free_darray5d(double *****m, long n1l, long n1h, long n2l, long n2h,
            		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h);
void free_iarray5d(int *****m, long n1l, long n1h, long n2l, long n2h,
		   long n3l, long n3h, long n4l, long n4h, long n5l, long n5h);

