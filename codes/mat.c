/*
  mat.c - Matrix Routines
*/

#ifdef ANSI
#define VOID void
#else
#define VOID
#endif

#include "mat.h"

/*... Function vec_alloc ...*/
/* allocates space for a vector V of n double precision values */
/* returns the pointer to V, access elements V[1..n] */

double *vec_alloc (int n)
{
 double *V ;
 V = (double *) calloc(n, sizeof(double)) ;
 if (V==NULL)
  {
   fprintf (stderr, "could not allocate memory") ;
   exit (1) ;
  }
 return (V-1) ;
}

/*... Function vec_free ...*/
/* free the space used by a vector V of n double precision values */

void vec_free(double *V)
{
 free(V+1);
}

/*... Function ivec_alloc ...*/
/* allocates space for a vector IV of n integer values */
/* returns the pointer to IV, access elements IV[1..n] */

int *ivec_alloc (int n)
{
 int *IV ;
 IV = (int *) calloc(n, sizeof(int)) ;
 if (IV==NULL)
  {
   fprintf (stderr,"could not allocate memory") ;
   exit (1) ;
  }
 return (IV-1) ;
}

/*... Function ivec_free ...*/
/* free the space used by a vector IV of n integer values */

void ivec_free(int *IV)
{
 free(IV+1);
}

/*... Function mat_alloc ...*/
/* allocates space for a matrix A of mxn double precision values */
/* returns the pointer to A, access elements A[1..m][1..n] */

double **mat_alloc (int m, int n)
{
 int i ;
 double **A ;

 A = (double **) calloc (m, sizeof (double *)) ;
 if (A==NULL)
  {
   fprintf (stderr, "could not allocate memory") ;
   exit (1) ;
  }
 A -= 1 ;
 for (i=1; i<=m; i++)
  {
   A[i] = (double *) calloc (n, sizeof (double)) ;
   if (A[i]==NULL)
    {
     fprintf (stderr, "could not allocate memory") ;
     exit (1) ;
    }
   A[i] -= 1 ;
  }
 return A ;
}

/*... Function mat_free ...*/
/* free space used by a matrix A of mxn double precision values */
/* allocated by function mat_alloc */

void mat_free(double **A, int nrows)
{
 int i;

 for (i=1; i<=nrows; i++)
   free(A[i]+1);
 free(A+1);
}

/*... Function mat_alloc_0 ...*/
/* allocates space for a matrix A of mxn float values */
/* returns the pointer to A, access elements A[0..m][0..n] */

float **mat_alloc_0 (int m, int n)
{
 int i ;
 float **A ;

 A = (float **) calloc (m, sizeof (float *)) ;
 if (A==NULL)
  {
   fprintf (stderr, "could not allocate memory") ;
   exit (1) ;
  }
 for (i=0; i<m; i++)
  {
   A[i] = (float *) calloc (n, sizeof (float)) ;
   if (A[i]==NULL)
    {
     fprintf (stderr, "could not allocate memory") ;
     exit (1) ;
    }
  }
 return A ;
}

/*... Function mat_free_0 ...*/
/* free space used by a matrix A of mxn float values */
/* allocated by function mat_alloc_0 */

void mat_free_0(float **A, int nrows)
{
 int i;

 for (i=0; i<nrows; i++)
   free(A[i]);
 free(A);
}

/*... Function cmat_alloc ...*/
/* allocates space for a matrix A of mxn unsigned char values */
/* returns the pointer to A, access elements A[0..m-1][0..n-1] */

unsigned char **cmat_alloc (int m, int n)
{
 int i ;
 unsigned char **A ;

 A = (unsigned char **) calloc (m, sizeof (unsigned char *)) ;
 if (A==NULL)
  {
   fprintf (stderr, "could not allocate memory") ;
   exit (1) ;
  }
 for (i=0; i<m; i++)
  {
   A[i] = (unsigned char *) calloc (n, sizeof (unsigned char)) ;
   if (A[i]==NULL)
    {
     fprintf (stderr, "could not allocate memory") ;
     exit (1) ;
    }
  }
 return A ;
}

/*... Function cmat_free ...*/
/* free space used by a matrix A of mxn unsigned char values */
/* allocated by function mat_alloc */

void cmat_free(unsigned char **A, int nrows)
{
 int i;

 for (i=0; i<nrows; i++)
   free(A[i]); /* both of these were with +1 */
 free(A);
}

/* ... Function AtA ... */
   void AtA(double **A, int m, int n, double **L)
   {
   int i,j,k;
   double sum;

   for (i=1; i<=n; i++)
    {
     for (j=1; j<=n; j++)
      {
      sum = 0.0;
      for(k=1; k<=m; k++) sum=sum + A[k][i]*A[k][j];
      L[i][j]=sum;
      }
    }
   }

/* ... Function AtB ... */
   void AtB(double **A, double **B, int m, int n, int l, double **L)
   {
   int i,j,k;
   double sum;

   for (i=1; i<=n; i++)
    {
     for (j=1; j<=l; j++)
      {
      sum = 0.0;
      for(k=1; k<=m; k++) sum=sum + A[k][i]*B[k][j];
      L[i][j]=sum;
      }
    }
   }


/* ... Function Atb ... */
   void Atb(double **A, double *b, int m, int n, double *L)
   {
   int i,j,k;
   double sum;

   for (i=1; i<=n; i++)
     {
     sum = 0.0;
     for (k=1; k<=m; k++) sum=sum + A[k][i]*b[k];
     L[i]=sum;
     }
   }

/* ... Function Ab ... */
   void Ab(double **A, double *b, int m, int n, double *L)
   {
   int i,j,k;
   double sum;

   for (i=1; i<=m; i++)
     {
     sum = 0.0;
     for (k=1; k<=n; k++) sum=sum + A[i][k]*b[k];
     L[i]=sum;
     }
   }

/*... Function AtWA ...*/
/* computes the product of the matrices At.W.A, the LHS of normal equations */
/* A[1..m][1..n], W[1..m][1..m] */
/* returns the pointer to the result, L, access elements L[1..n][1..n] */

   void AtWA(double **A, double **W, int m, int n, double **L)
   {
   double *V;
   int i,j,k,l;

   V = vec_alloc (m) ;

   for (i=1; i<=n; i++)
    {
     for (k=1; k<=m; k++)
      {
       V[k] = 0.0;
       for (l=1; l<=m; l++)
	 V[k] += A[l][i]*W[l][k];
      }
     for (j=1; j<=n; j++)
      {
	L[i][j] = 0.0;
	for ( k=1; k<=m; k++)
	 L[i][j] += V[k]*A[k][j];
      }
    }
   vec_free(V);
   }

/*... Function AtWB ...*/
/* computes the product of the matrices At.W.B, the RHS of normal equations */
/* A[1..m][1..n], W[1..m][1..m], B[1..m] */
/* returns the pointer to the result, R, access elements R[1..n] */

   void AtWB(double **A, double **W, double *B, int m, int n, double *R)
   {
   double *V;
   int i,j,k,l;

   V = vec_alloc (m) ;

   for (i=1; i<=n; i++)
    {
     R[i] = 0.0;
     for (k=1; k<=m; k++)
      {
       V[k] = 0.0;
       for (l=1; l<=m; l++)
	 V[k] += A[l][i]*W[l][k];
      }
     for ( k=1; k<=m; k++)
      R[i] += V[k]*B[k];
    }
   vec_free(V);
   }


/*... Function AQAt ...*/
/* computes the product of the matrices A.Q.At, A[1..m][1..n], Q[1..n][1..n] */
/* returns the pointer to the result, Qe, access elements Qe[1..m][1..m] */

   void AQAt(double **A, double **Q, int m, int n, double **Qe)
   {
    int i, j, k, l;
    double *V;

    V = vec_alloc (n);

    for (i=1; i<=m; i++)
     {
      for (k=1; k<=n; k++)
       {
	V[k] = 0.0;
	for (l=1; l<=n; l++)
	  V[k] += A[i][l]*Q[l][k];
       }
      for (j=1; j<=m; j++)
       {
	Qe[i][j] = 0.0;
	for (l=1; l<=n; l++)
	  Qe[i][j] += V[l]*A[j][l];
       }
     }
    vec_free(V);
   }


/*... Function MM ...*/
/* multiplies two matrices A[1..m][1..n] B[1..n][1..l] */
/* returns the pointer to AB, access elements AB[1..m][1..l] */

   void MM(double **A, double **B, int m, int n, int l, double **AB)
  {
   int i,j,k;

   for (i=1; i<=m; i++)
     for (j=1; j<=l; j++)
	 AB[i][j] = 0.0;
   for (i=1; i<=m; i++)
     for (j=1; j<=n; j++)
       for (k=1; k<=l; k++)
	 AB[i][k] += A[i][j] * B[j][k] ;

   }


/*... Function ROTM ...*/
/* computes the omega, phi, kappa rotation matrix, M */
/* returns the pointer to M, access elements M[1..3][1..3] */

   void ROTM(double **M, double omega, double phi, double kappa)
   {
   double sw, sp, sk;
   double cw, cp, ck;

   sw = sin (omega);
   cw = cos (omega);
   sp = sin (phi);
   cp = cos (phi);
   sk = sin (kappa);
   ck = cos (kappa);

   M[1][1] = cp*ck;
   M[2][1] = -cp*sk;
   M[3][1] = sp;
   M[1][2] = cw*sk + sw*sp*ck;
   M[2][2] = cw*ck - sw*sp*sk;
   M[3][2] = -sw*cp;
   M[1][3] = sw*sk - cw*sp*ck;
   M[2][3] = sw*ck + cw*sp*sk;
   M[3][3] = cw*cp;

   }

/*... Function Print_Matrix ...*/

void Print_Matrix(double **A, int m, int n, char *s)
{
int i,j,k;
int cc;
div_t qr;

printf("%s\n\n",s);
qr=div(n,7);
cc=1;
for(k=1; k<=qr.quot; k++)
  {
  for(i=1; i<=m; i++)
    {
    for(j=cc; j<=cc+6; j++) printf(" %10.4lf",A[i][j]);
    printf("\n");
    }
  printf("\n");
  cc=cc + 7;
  }

if(qr.rem != 0)
  {
  for(i=1; i<=m; i++)
    {
    for(j=cc; j<=n; j++) printf(" %10.4lf",A[i][j]);
    printf("\n");
    }
  printf("\n");
  }
}

void Print_Vector(double *b, int n, char *s)
{
int i;

printf("%s\n\n",s);
for(i=1; i<=n; i++)
  {
  printf(" %10.4lf\n",b[i]);
  }
printf("\n");
}

void zerom(double **mx, int m, int n)
{
int i,j;

for(i=1; i<=m; i++)
  {
  for(j=1; j<=n; j++)
    {
    mx[i][j]=0.0;
    }
  }
}

/*
  fills the given m by n matrix with 0
  */
int mat_zero(double **mat, int m, int n)
{
  int i, j;

  for (i = 1; i <= m; i++){
    for (j = 1; j <= n; j++){
      mat[i][j] = 0.0;
    }
  }

  return(1);
}

/*
  fills the given m vector with 0
  */
int vec_zero(double *vec, int m)
{
  int i;

  for (i = 1; i <= m; i++){
    vec[i] = 0.0;
  }

  return(1);
}

/*
  adds matrices A and B, places result in C
  (C can be the same as A or B) 
  all matrices are dimensioned m by n 
  */
int mat_add(double **A, double **B, double **C, int m, int n)
{
  int i, j;

  for (i = 1; i <= m; i++){
    for (j = 1; j <= n; j++){
      C[i][j] = A[i][j] + B[i][j];
    }
  }

  return(1);
}

/*
  adds vectors A and B, places result in C
  (C can be the same as A or B) 
  all vectors are dimensioned m by 1 
  */
int vec_add(double *A, double *B, double *C, int m)
{
  int i;

  for (i = 1; i <= m; i++){
      C[i] = A[i] + B[i];
  }

  return(1);
}

/*
  calculates determinant of a 3 by 3 matrix 
  return error (0) if determinant is smaller than tolerance
 */
#define DETERM_TOLERANCE 1.0e-10
int determ_3x3(double **matrix, double *det)
{
  double tmp1, tmp2, tmp3; 

  tmp1 = matrix[2][2] * matrix[3][3] -
    matrix[2][3] * matrix[3][2];
  tmp2 = matrix[2][1] * matrix[3][3] -
    matrix[2][3] * matrix[3][1];
  tmp3 = matrix[2][1] * matrix[3][2] -
    matrix[2][2] * matrix[3][1];

  *det = matrix[1][1] * tmp1 - matrix[1][2] * tmp2
    + matrix[1][3] * tmp3;

  if (fabs(*det) < DETERM_TOLERANCE){
    fprintf(stderr, "determ_3X3: zero determinant \n");
    return(0);
  } 
  return(1);
}

int determ_2x2(double **matrix, double *det)
{
  *det = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];

  if (fabs(*det) < DETERM_TOLERANCE){
    fprintf(stderr, "determ_2X2: zero determinant \n");
    return(0);
  } 
  return(1);
}
