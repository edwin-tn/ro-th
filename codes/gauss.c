/*
  Bethel - 15-JAN-93

  gauss.c - Gauss equation solver and matrix inversion
*/

#include "mat.h"


void Triangular_Decomposition (double **A, int n, double **UL, int *IPS)
{
 double *Scales;
 int i,j,k;
 int kp,kp1, nm1;
 int ip, idxpiv;
 double row_norm, big, size, pivot, em;

 Scales = vec_alloc(n);

 for (i=1; i<=n; i++)
  {
   IPS[i] = i ;
   row_norm = 0.0 ;
   for (j=1; j<=n; j++)
    {
     UL[i][j] = A[i][j];
     if (fabs(UL[i][j]) > row_norm)
       row_norm = fabs(UL[i][j]);
    }
   if (row_norm == 0.0)
     {
      printf ("error: matrix with zero row in Triangular_Decomposition\n");
      Scales[i] = 0.0 ;
     }
   else Scales[i] = 1.0/row_norm ;
  }

 nm1 = n - 1 ;
 for (k=1; k<=nm1; k++)
  {
   big = 0.0 ;
   for (i=k; i<=n; i++)
    {
     ip = IPS[i] ;
     size = fabs(UL[ip][k])*Scales[ip] ;
     if (size > big)
      {
       big = size ;
       idxpiv = i ;
      }
    }
  if (big == 0.0)
    printf ("error: singular matrix in Triangular_Decomposition\n");
  else
   {
    if (idxpiv != k)
     {
      j = IPS[k] ;
      IPS[k] = IPS[idxpiv] ;
      IPS[idxpiv] = j ;
     }
    kp = IPS[k] ;
    pivot = UL[kp][k] ;
    kp1 = k + 1 ;
    for (i=kp1; i<=n; i++)
     {
      ip = IPS[i] ;
      em = -UL[ip][k]/pivot ;
      UL[ip][k] = -em ;
      for (j=kp1; j<=n; j++)
	UL[ip][j] += em*UL[kp][j] ;
     }
   }
 }

 kp = IPS[n] ;
 if (UL[kp][n] == 0.0)
   printf ("error: singular matrix in Triangular_Decomposition\n") ;

 vec_free (Scales);
}

void Solve (double **UL, double *B, int n, int *IPS, double *X)
{
 int i, j;
 int np1;
 int ip, ip1;
 int iback, im1;
 double sum;

 np1 = n + 1 ;
 ip = IPS[1] ;
 X[1] = B[ip] ;
 for (i=2; i<=n; i++)
  {
   ip = IPS[i] ;
   im1 = i - 1 ;
   sum = 0.0 ;
   for (j=1; j<=im1; j++)
     sum += UL[ip][j]*X[j] ;
   X[i] = B[ip] - sum ;
  }

 ip = IPS[n] ;
 X[n] = X[n]/UL[ip][n] ;

 for (iback=2; iback<=n; iback++)
  {
   i = np1 - iback ;
   ip = IPS[i] ;
   ip1 = i + 1 ;
   sum = 0.0 ;
   for (j=ip1; j<=n; j++)
     sum += UL[ip][j]*X[j] ;
   X[i] = (X[i] - sum)/UL[ip][i] ;
  }
}

void Gauss (double **A, double *B, int n, double *X)
{
 double **UL;
 int *IPS;
 int i,j;

 UL = mat_alloc (n,n);
 IPS = ivec_alloc (n);

 Triangular_Decomposition (A, n, UL, IPS) ;
 Solve (UL, B, n, IPS, X) ;

 ivec_free (IPS);
 mat_free  (UL,n);
}


void Gauss_Inverse (double **A, int n)
{
 double **UL;
 double *Baux;
 double *X;
 int *IPS;
 int i,j;

 UL = mat_alloc (n,n);
 IPS = ivec_alloc (n);
 Baux = vec_alloc (n);
 X = vec_alloc(n);

 Triangular_Decomposition (A, n, UL, IPS) ;

 for (i=1; i<=n; i++)
  {
   for (j=1; j<=n; j++)
     Baux[j] = 0.0;
   Baux[i] = 1.0;
   Solve (UL, Baux, n, IPS, X);
   for (j=1; j<=n; j++)
     A[j][i] = X[j];
  }

  vec_free (X);
  vec_free (Baux);
  ivec_free (IPS);
  mat_free (UL,n);

}
