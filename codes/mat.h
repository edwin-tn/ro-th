/*
   mat.h - header file for "mat.c" and "gauss.c"
*/

#ifdef ANSI
#define VOID void
#else
#define VOID
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double *vec_alloc (int n);
void vec_free (double *V);
int *ivec_alloc (int n);
void ivec_free (int *V);
double **mat_alloc (int m, int n);
void mat_free(double **A, int nrows);
float **mat_alloc_0 (int m, int n);
void mat_free_0(float **A, int nrows);
unsigned char **cmat_alloc(int m, int n);
void cmat_free(unsigned char **A, int nrows);
void AtA(double **A, int m, int n, double **L);
void AtB(double **A, double **B, int m, int n, int l, double **L);
void Atb(double **A, double *b, int m, int n, double *L);
void Ab(double **A, double *b, int m, int n, double *L);
void AtWA(double **A,double **W,int m,int n, double **L);
void AtWB(double **A,double **W,double *B,int m,int n, double *R);
void AQAt(double **A,double **Q,int m, int n, double **Qe);
void MM(double **A,double **B, int m,int n,int l,double **AB);
void ROTM(double **M, double omega, double phi, double kappa);
void Gauss (double **A, double *B, int n, double *X);
void Triangular_Decomposition (double **A, int n, double **UL, int *IPS);
void Solve (double **UL, double *B, int n, int *IPS, double *X);
void Gauss_Inverse (double **A, int n);
void HouseholderReduction(double **a, int m, int n);
void Qv(double **a, double *b, int m, int n);
void TriSolve(double **a, double *b, double *x, int m, int n);
void h1(int p, int l, int m, double *v, double *h);
void h2(int p, int l, int m, double *v, double h, double *c);
void Print_Matrix(double **A, int m, int n, char *s);
void Print_Vector(double *b, int n, char *s);
void zerom(double **mx, int m, int n);
int mat_zero(double **mat, int m, int n);
int vec_zero(double *vec, int m);
int mat_add(double **A, double **B, double **C, int m, int n);
int vec_add(double *A, double *B, double *C, int m);
int determ_3x3(double **matrix, double *det);
int determ_2x2(double **matrix, double *det);
