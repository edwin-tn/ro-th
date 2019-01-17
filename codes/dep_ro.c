/*

  dependent relative orientation using the coplanarity condition

  PC coordinates of left image taken as 0, 0, 0
  orientation of left image taken as omega = 0, phi=0, kappa = 0
  focal length of both images assumed to be the same 
  input image coordinates assumed to be relative to principal point and 
  corrected for systematic errors 

  Source: "Introduction to Modern Photogrammetry"", 2001 -- Edward M. Mikhail, James S. Bethel, J. Chris McGlone
  
 */


#include <stdio.h>
#include <math.h>
#include "mat.c"
#include "gauss.c"


#define DEBUG 0

int main(int argc, char **argv)
{
  double *params;    /* parameter values bx, by, omega2, phi2, kappa2 */
  double **a;             /* partials wrt observations */
  double **b;             /* partials wrt parameters */
  double *f;              /* discrepancy vector */
  double **n;             /* normal eqn coefficient matrix */
  double *t;              /* normal eqn constant matrix */
  double **n_temp;
  double *t_temp;
  double *temp;
  double *temp2;
  double *temp3;
  double **temp4; 
  double *delta;          /* corrections to parameters */
  double **transmat;      /* rotation matrix */
  double **left;          /* point coordinates on left image  */
  double **right;         /* point coordinates on right image  */
  double *resids;         /* residuals */
  double *uvw1;           /* image coordinates on left image, 
			     in world coord system */
  double *uvw2;           /* image coordinates on right image, 
			     in world coord system */
  double we;              /* equivalent weight matrix (A At) inverse */
  double **determ;        /* determinant */
  double **determ2;       /* determinant */

  double Bx;              /* X component of base.  This is fixed, to set model scale */
  double focal;           /* focal length, mm, assumed to be the same for both images */
  double sin_w, cos_w; 
  int error, iteration, i, j, k, kk;
  int nbr_pts; 
  
  double cosk, sink;
  FILE *datafile; 

  if (argc < 2){
    fprintf(stderr, "dep_ro <datafile>  \n"); 
    exit(0);
  }

  /* read input datafile 
     nbr_pts
     left x, left y, right x, right y
     ...
     ...
     ...
     approx by, approx bz, approx omega2, approx phi2, approx kappa2 
     Bx, f
     */
     
  params = vec_alloc(5);
  datafile = fopen(argv[1], "r");
  if (!datafile){
    fprintf(stderr, "Unable to open datafile '%s' \n", argv[1]);
    exit(0);
  }
  fscanf(datafile, "%d", &nbr_pts);
  if (nbr_pts < 5){
    fprintf(stderr, "Only %d points given: at least 5 points are required \n", nbr_pts);
    exit(0);
  }

  /* allocate point coordinate arrays and read */
  left = mat_alloc(nbr_pts, 2);
  right = mat_alloc(nbr_pts, 2);
  for (i = 1; i <= nbr_pts; i++){
    fscanf(datafile, "%lf %lf %lf %lf ", 
	   &(left[i][1]), &(left[i][2]), &(right[i][1]), &(right[i][2]));
  }

  /* read approximate parameters */
  fscanf(datafile, "%lf %lf %lf %lf %lf ", 
	 &(params[1]), &(params[2]), &(params[3]), 
	 &(params[4]), &(params[5]));

  /* read Bx and focal length */
  fscanf(datafile, "%lf %lf  ", &Bx, &focal);

  /* set up solution matrices */
  transmat = mat_alloc(3, 3);
  a = mat_alloc(1, 4);  
  b = mat_alloc(1, 5);  
  f = vec_alloc(1);  
  n = mat_alloc(5, 5);  
  n_temp = mat_alloc(5, 5); 
  t = vec_alloc(5);  
  t_temp = vec_alloc(5 );
  delta = vec_alloc(5);  
  temp = vec_alloc(1);   
  temp2 = vec_alloc(3);   
  temp3 = vec_alloc(3);   
  temp4 = mat_alloc(3,3);   
  uvw1 = vec_alloc(3);
  uvw2 = vec_alloc(3);
  determ = mat_alloc(3, 3);
  determ2 = mat_alloc(2, 2);

  iteration = 0;
  while(iteration < 10){
    iteration++;
    vec_zero(t, 5);
    mat_zero(n, 5, 5);

    /* calculate orientation matrix for right image, 
       based on current parameter values 
       (left is always identity matrix) */
    ROTM(transmat, params[3], params[4], params[5]);

    /* loop over points, form A, B, and f for each, 
       and add BtWeB and BtWef to the normal equations
    */
    for (i = 1; i <= nbr_pts; i++){
      
      mat_zero(a, 1, 4);
      mat_zero(b, 1, 5);
      
      /* calculate UVW for left and right image coordinates */
      uvw1[1] = left[i][1];  
      uvw1[2] = left[i][2]; 
      uvw1[3] = -focal;
      temp2[1] = right[i][1];       
      temp2[2] = right[i][2];
      temp2[3] = -focal;
      Atb(transmat, temp2, 3, 3, uvw2);

      determ[1][1] = Bx; 
      determ[1][2] = params[1];
      determ[1][3] = params[2];

      determ[3][1] = uvw2[1];
      determ[3][2] = uvw2[2];
      determ[3][3] = uvw2[3];

      determ[2][1] = 1.0;
      determ[2][2] = 0.0;
      determ[2][3] = 0.0;

      determ_3x3(determ, &(a[1][1]));

      determ[2][1] = 0.0;
      determ[2][2] = 1.0;
      determ[2][3] = 0.0;

      determ_3x3(determ, &(a[1][2]));

      determ[2][1] = uvw1[1];
      determ[2][2] = uvw1[2];
      determ[2][3] = uvw1[3];

      determ[3][1] = transmat[1][1];
      determ[3][2] = transmat[1][2];
      determ[3][3] = transmat[1][3];
      determ_3x3(determ, &(a[1][3]));

      determ[3][1] = transmat[2][1];
      determ[3][2] = transmat[2][2];
      determ[3][3] = transmat[2][3];
      determ_3x3(determ, &(a[1][4]));

      /* now have A, calc we */
      we = 0.0;
      for (j = 1; j <= 4; j++)
	we += a[1][j] * a[1][j];  /* this is actually Qe */
      we = 1.0 / we;              /* invert Qe to get We */

      /* now form B matrix for this point */

      determ2[1][1] = uvw1[1];
      determ2[1][2] = uvw1[3];
      determ2[2][1] = uvw2[1];
      determ2[2][2] = uvw2[3];
      determ_2x2(determ2, &(b[1][1]));
      b[1][1] = -b[1][1];

      determ2[1][2] = uvw1[2];
      determ2[2][2] = uvw2[2];
      determ_2x2(determ2, &(b[1][2]));
      
      /* wrt omega 2 */
      determ[1][1] = Bx; 
      determ[1][2] = params[1];
      determ[1][3] = params[2];

      determ[2][1] = uvw1[1];
      determ[2][2] = uvw1[2];
      determ[2][3] = uvw1[3];

      determ[3][1] = 0.0;
      determ[3][2] = -uvw2[3];
      determ[3][3] = uvw2[2];
      determ_3x3(determ, &(b[1][3]));

      /* wrt phi2 */
      sin_w = sin(params[3]);  /* sin omega */ 
      cos_w = cos(params[3]);  /* cos omega */
      mat_zero(temp4, 3, 3); 
      temp4[1][2] = -sin_w;
      temp4[1][3] =  cos_w;
      temp4[2][1] =  sin_w;
      temp4[3][1] = -cos_w;
      Ab(temp4, uvw2, 3, 3, temp3); 

      determ[3][1] = temp3[1];
      determ[3][2] = temp3[2];
      determ[3][3] = temp3[3];
      determ_3x3(determ, &(b[1][4]));

      /* wrt kappa2 */
      temp2[1] = -right[i][2];       
      temp2[2] = right[i][1];
      temp2[3] = 0.0;
      Atb(transmat, temp2, 3, 3, temp3);

      determ[3][1] = temp3[1];
      determ[3][2] = temp3[2];
      determ[3][3] = temp3[3];
      determ_3x3(determ, &(b[1][5]));
      

      determ[3][1] = uvw2[1];
      determ[3][2] = uvw2[2];
      determ[3][3] = uvw2[3];
      determ_3x3(determ, &(f[1]));
      f[1] = -f[1];

      /* 
	 since we assume that the equations for each point are independent,
	 the normal equations can be formed by summing the BtB and Btf
	 for each point, instead of forming the complete B and f matrices
	 and multiplying them to form the final N and t matrices
       */

      AtA(b, 1, 5, n_temp);
      Atb(b, f, 1, 5, t_temp);      

      /* since We is a 1 by 1 matrix, we can treat it as a scalar 
	 and just multiply each element of N_temp and t_temp by it */
      for (k = 1; k <= 5; k++){
	t_temp[k] *= we; 
	for (kk = 1; kk <=5; kk++){
	  n_temp[k][kk] *= we; 
	}
      }

      mat_add(n, n_temp, n, 5, 5);
      vec_add(t, t_temp, t, 5);
    }
    /* calculate parameter corrections and update parameter approximations */
    Gauss(n, t, 5, delta);
    vec_add(params, delta, params, 5);

    printf("Parameter values and corrections, iteration %d \n", iteration);
    printf("By           :  %12.6f %12.6f \n", params[1], delta[1]);
    printf("Bz           :  %12.3f %12.6f \n", params[2], delta[2]);
    printf("omega2 (rad) :  %12.8f %12.8f \n", params[3], delta[3]);
    printf("phi2 (rad)   :  %12.8f %12.8f \n", params[4], delta[4]);
    printf("kappa2 (rad) :  %12.8f %12.8f \n", params[5], delta[5]);

    /* calculate the point residuals */
    printf("\nPoint image residuals \n");
    resids = vec_alloc(4);
    for (i = 1; i <= nbr_pts; i++){
      mat_zero(a, 1, 4);
      mat_zero(b, 1, 5);
      
      /* calculate UVW for left and right image coordinates */
      uvw1[1] = left[i][1];  
      uvw1[2] = left[i][2]; 
      uvw1[3] = -focal;
      temp2[1] = right[i][1];       
      temp2[2] = right[i][2];
      temp2[3] = -focal;
      Atb(transmat, temp2, 3, 3, uvw2);

      determ[1][1] = Bx; 
      determ[1][2] = params[1];
      determ[1][3] = params[2];

      determ[3][1] = uvw2[1];
      determ[3][2] = uvw2[2];
      determ[3][3] = uvw2[3];

      determ[2][1] = 1.0;
      determ[2][2] = 0.0;
      determ[2][3] = 0.0;

      determ_3x3(determ, &(a[1][1]));

      determ[2][1] = 0.0;
      determ[2][2] = 1.0;
      determ[2][3] = 0.0;

      determ_3x3(determ, &(a[1][2]));

      determ[2][1] = uvw1[1];
      determ[2][2] = uvw1[2];
      determ[2][3] = uvw1[3];

      determ[3][1] = transmat[1][1];
      determ[3][2] = transmat[1][2];
      determ[3][3] = transmat[1][3];
      determ_3x3(determ, &(a[1][3]));

      determ[3][1] = transmat[2][1];
      determ[3][2] = transmat[2][2];
      determ[3][3] = transmat[2][3];
      determ_3x3(determ, &(a[1][4]));

      /* now have A, calc we */
      we = 0.0;
      for (j = 1; j <= 4; j++)
	we += a[1][j] * a[1][j];  /* this is actually Qe */
      we = 1.0 / we;              /* invert Qe to get We */

      /* now form B matrix for this point */

      determ2[1][1] = uvw1[1];
      determ2[1][2] = uvw1[3];
      determ2[2][1] = uvw2[1];
      determ2[2][2] = uvw2[3];
      determ_2x2(determ2, &(b[1][1]));
      b[1][1] = -b[1][1];

      determ2[1][2] = uvw1[2];
      determ2[2][2] = uvw2[2];
      determ_2x2(determ2, &(b[1][2]));
      
      /* wrt omega 2 */
      determ[1][1] = Bx; 
      determ[1][2] = params[1];
      determ[1][3] = params[2];

      determ[2][1] = uvw1[1];
      determ[2][2] = uvw1[2];
      determ[2][3] = uvw1[3];

      determ[3][1] = 0.0;
      determ[3][2] = -uvw2[3];
      determ[3][3] = uvw2[2];
      determ_3x3(determ, &(b[1][3]));

      /* wrt phi2 */
      sin_w = sin(params[3]);  /* sin omega */ 
      cos_w = cos(params[3]);  /* cos omega */
      mat_zero(temp4, 3, 3); 
      temp4[1][2] = -sin_w;
      temp4[1][3] =  cos_w;
      temp4[2][1] =  sin_w;
      temp4[3][1] = -cos_w;
      Ab(temp4, uvw2, 3, 3, temp3); 

      determ[3][1] = temp3[1];
      determ[3][2] = temp3[2];
      determ[3][3] = temp3[3];
      determ_3x3(determ, &(b[1][4]));

      /* wrt kappa2 */
      temp2[1] = -right[i][2];       
      temp2[2] = right[i][1];
      temp2[3] = 0.0;
      Atb(transmat, temp2, 3, 3, temp3);

      determ[3][1] = temp3[1];
      determ[3][2] = temp3[2];
      determ[3][3] = temp3[3];
      determ_3x3(determ, &(b[1][5]));
      
      determ[3][1] = uvw2[1];
      determ[3][2] = uvw2[2];
      determ[3][3] = uvw2[3];
      determ_3x3(determ, &(f[1]));
      f[1] = -f[1];

      /* residuals are Q At We (f - B delta) 
	 assume that Q = I */
      Ab(b, delta, 1, 5, temp);
      printf("we %20.12f f %12.6f temp %12.6f \n",
	     we, f[1], temp[1]);
      printf("a %12.6f  %12.6f %12.6f  %12.6f \n",
	     a[1][1], a[1][2], a[1][3], a[1][4]);

      for (k = 1; k <= 4; k++){
	resids[k] = a[1][k] * we * (f[1] - temp[1]);
      }
      printf("Point %4d  x1 %8.3f y1 %8.3f x2 %8.3f y2 %8.3f \n",
	     i, resids[1], resids[2], resids[3], resids[4]);
    } /* end point loop */

  } /* end iteration loop */

  exit(1);
}





