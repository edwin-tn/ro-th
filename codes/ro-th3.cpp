// ro-th3.cpp: interface for direct relative orientation using Thompson's Parameterization.
// if the baseline's length is known from the coordinates tagged on the EXIF header.
// assume that geotagged coordinates are available and the baseline vector componenets are precalculated in advance.
// The Thompson RO computation is constrained using baseline length
// only 3 rotaninnal parameters are involved
//
// Copyright 2018 Martinus Edwin Tjahjadi
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include "mat.c"
#include "gauss.c"

using namespace std;

void form_B(double** b, double& x1, double& y1, double& x2, double& y2);
void form_f_linear(double* f, double& by, double& bz, double& x1, double& y1, double& x2, double& y2);
void form_f_secondTerm(double* f, double& by, double& bz, double& a1, double& a2, double& a3, 
	double& x1, double& y1, double& x2, double& y2);
void form_f_thirdTerm(double* f, double& by, double& bz, double& a1, double& a2, double& a3,
	double& b1, double& b2, double& x1, double& y1, double& x2, double& y2);
void PrintMatrix(double **A, int m, int n, char *s);
void PrintVector(double *b, int n, char *s);
void form_a1a2a3(double& a1, double& a2, double& a3, double* param, double& x1, double& y1, double&x2, double& y2);
void form_b1b2(double& b1, double& b2, double* param, double& x1, double& y1, double& x2, double& y2);

int main(int argc, char **argv)
{
	double *params = nullptr;		  /* parameter values: omega2, phi2, kappa2 */
	double **left = nullptr;          /* feature point coordinates on left image  */
	double **right = nullptr;         /* feature point coordinates on right image  */
	double **normleft = nullptr;      /* normalized feature point coordinates on left image  */
	double **normright = nullptr;     /* normalized feature point coordinates on right image  */
	int nbr_pts = 0;				  /* Number of correspondent points*/
	const double Bx = 1.0;			  /* X component of base.  This is fixed, to set model scale */
	double by = 0.0;				  /* fixed b'y, obtained from file*/
	double bz = 0.0;			  	  /* fixed b'z, obtained from file*/
	double focal = 0.0;				  /* focal length in mm, assumed to be the same for both images */
	double **b = nullptr;             /* partials wrt parameters */
	double *f = nullptr;              /* discrepancy vector for linear terms, second order and third order terms*/
	double **bTb = nullptr;			  /* matrix b Transpose b */
	double *bTf = nullptr;			  /* vector of b Transpose f */
	double **bTb_temp = nullptr;
	double *bTf_temp = nullptr;
	double a1, a2, a3 = 0.0;		  /* discrepancy vector for second order terms*/
	double delta2, b1, b2 = 0.0;	  /* discrepancy vector for third order terms*/

	FILE *datafile;
	errno_t err;

	if (argc < 2) {
		fprintf(stderr, "ro-th3 <datafile>  \n");
		exit(0);
	}

	/* read input datafile
	nbr_pts
	left x, left y, right x, right y
	...
	...
	...
	known b'y, b'z, approx omega2=0, phi2=0, kappa2=0	//dummy input - not in use
	Bx, focallength
	*/

	params = vec_alloc(3);	//3 params
	err = fopen_s(&datafile, argv[1], "r");
	//datafile = fopen(argv[1], "r");
	if (err) {
		fprintf_s(stderr, "Unable to open datafile '%s' \n", argv[1]);
		exit(0);
	}
	else {
		// Set pointer to beginning of file:
		fseek(datafile, 0L, SEEK_SET);

		//read number of feature points
		fscanf_s(datafile, "%d", &nbr_pts);
		if (nbr_pts < 3) {
			fprintf_s(stderr, "Only %d points given: at least 3 points are required \n", nbr_pts);
			exit(0);
		}

		/* allocate feature point coordinate arrays and read */
		left = mat_alloc(nbr_pts, 2);		//read feature coords on left image
		right = mat_alloc(nbr_pts, 2);		//read feature coords on right image
		for (int i = 1; i <= nbr_pts; i++) {
			fscanf_s(datafile, "%lf %lf %lf %lf ",
				&(left[i][1]), &(left[i][2]), &(right[i][1]), &(right[i][2]));
		}

		/* read approximate parameters -- known b'y and b'z only , zeroes for rotations*/
		fscanf_s(datafile, "%lf %lf %lf %lf %lf ",
			&(by), &(bz), &(params[1]), &(params[2]), &(params[3]));

		/* read Bx and focal length */
		fscanf_s(datafile, "%lf %lf  ", &Bx, &focal);

		/* calculate normalized image coordinates == homogenous coordinates*/
		normleft = mat_alloc(nbr_pts, 2);		//read normalized coords on left image
		normright = mat_alloc(nbr_pts, 2);		//read normalized coords on right image
		for (int i = 1; i <= nbr_pts; i++) {
			normleft[i][1] = left[i][1] / -focal; normleft[i][2] = left[i][2] / -focal;
			normright[i][1] = right[i][1] / -focal; normright[i][2] = right[i][2] / -focal;
		}

	}//else

	/* loop over points, form b and f for each to solve bx=f */
	/* set up solution matrices */
	b = mat_alloc(1, 3);
	f = vec_alloc(1);
	bTb = mat_alloc(3, 3);		mat_zero(bTb, 3, 3);
	bTf = vec_alloc(3);			vec_zero(bTf, 3);
	bTb_temp = mat_alloc(3, 3);
	bTf_temp = vec_alloc(3);
	vec_zero(params, 3);

	/* setup a linear terms of Equation 11 and 18 */
	for (int i = 1; i <= nbr_pts; i++) {
		mat_zero(b, 1, 3);
		vec_zero(f, 1);
		mat_zero(bTb_temp, 3, 3);
		vec_zero(bTf_temp, 3);

		form_B(b, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		form_f_linear(f, by, bz, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		PrintMatrix(b, 1, 3, "Matrix M: ");
		PrintVector(f, 1, "Vector l: ");

		/* since we assume that the equations for each point are independent, the normal equations
		can be formed 	by summing the BtB and Btf 	for each point, instead of forming the complete
		B and f matrices and multiplying them to form the final N and t matrices 	*/

		AtA(b, 1, 3, bTb_temp);		Atb(b, f, 1, 3, bTf_temp);
		PrintMatrix(bTb_temp, 3, 3, "bTb_temp: ");
		PrintVector(bTf_temp, 3, "bTf_temp: ");
		mat_add(bTb, bTb_temp, bTb, 3, 3);
		vec_add(bTf, bTf_temp, bTf, 3);

	}//for linear terms

	PrintMatrix(bTb, 3, 3, "bTb linear terms: ");
	PrintVector(bTf, 3, "bTf linear terms: ");

	/* calculate parameter and update parameter approximations */
	Gauss(bTb, bTf, 3, params);
	PrintVector(params, 3, "Linear solution: ");

	/* setup a second order constant terms*/
	mat_zero(bTb, 3, 3);
	vec_zero(bTf, 3);
	for (int i = 1; i <= nbr_pts; i++) {
		mat_zero(b, 1, 3);
		vec_zero(f, 1);
		mat_zero(bTb_temp, 3, 3);
		vec_zero(bTf_temp, 3);

		form_a1a2a3(a1, a2, a3, params, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);

		form_B(b, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		//form_f_secondTerm(f, normleft[i][2], normright[i][2], a1, a2, a3);
		form_f_secondTerm(f, by, bz, a1, a2, a3, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		Print_Vector(f, 1, "Vector l of second terms: ");

		AtA(b, 1, 3, bTb_temp);
		Atb(b, f, 1, 3, bTf_temp);
		mat_add(bTb, bTb_temp, bTb, 3, 3);
		vec_add(bTf, bTf_temp, bTf, 3);
	}//for second order term

	PrintMatrix(bTb, 3, 3, "bTb of second order term: ");
	PrintVector(bTf, 3, "bTf of second order term: ");
	/* calculate parameter and update parameter approximations */
	Gauss(bTb, bTf, 3, params);
	PrintVector(params, 3, "2nd order solution: ");

	/* setup a third order constant terms*/
	mat_zero(bTb, 3, 3);
	vec_zero(bTf, 3);

	for (int i = 1; i <= nbr_pts; i++) {
		mat_zero(b, 1, 3);
		vec_zero(f, 1);
		mat_zero(bTb_temp, 3, 3);
		vec_zero(bTf_temp, 3);

		form_b1b2(b1, b2, params, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);

		form_B(b, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		//form_f_thirdTerm(f, normleft[i][2], normright[i][2], a1, a2, a3, b1, b2);
		form_f_thirdTerm(f, by, bz, a1, a2, a3, b1, b2, normleft[i][1], normleft[i][2], normright[i][1], normright[i][2]);
		Print_Vector(f, 1, "Vector l of third order terms: ");

		AtA(b, 1, 3, bTb_temp);
		Atb(b, f, 1, 3, bTf_temp);
		mat_add(bTb, bTb_temp, bTb, 3, 3);
		vec_add(bTf, bTf_temp, bTf, 3);
	}//for third order terms

	PrintVector(bTf, 3, "bTf of third order term: ");
	/* calculate parameter and update parameter approximations */
	Gauss(bTb, bTf, 3, params);
	PrintVector(params, 3, "3rd order solution: ");

	/* clean all allocated memories*/
	mat_free(left, nbr_pts); 		mat_free(right, nbr_pts);
	mat_free(normleft, nbr_pts); 	mat_free(normright, nbr_pts);
	mat_free(b, 1);					vec_free(f);
	mat_free(bTb, 3);				vec_free(bTf);
	mat_free(bTb_temp, 3);			vec_free(bTf_temp);
	vec_free(params);


	return 0;
}

void form_B(double** b, double& x1, double& y1, double& x2, double& y2)
{
	b[1][1] = 1 + y1 * y2;			/* omega term: (1+y'1*y'2) */
	b[1][2] = -y1 * x2;				/* phi term: -y'1*x'2 */
	b[1][3] = -x2;					/* kappa term: -x'2 */
}

void form_f_linear(double* f, double& by, double& bz, double& x1, double& y1, double& x2, double& y2)
{
	/* l term: -(y'1 - y'2) + (x1-x2)by - (x1y2 - x2y1)bz*/
	f[1] = -(y1 - y2) + (x1 - x2) * by - (x1 * y2 - x2 * y1) * bz;	
}

void form_f_secondTerm(double* f, double& by, double& bz, double& a1, double& a2, double& a3, 
	double& x1, double& y1, double& x2, double& y2)
{
	/* 2nd term: -(y'1 - y'2) + (x1-x2)by - (x1y2 - x2y1)bz - a1 - a2 - a3 */
	f[1] = -(y1 - y2) + (x1 - x2) * by - (x1 * y2 - x2 * y1) * bz - a1 - a2 - a3;
}

void form_f_thirdTerm(double* f, double& by, double& bz, double& a1, double& a2, double& a3,
	double& b1, double& b2, double& x1, double& y1, double& x2, double& y2)
{
	/* 3rd term: -(y'1 - y'2) + (x1-x2)by - (x1y2 - x2y1)bz - a1 - a2 - a3 - b1 - b2 */
	f[1] = -(y1 - y2) + (x1 - x2) * by - (x1 * y2 - x2 * y1) * bz - a1 - a2 - a3 - b1 - b2;
}

void PrintMatrix(double **A, int m, int n, char *s)
{
	int i, j, k;
	int cc;
	div_t qr;

	printf("%s\n\n", s);
	qr = div(n, 7);
	cc = 1;
	for (k = 1; k <= qr.quot; k++)
	{
		for (i = 1; i <= m; i++)
		{
			for (j = cc; j <= cc + 6; j++) printf(" %10.7lf", A[i][j]);
			printf("\n");
		}
		printf("\n");
		cc = cc + 7;
	}

	if (qr.rem != 0)
	{
		for (i = 1; i <= m; i++)
		{
			for (j = cc; j <= n; j++) printf(" %10.7lf", A[i][j]);
			printf("\n");
		}
		printf("\n");
	}
}

void PrintVector(double *b, int n, char *s)
{
	int i;

	printf("%s\n\n", s);
	for (i = 1; i <= n; i++)
	{
		printf(" %10.7lf\n", b[i]);
	}
	printf("\n");
}


void form_a1a2a3(double& a1, double& a2, double& a3, double* param, double& x1, double& y1, double&x2, double& y2)
{
	double by = param[1];	double bz = param[2];	double omg = param[3];	double phi = param[4];	double kappa = param[5];

	a1 = -0.5  *x2 * phi + 0.5 * x1 *y2 * kappa - x1 * y2 * by - x1 * bz; //Eq.12 a_1=-1⁄2 (x2 ϕ2 )+1⁄2 (x1 y2 κ2 ) -x1 y2 by -x1 bz

	a2 = -0.5 * y2 * phi - 0.5 * (1 - y1 * y2) * kappa + (1 + x1 * x2) * by - y1 * bz;	//Eq.13

	a3 = -y2 * by + (x1 * x2 + y1 * y2) * bz;		//Eq.14

	/* taking care of operating on the coefficients of the second order terms as if they were part of the constant terms */
	/* according to Equation (18): */
	a1 = a1 * omg;	a2 = a2 * phi;	a3 = a3 * kappa;

}

void form_b1b2(double & b1, double & b2, double * param, double & x1, double & y1, double & x2, double & y2)
{
	double by = param[1];	double bz = param[2];	double omg = param[3];	double phi = param[4];	double kappa = param[5];

	double delta2 = 0.25 * (omg * omg + phi * phi + kappa * kappa);

	b1 = 0.5 * (x2 * omg + y2 * phi + kappa) * ((omg - x1 * kappa) * by - (y1 * omg - x1 * phi) * bz);

	b2 = -delta2 * ((y1 - y2) - (x1 - x2) * by - (x1 * y2 - x2 * y1) * bz);
}
