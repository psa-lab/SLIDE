#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

/** Matrix Multiplication: c=a.b **/
void mat_mult(double a[DIM][DIM], double b[DIM], double c[DIM])
{
	int i, k;
	for (i = 0; i < DIM; i++)
	{
		c[i] = 0;
		for (k = 0; k < DIM; k++)
			c[i]  +=  a[i][k] * b[k];
	}
}

/** Print out index column of the matrices **/
void print_mat(double matr1[DIM][MAXSIZE], double matr2[DIM][MAXSIZE], double matr3[DIM][MAXSIZE], int index)
{
	int i=0;

	for(i=0;i<DIM;i++)
		printf("%f ", matr1[i][index]);
	for(i=0;i<DIM;i++)
		printf("%f ", matr2[i][index]);
	for(i=0;i<DIM;i++)
		printf("%f ", matr3[i][index]);
	printf("\n");
}


/** Debug stub: print only one matrix **/
void print_one_mat(double matr[DIM][DIM])
{
	int i=0,j=0;

	for(i=0;i<DIM;i++)
	{
		for(j=0;j<DIM;j++)
			printf("%f ", matr[i][j]);
		printf("\n");
	}
	printf("\n");
}

/** calculates inverse of a matrix **/
void invert(double a[DIM][DIM], double res[DIM][DIM])
{
	double det;
	int i,j;
	double temp[DIM][DIM];
	double D[DIM][DIM];

	D[0][0] = a[1][1]*a[2][2]-a[1][2]*a[2][1];
	D[1][0] = -(a[0][1]*a[2][2]-a[0][2]*a[2][1]);
	D[2][0] = a[0][1]*a[1][2]-a[0][2]*a[1][1];

	D[0][1] = -(a[1][0]*a[2][2]-a[1][2]*a[2][0]);
	D[1][1] = a[0][0]*a[2][2]-a[0][2]*a[2][0];
	D[2][1] = -(a[0][0]*a[1][2]-a[0][2]*a[1][0]);

	D[0][2] = a[1][0]*a[2][1]-a[1][1]*a[2][0];
	D[1][2] = -(a[0][0]*a[2][1]-a[0][1]*a[2][0]);
	D[2][2] = a[0][0]*a[1][1]-a[0][1]*a[1][0];

	det = (a[0][0]*D[0][0]) + (a[1][0]*D[1][0]) + (a[2][0]*D[2][0]);

	if(det==0)
	{
		fprintf(stderr, "ERROR: Determinant zero while inverting matrix.  Exiting...\n");
		exit(1);
	}

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			temp[i][j]=D[i][j]/det;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			res[i][j]=temp[j][i];
}

/** transforms matrix A0[], B0[], C0[] in a series of steps as described 
    below. The resulting matrix is stored in R **/

void mat_trans(double R[DIM][MAXSIZE], double A0[DIM][MAXSIZE], double B0[DIM][MAXSIZE], double C0[DIM][MAXSIZE], double alp, double bet, double bl, int n)
{
	double theta1, theta2, theta3;
	double T1[DIM], T2[DIM][DIM], T3[DIM][DIM], T4[DIM][DIM];
	double T2i[DIM][DIM], T3i[DIM][DIM], T4i[DIM][DIM];
	double P4[DIM], P0[DIM];
	double Temp1[DIM], Temp2[DIM];
	double A1[DIM], B1[DIM], C1[DIM];
	double A2[DIM], B2[DIM], C2[DIM];
	double A3[DIM], B3[DIM], C3[DIM];
	double A4[DIM], B4[DIM], C4[DIM];
	double A6[DIM], B6[DIM], C6[DIM];

	int i,j;

	/** Step 1. Translation of atoms to get them to Origin **/
	for(i=0; i<DIM; i++)
		T1[i] = B0[i][n-1];

	for(i=0; i<DIM; i++)
	{
		A1[i]=A0[i][n-1]-T1[i];
		C1[i]=C0[i][n-1]-T1[i];
		B1[i]=B0[i][n-1]-T1[i];
	}
	/**print_mat(A1, B1, C1, n-1);**/

	/** Step 2. Rotation of the coordinate system around z-axis, 
			in order to bring the projection of A on the x-axis **/

	if(A1[0]>0)
		theta1 = atan(A1[1]/A1[0]);
	else if(A1[0]==0)
		theta1 = PI/2;
	else
		theta1 = PI + atan(A1[1]/A1[0]);
	if(isnan(theta1))
	{ 
		fprintf(stderr, "Panic: NaN in Step 2.\n"); 
		exit(1); 
	}

	T2[0][0]=cos(theta1);
	T2[0][1]=sin(theta1);
	T2[1][0]=-sin(theta1);
	T2[1][1]=cos(theta1);
	T2[0][2]=T2[1][2]=T2[2][0]=T2[2][1]=0;
	T2[2][2]=1;

	/**print_one_mat(T2);**/ 

	mat_mult(T2, A1, A2);
	mat_mult(T2, B1, B2);
	mat_mult(T2, C1, C2);

	/**print_mat(A2, B2, C2, n-1);**/

	/** Step 3. Rotation of the coordinate system around y-axis, 
			in order to bring the projection of A on the x-axis **/

	if(A2[0]==0)
		theta2 = PI/2;
	else
		theta2 = atan(A2[2]/A2[0]);

	if(isnan(theta2))
	{ 
		fprintf(stderr, "Panic: NaN in Step 3.\n"); 
		exit(1); 
	}

	T3[0][0]=cos(theta2);
	T3[0][2]=sin(theta2);
	T3[2][0]=-sin(theta2);
	T3[2][2]=cos(theta2);
	T3[0][1]=T3[1][0]=T3[1][2]=T3[2][1]=0;
	T3[1][1]=1;

	/*print_one_mat(T3);*/ 

	mat_mult(T3, A2, A3);
	mat_mult(T3, B2, B3);
	mat_mult(T3, C2, C3);

	/**print_mat(A3, B3, C3, n-1);**/

	/** Step 4. Rotation of the coordinate system, 
			in order to bring C in the xOy-axis **/

	if(C3[1] > 0)
		theta3 = atan(C3[2]/C3[1]);
	else if(C3[1] > 0)
		theta3 = PI/2;
	else
		theta3 = PI + atan(C3[2]/C3[1]);

	if(isnan(theta3))
	{ 
		fprintf(stderr, "Panic: NaN in Step 4.\n"); 
		exit(1); 
	}

	T4[1][1]=cos(theta3);
	T4[1][2]=sin(theta3);
	T4[2][1]=-sin(theta3);
	T4[2][2]=cos(theta3);
	T4[0][1]=T4[0][2]=T4[1][0]=T4[2][0]=0;
	T4[0][0]=1;

	/*print_one_mat(T4); */

	mat_mult(T4, A3, A4);
	mat_mult(T4, B3, B4);
	mat_mult(T4, C3, C4);

	/**print_mat(A4, B4, C4, n-1);**/

	/** Step 5. Coordinates of point P in the transformed coordinate
			system, at distance d from B, at angle alpha between BP' and BA 
			(P' is the projection of point P in the xOy plane), and at angle 
			beta between BP and BP' **/
      
        /* This is incorrect -- what we want is all of the points to be the 
         * same distance from the plane and 20 degrees between them.  
         * Preferably they should be spaced a bit over 1.0 (A) apart, but we
         * are ignoring the spacing since we dont' want to do lots of testing
         * at this time.
         *
         * The "old" method is to effectively space the points 20 degrees apart
         * in the plane and then proceed to rotate that plane 20 degrees.  The
         * effect is the points closer to the axis of rotation are rotated less
         * than those closer to the X axis.
         *
         * What we really want is to reverse the order of operations to 
         * achieve the correct global effect.  That is, we first rotate about
         * the Y-axis by beta ( to get the correct altitude-- Z component), 
         * then rotate the "new" points about the global Z-axis by the angle 
         * alpha to get the desired spread (again this is ok, since rotation 
         * about the Z-axis will keep the Z component constant).
         *
         * As alluded to above, it is known that the points in and out of plane
         * will be closer than 1.0 (A) and result in some clustering of points.
         * At the present, there is limited resources in the lab and it was 
         * deemed more effective to pursue other ideas.
         */

	P4[0] = bl*cos(alp);
	P4[1] = bl*cos(bet)*sin(alp);
	P4[2] = bl*sin(bet)*sin(alp);

        /* toneroma 2009_04_14 -  this code change did not perform as well as 
         * previous code with current scoring */
        /* vanvoor4 June 18, 2009 -- That is because I messed up the rotation
         * matrices and didn't realize that the rotations on mathworld are to
         * rotate the coordinate system.  We want to rotate the vector, and
         * keep the local coordinate system fixed.  That is, use the transpose
         * of the rotation matrices presented on Mathworld. */
	/* 
        P4[0] = bl*cos(bet)*cos(alp);  
        P4[1] = bl*cos(bet)*sin(alp);      
        P4[2] = -1.0 *bl*sin(bet);
        */

	/**print_mat(P4, B4, C4, n-1); **/

	invert(T2, T2i);
	invert(T3, T3i);
	invert(T4, T4i);

	mat_mult(T4i, P4, Temp1);
	mat_mult(T3i, Temp1, Temp2);
	mat_mult(T2i, Temp2, P0);

	for(i=0;i<DIM;i++)
	{
		P0[i]+=T1[i];
		R[i][n-1] = P0[i];
	}

	/** Step 6. Verification **/

	/**mat_mult(T4i, A4, Temp1);
		mat_mult(T3i, Temp1, Temp2);
		mat_mult(T2i, Temp2, A6);
		for(i=0;i<DIM;i++)
		{
			A6[i][n-1]+=T1[i][n-1];
			printf("%f %f ",A0[i][n-1], A6[i][n-1]);
		}
		mat_mult(T4i, B4, Temp1);
		mat_mult(T3i, Temp1, Temp2);
		mat_mult(T2i, Temp2, B6);
		for(i=0;i<DIM;i++)
		{
			B6[i][n-1]+=T1[i][n-1];
			printf("%f %f ",B0[i][n-1], B6[i][n-1]);
		}
		mat_mult(T4i, C4, Temp1);
		mat_mult(T3i, Temp1, Temp2);
		mat_mult(T2i, Temp2, C6);
		for(i=0;i<DIM;i++)
		{
			C6[i][n-1]+=T1[i][n-1];
			printf("%f %f ",C0[i][n-1], C6[i][n-1]);
		}
		printf("\n");
		**/

	/**printf("Total = %d\n", n);**/
}
