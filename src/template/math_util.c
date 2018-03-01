/************************************************************/
/*          util.c Utility functions for SI.c               */
/************************************************************/

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

/** Routine to sqaure a number **/
double square2(double number) {
  double ansr;
  ansr = number * number;
  return(ansr);
} 

/** Euclidian distance between two points in 3-space **/
double distance2( double pt1[3], double pt2[3] ) {
  double dis=0;
  dis = sqrt( square2( pt2[0] - pt1[0] ) + square2( pt2[1] - pt1[1] ) + 
	      square2( pt2[2] - pt1[2] ) );
  return( dis );
}

/**      Return the smallest of four numbers.       **/
double get_shortest( double a, double b, double c, double d ) {
  double temp1=0, temp2=0, shortest=0;
  temp1 = ( a < b ) ? a : b;
  temp2 = ( temp1 < c ) ? temp1 : c;
  shortest = ( temp2 < d ) ? temp2 : d;
  return( shortest );
}

/** Calculate the angle between two vectors. Returns the acos **/
double angle( double start1[3], double start2[3], double length1, 
	      double length2, double end1[3], double end2[3] ) {
  double vec_product=0, AnGle=0, move1[3], move2[3];
  move1[0] = end1[0] - start1[0];
  move1[1] = end1[1] - start1[1];
  move1[2] = end1[2] - start1[2];
  move2[0] = end2[0] - start2[0];
  move2[1] = end2[1] - start2[1];
  move2[2] = end2[2] - start2[2];
  vec_product = (move1[0] * move2[0]) + (move1[1] * move2[1]) + 
    (move1[2] * move2[2]);
  AnGle = vec_product / ( length1 * length2 );
  return(AnGle);
}
 
/* Calculate vector product of A transpose B */
double A_trans_B( double vec1[3], double vec2[3] ){
  int a;
  double ansr=0.0;
  for(a=0;a<3;a++)
    ansr+=vec1[a]*vec2[a];
  return(ansr);
}

/* Calculate closest dis ( projection ) of a point onto a vector */
double pt_projection( double point[3], double vector_s[3], 
		      double vector_e[3], double *scalar ) {
  int a=0;
  double aTb, aTa, bTb, cd, vec[3], pt[3];
  aTb=aTa=bTb=cd=0;
  while(a<3){
    vec[a] = vector_e[a] - vector_s[a];
    pt[a]  = point[a] - vector_s[a];
    a++;
  }
  aTb = A_trans_B( pt, vec );
  aTa = A_trans_B( vec, vec );
  bTb = A_trans_B( pt, pt );
  *scalar = aTb/aTa;
  cd = sqrt(( (bTb*aTa) - square2(aTb) )/aTa);
  return(cd);
}


void cross_product( double A[3], double B[3], double product[3] ) {
  
  product[2] = (A[0] * B[1]) - (A[1] * B[0]);
  product[1] = (A[2] * B[0]) - (A[0] * B[2]);
  product[0] = (A[1] * B[2]) - (A[2] * B[1]);

}




