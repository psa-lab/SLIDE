#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "math_util.h"

/**********************************************************************/
/* Given three initial vectors, and three final vectors, this routine */
/* will determine the required rotation matrices to overlay the two   */
/* sets of vectors. Requires that the vectors have been translated to */
/* the origin.                                                        */
/**********************************************************************/
void get_rotation_matrix( double X[3], double Y[3], double Z[3],
			  double rot_matrix[3][3] ) {

  int 
    a = 0; /* counters */

  double 
    Origin[3],                               /* coords for the origin */
    mag_X = 0.0, mag_Y = 0.0, mag_Z = 0.0;   /* magnitude of incoming vectors */

  /****************************************/
  /****************************************/
  
  for( a = 0; a < 3; a++ )
    Origin[a] = 0.0;
  
  /****************************************/
  /* Calculate the magnitude of each vect.*/
  /****************************************/
  mag_X = distance2( X, Origin );
  mag_Y = distance2( Y, Origin );
  mag_Z = distance2( Z, Origin );

  /****************************************/
  /* Normalize all the coordinates.       */
  /****************************************/
  for( a = 0; a < 3; a++ ) {
    X[a] = X[a] / mag_X;
    Y[a] = Y[a] / mag_Y;
    Z[a] = Z[a] / mag_Z;
  }

  for( a = 0; a < 3; a++ ){
    rot_matrix[a][0] = X[a];
    rot_matrix[a][1] = Y[a];
    rot_matrix[a][2] = Z[a];
  }
  
}






