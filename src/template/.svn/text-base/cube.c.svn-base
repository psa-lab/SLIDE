/* Functions for checking if a point is within a cube */
/*          -- RSK                                    */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

/* check which side of the plane, a point X, Y, Z is */
/* Plane is specified by three points P1, P2, and P3 */
/* which are index into the NBox array.              */

int point_on_plane(double Nbox[8][3], int p1, int p2, int p3, double x, double y, double z)
{
	double Ac, Bc, Cc, Dc;
	double s;

	Ac=Nbox[p1][1]*(Nbox[p2][2]-Nbox[p3][2]) + Nbox[p2][1]*(Nbox[p3][2]-Nbox[p1][2]) + Nbox[p3][1]*(Nbox[p1][2]-Nbox[p2][2]);
	Bc=Nbox[p1][2]*(Nbox[p2][0]-Nbox[p3][0]) + Nbox[p2][2]*(Nbox[p3][0]-Nbox[p1][0]) + Nbox[p3][2]*(Nbox[p1][0]-Nbox[p2][0]);
	Cc=Nbox[p1][0]*(Nbox[p2][1]-Nbox[p3][1]) + Nbox[p2][0]*(Nbox[p3][1]-Nbox[p1][1]) + Nbox[p3][0]*(Nbox[p1][1]-Nbox[p2][1]);
	Dc=-(Nbox[p1][0]*(Nbox[p2][1]*Nbox[p3][2]-Nbox[p3][1]*Nbox[p2][2]) + Nbox[p2][0]*(Nbox[p3][1]*Nbox[p1][2]-Nbox[p1][1]*Nbox[p3][2]) + Nbox[p3][0]*(Nbox[p1][1]*Nbox[p2][2]-Nbox[p2][1]*Nbox[p1][2]));

	s =  (Ac*x + Bc*y + Cc*z +Dc);
	if(s<=-1 ) return -1;
	else if(s>=1) return 1;
	else return 0;
}

/** check if a point is inside a cube                                   **/
/** checks if a point (x,y,z) is on the correct side of all the planes. **/

int point_in_cube(double Nbox[8][3], double x, double y, double z)
{
	if((point_on_plane(Nbox, 0, 1, 2, Nbox[6][0], Nbox[6][1], Nbox[6][2]) == point_on_plane(Nbox, 0, 1, 2, x, y, z)) &&
	    (point_on_plane(Nbox, 4, 5, 6, Nbox[2][0], Nbox[2][1], Nbox[2][2]) == point_on_plane(Nbox, 4, 5, 6, x, y, z)) &&
	    (point_on_plane(Nbox, 0, 2, 4, Nbox[7][0], Nbox[7][1], Nbox[7][2]) == point_on_plane(Nbox, 0, 2, 4, x, y, z)) &&
	    (point_on_plane(Nbox, 3, 5, 7, Nbox[2][0], Nbox[2][1], Nbox[2][2]) == point_on_plane(Nbox, 3, 5, 7, x, y, z)) &&
	    (point_on_plane(Nbox, 2, 3, 7, Nbox[5][0], Nbox[5][1], Nbox[5][2]) == point_on_plane(Nbox, 2, 3, 7, x, y, z)) &&
	    (point_on_plane(Nbox, 0, 1, 4, Nbox[2][0], Nbox[2][1], Nbox[2][2]) == point_on_plane(Nbox, 0, 1, 4, x, y, z)))
		return 1;
	else return 0;
}
