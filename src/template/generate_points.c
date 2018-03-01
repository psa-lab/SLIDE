#include <stdio.h>
#include <stdlib.h>
/*#include <sys/times.h> */
#include <time.h>
#include "defs.h"
#include "types.h"
#include "basics.h"
#include <mymalloc.h>

#define TRAC

void  generate_points ( point_pt surface_points,
			int      number_of_surface_points,
			point_pt center_points,
			int      number_of_center_points,
			FILE     *fp )
{
  double  min[3],
          max[3];
  double  density,
          volume;
  int     interval[3];
  int     i, j;

  min[X] = min[Y] = min[Z] = 99999.9;
  max[X] = max[Y] = max[Z] = -99999.9;
  for ( i = 0; i < number_of_center_points; i++ )
    for ( j = 0; j < 3; j++ )
      {
	if ( center_points[i].pos[j] < min[j] )
	  min[j] = center_points[i].pos[j];
	if ( center_points[i].pos[j] > max[j] )
	  max[j] = center_points[i].pos[j];
      }
  volume = 1.0;
  for ( i = 0; i < 3; i++ )
    {
      interval[i] = (int) ( 1000 * ( max[i] - min[i] ) );
      volume *= max[i] - min[i];
    }
  density = number_of_surface_points / volume;
  printf ( "%d random points generated\n", 
	   number_of_surface_points );
  fprintf ( fp,
	    "%d random points generated\n", 
	    number_of_surface_points );
  printf ( "%5.2f points per A^3\n", 
	   density );
  fflush ( stdout );
  fprintf ( fp, 
	    "%5.0f points per A^3\n", 
	    density );
  srandom ( (int) time ( 0 ) );
  for ( i = 0; i < number_of_surface_points; i++ )
    {
      for ( j = 0; j < 3; j++ )
	surface_points[i].pos[j] = 
	  min[j] + ( (double) ( random ( ) % interval[j] ) / 1000 );
      surface_points[i].type = NOTHING;
    }
}

point_pt  generate_grid_points ( int      *number_of_surface_points,
				 point_pt center_points,
				 int      number_of_center_points,
				 double   grid_spacing,
				 FILE     *fp )
{
  point_pt  surface_points;
  double    min[3],
            max[3],
            x, y, z;
  int       i, j;
/* ADDED BY PCS -- 05-Apr-2000 */
  double temp;
/*******************************/

  min[X] = min[Y] = min[Z] = 99999.9;
  max[X] = max[Y] = max[Z] = -99999.9;

for ( i = 0; i < number_of_center_points; i++ ) {
/* ADDED BY PCS -- 05-Apr-2000 */
#ifdef TRACE
    printf("C Pt %3d: %5.2f  %5.2f  %5.2f\n",i+1,center_points[i].pos[X],
	   center_points[i].pos[Y],center_points[i].pos[Z]);
#endif
/*******************************/
    for ( j = 0; j < 3; j++ )
      {
	if ( center_points[i].pos[j] < min[j] )
	  min[j] = center_points[i].pos[j];
	if ( center_points[i].pos[j] > max[j] )
	  max[j] = center_points[i].pos[j];
      }
  }
/* ADDED BY PCS -- 05-Apr-2000 */
#ifdef TRACE
    for (i=0;i<=Z;i++) {
      printf("%c:  %5.2f  %5.2f\n",i+88,min[i],max[i]);
    }
#endif
    /* For some reason, setting number_of_surface_points directly to the typecast
       function gives an erroneous result. Instead, a the results is assigned to
       a temporary double, which is then typecast to assign to n_o_s_p, which
       gives the correct number -- PCS, 05-Apr-00 */
    temp = ( ( max[0] - min[0] + grid_spacing ) / grid_spacing )  
    * ( ( max[1] - min[1] + grid_spacing ) / grid_spacing ) 
    * ( ( max[2] - min[2] + grid_spacing ) / grid_spacing );
  *number_of_surface_points = (int) temp;
  /*    (int) ( ( max[0] - min[0] + grid_spacing ) / grid_spacing )  
    * ( ( max[1] - min[1] + grid_spacing ) / grid_spacing ) 
    * ( ( max[2] - min[2] + grid_spacing ) / grid_spacing );*/
/*******************************/  
surface_points = 
    (point_pt) mymalloc ( *number_of_surface_points * sizeof (point_t) );
  

  i = 0;
  /** changed < to <= -- RSK (may 30, 2002) **/
  for ( x = min[0]; x <= max[0]; x += grid_spacing )
    for ( y = min[1]; y <= max[1]; y += grid_spacing )
      for ( z = min[2]; z <= max[2]; z += grid_spacing )
	{
	  if (i < *number_of_surface_points) {
	  surface_points[i].pos[0] = x;
	  surface_points[i].pos[1] = y;
	  surface_points[i].pos[2] = z;
	  surface_points[i].type = NOTHING;
	  i++;
	}}
  
  /** added by RSK (may 30, 2002) to eliminate extra points **/
  *number_of_surface_points = i;
  /*  printf("number of surface points is %d\n", i);*/

  printf ( "%d grid points generated\n", 
	   *number_of_surface_points);
  fprintf ( fp,
	    "%d grid points generated\n", 
	    *number_of_surface_points );
  printf ( "%5.3f A grid_spacing\n", 
	   grid_spacing );
  fflush ( stdout );
  return surface_points;
}
