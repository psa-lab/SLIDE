/*
 *
 *    refind.c     Litian He     06/15/2003
 *
 *    refine the coordinates that are too small
 *    
 */

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "defs.h"

void  refine ( molecule_pt molecule )
{
  int     i, j;

  for ( i = 0; i < molecule->number_of_atoms; i++ )
      {
	  if ( molecule->atoms[i].hyd != NOTHING )
	      {
		  for ( j = 0; j < 3; j ++ )
		      if ( fabs( molecule->atoms[i].pos[j] ) < 0.0001 )
			  molecule->atoms[i].pos[j] = 0.0;
	      }
      }

  for ( i = 0; i < molecule->number_of_carbon_rings; i++ )
      {
	  for ( j = 0; j < 3; j ++ )
	      if ( fabs( molecule->carbon_ring_centers[i][j] ) < 0.0001 )
		  molecule->carbon_ring_centers[i][j] = 0.0;
      }
}
