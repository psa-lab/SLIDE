#include "defs.h"
#include "types.h"
#include "distance.h"

void  filter_points ( point_pt     points,
		      int          number_of_points,
		      pdb_atom_pt  atoms,
		      int          number_of_atoms )
{
  double dist;
  int    point_out;
  int    i, j;

  for ( i = 0; i < number_of_points; i++ )
    for ( j = 0, point_out = FALSE; 
	  j < number_of_atoms && point_out == FALSE; 
	  j++ )
      {
	dist = distance ( points[i].pos,
			  atoms[j].pos );
	if ( ( atoms[j].act == NOTHING && dist < MIN_VDW_DIST )
	     || ( atoms[j].act != NOTHING && dist < MIN_HBOND_LENGTH ))
	  {
	    point_out = TRUE;
	    points[i].type = TOO_CLOSE;
	  }
	else
	  {
	    if ( dist < MAX_HYDRO_DIST )
	      points[i].type = OK;
	  }
      }
}

void  filter_ligand_points ( point_pt    points,
			     int         number_of_points,
			     molecule_pt ligand,
			     double      threshold )
{
  int  point_in;
  int  i, j;

  for ( i = 0; i < number_of_points; i++ )
    if ( points[i].type == OK )
      {
	for ( j = 0, point_in = FALSE;
	      j < ligand->number_of_atoms && point_in == FALSE;
	      j++ )
	  if ( distance ( points[i].pos,
			  ligand->atoms[j].pos ) < threshold ) 
	    point_in = TRUE;
	if ( point_in == FALSE )
	  points[i].type = TOO_FAR;
      }
}
	    
