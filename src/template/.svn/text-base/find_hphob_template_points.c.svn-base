#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "distance.h"

/* Modified by P. Sanschagrin -- 5-Jun-00 */
/* find_hphob_template_points modified to convert from old hydrophobic
 * point criteria of aphil < 1 (100/100) to new criteria that the surface
 * point must have at least 3 hydrophobic atoms (any C attached only to C
 * or H and any S) and no more than 1 hydrophilic atoms (N or O) within 
 * 5.8 angstroms -- these modifications are noted by PCS-5Jun00
 */

int  find_hphob_template_points ( pdb_atom_pt atoms,
				  int         number_of_atoms,
				  point_pt    surface_points,
				  int         number_of_surface_points,
				  point_pt    template_points,
				  FILE        *fp)
{
  double  dist;
  float   avg_hphil;
  int     sum,
          number,
          number_of_template_points,
          too_close;
  int     i, j, k;
  int     atom_class[NUM_ATOM_CLASSES]; /*PCS-05Jun00 */


  printf("Using plusminus hphobic method\n");
  fflush ( stdout );
  fprintf ( fp,
	    "## Using plusminus hphobic method ##\n");
  number_of_template_points = 0;
  for ( i = 0; i < number_of_surface_points; i++ )
    if ( surface_points[i].type == OK )
      {
	/* Initialize the atom class array -- PCS-05Jun00 */
	for (k=0; k<NUM_ATOM_CLASSES; k++)
	  atom_class[k] = 0;
	sum = 0.0;
	number = 0;
	for ( j = 0, too_close = FALSE; 
	      j < number_of_atoms && !too_close; 
	      j++ )
	  {
	    dist = distance ( surface_points[i].pos, 
			      atoms[j].pos );
	    if ( dist < MAX_HYDRO_DIST )
	      {
		
		atom_class[atoms[j].class]++; /* PCS-05Jun00 */
		number++;
		/* Atoms with hphil < HPHIL_CUTOFF = +1, 
		   Atoms with hphil > HPHIL_CUTOFF = -1 */
		/* PCS-21Jul00 */
		if (atoms[j].hphil < HPHIL_CUTOFF) {
		  sum += HPHOB_VALUE;
		}
		else {
		  sum += HPHIL_VALUE;
		}
	      }	 
	    if ( dist < MIN_HYDRO_DIST )
	      {
		number = 0;
		too_close = TRUE;
	      }
	  }
	if ( number > 0 )
	  {
	    /* consider only those template points which are in contact with
	       at least one protein atom */
	    /* avg_hphil = (float) sum / number;
	       if ( avg_hphil < HPHIL_CUTOFF.0 ) */
	    /* New criteria for a point to be hydrophobic:
	       Atoms with hphil < HPHIL_CUTOFF = +1, Atoms with hphil 
	       > HPHIL_CUTOFF = -1 
	       Points are hydrophobic of sum of neighboring atoms is >= 3 */
	    /* PCS-21Jul00 */
	    if (sum >= HPHOB_SCORE_CUTOFF) {
	      /* the neighborhood of this point is quite hydrophobic */
	      {
		for ( j = 0; j < 3; j++ )
		  template_points[number_of_template_points].pos[j] =
		    surface_points[i].pos[j];
		number_of_template_points++;
	      }
	    }
	  }
  
      }
  return number_of_template_points;
}

