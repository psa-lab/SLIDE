
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <limits.h>
#include "defs.h"
#include "types.h"
#include "basics.h"
#include "read_mol2.h"
#include "bitstrings.h"
#include <mymalloc.h>


int  main ( int  argc, char  **argv )
{

	int i,j;
	molecule_t molecule;
	double xmin, ymin, zmin, xmax, ymax, zmax;

	xmin=ymin=zmin=DBL_MAX;
	xmax=ymax=zmax=-DBL_MAX;

/** Get all points from molecules **/
  molecule.atoms =
    (atom_pt) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (atom_t) );
  molecule.bonds =
    (bond_pt) mymalloc ( MAX_NUMBER_OF_MOL2_BONDS * sizeof (bond_t) );
  molecule.flexible_bonds =
    (int *) mymalloc ( MAX_NUMBER_OF_FLEXIBLE_BONDS * sizeof (int) );
  molecule.carbon_ring_centers =
    (float **) mymalloc ( MAX_NUMBER_OF_CARBON_RINGS * sizeof (float *) );
  molecule.carbon_ring_atom =
    (int *) mymalloc ( MAX_NUMBER_OF_CARBON_RINGS * sizeof (int) );
  for ( i = 0; i < MAX_NUMBER_OF_CARBON_RINGS; i++ )
    molecule.carbon_ring_centers[i] =
      (float *) mymalloc ( 3 * sizeof (float) );
  molecule.neighbors =
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    {
      molecule.neighbors[i] =
    (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
      molecule.atoms[i].fragments = bitstring_create ( MAX_NUMBER_OF_FLEXIBLE_BONDS );
      bitstring_clear_all ( molecule.atoms[i].fragments );
    }
  molecule.bond_to_neighbor =
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    molecule.bond_to_neighbor[i] =
      (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
  molecule.number_of_neighbors =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );


	for(i=6; i<argc;i++)
	/*for(i=1; i<argc;i++)*/
	{
		read_mol2(argv[i], &molecule);
		for(j=0; j<molecule.number_of_atoms; j++)
		{
			/** create a max-min box around it **/
			/**	printf( "%8.2f %8.2f %8.2f\n", molecule.atoms[j].pos[X], molecule.atoms[j].pos[Y], molecule.atoms[j].pos[Z]); **/

			xmin=(molecule.atoms[j].pos[X]<xmin)?molecule.atoms[j].pos[X]:xmin;
			ymin=(molecule.atoms[j].pos[Y]<ymin)?molecule.atoms[j].pos[Y]:ymin;
			zmin=(molecule.atoms[j].pos[Z]<zmin)?molecule.atoms[j].pos[Z]:zmin;

			xmax=(molecule.atoms[j].pos[X]>xmax)?molecule.atoms[j].pos[X]:xmax;
			ymax=(molecule.atoms[j].pos[Y]>ymax)?molecule.atoms[j].pos[Y]:ymax;
			zmax=(molecule.atoms[j].pos[Z]>zmax)?molecule.atoms[j].pos[Z]:zmax;
		}
	}

	printf("%8.2f %8.2f %8.2f\n", xmin-TOLERANCE, ymin-TOLERANCE, zmin-TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmax+TOLERANCE, ymin-TOLERANCE, zmin-TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmin-TOLERANCE, ymax+TOLERANCE, zmin-TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmax+TOLERANCE, ymax+TOLERANCE, zmin-TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmin-TOLERANCE, ymin-TOLERANCE, zmax+TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmax+TOLERANCE, ymin-TOLERANCE, zmax+TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmin-TOLERANCE, ymax+TOLERANCE, zmax+TOLERANCE);
	printf("%8.2f %8.2f %8.2f\n", xmax+TOLERANCE, ymax+TOLERANCE, zmax+TOLERANCE);
}
