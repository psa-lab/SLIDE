/******************************************************************************
~~~~~~~~~~~~~~~
SLIDE Licensing
~~~~~~~~~~~~~~~
                                                            
The SLIDE software (Screening for Ligands with Induced-Fit Docking,
Efficiently) was developed by Drs. Volker Schnecke and Leslie Kuhn at
Michigan State University.

Usage of this software or any part thereof is permitted only if the
SLIDE user (referred to hereafter as USER) has executed a license
agreement with Michigan State University (MSU), East Lansing, MI

If you are interested in licensing SLIDE, please contact:

Dr. Leslie Kuhn
Protein Structural Analysis and Design Lab
502C Biochemistry Building
Michigan State University
East Lansing, MI  48824-1319 USA
kuhn@agua.bch.msu.edu
(517) 353-8745 office phone
(517) 353-9334 fax

The License is by and between Michigan State University, East Lansing,
Michigan 48824 (MSU) and users of SLIDE (referred to as USER from
here).

1.  Description of Product.

As used in this agreement, Product means the full, integrated Version
2.10 of the SLIDE software for screening molecular databases for
ligands to proteins and any associated documentation, developed by
Michigan State University personnel and copyrighted by the Michigan
State University Board of Trustees.

2.  License

Neither the Product, Product copies, or derivatives from the Product
may be transferred, licensed, or sold.

3.  Acknowledgements.

USER agrees to acknowledge the use of SLIDE in publications or
presentations by citing the following references:

V. Schnecke, C. A. Swanson, E. D. Getzoff, J. A. Tainer, and L. A.
Kuhn (1998) "Screening a Peptidyl Database for Potential Ligands to
Proteins Including Side-Chain Flexibility", Proteins: Structure,
Function, and Genetics 33, 74-87.

V. Schnecke and L. A. Kuhn (2000) "Virtual Screening with Solvation
and Ligand-Induced Complementarity", Perspectives in Drug Discovery
and Design, 20, in press.

4.  Prohibited Uses of the Product.

USER may not make copies of the Product which do not contain the
notifications of copyright exactly as provided in the Product supplied
to USER by MSU.

USER may not transfer or assign its rights under this License without
the prior express written consent of MSU.

5.  Prohibited Uses of the University Name and Marks.

USER agrees that it will not use the MSU name or marks in publicity,
advertising, fund-raising, or similar activities without the prior
written approval of MSU.

6.  Intellectual Property.

Michigan State University retains title to Product. USER agrees to use
reasonable efforts to protect the Product from unauthorized use or
reproduction.  All rights not specifically granted in this License are
reserved by MSU.

7.   Warranty.

MSU MAKES NO OTHER WARRANTY, EXPRESS OR IMPLIED, TO USER OR ANY OTHER
PERSON OR ENTITY.  SPECIFICALLY, MSU MAKES NO WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OF PRODUCT.  MSU
WILL NOT BE LIABLE FOR SPECIAL, INCIDENTAL, CONSEQUENTIAL, INDIRECT OR
OTHER SIMILAR DAMAGES, EVEN IF MSU OR ITS EMPLOYEES HAVE BEEN ADVISED
OF THE POSSIBILITY OF SUCH DAMAGES.  IN NO EVENT WILL MSU LIABILITY
FOR ANY DAMAGES TO USER OR ANY PERSON EVER EXCEED THE FEE PAID FOR THE
LICENSE TO USE THE PRODUCT, REGARDLESS OF ANY FORM OF THE CLAIM.

Additional statements by employees of MSU, such as correspondence or
oral presentations, do not constitute warranties by MSU and should not
be relied upon.

8.   Supplementary Provisions.

This License represents the entire understanding and agreement between
MSU and USER regarding the Product, and supersedes any prior purchase
order, communications, advertising, or representations.  This License
may be modified only in a written amendment signed by an authorized
MSU officer.  If any provision of this License shall be unlawful,
void, or for any reason unenforceable, it shall be deemed severable
from, and shall in no way affect the validity or enforceability of the
remaining provision of this agreement.  This License shall be governed
by Michigan law.

9.   Termination.

In the event that either party hereto commits any breach of or default
in any of the terms or conditions of this Agreement, and also shall
fail to remedy such default or breach within ninety (90) days after
receipt of written notice of such breach or default, the party giving
notice may at its option and in addition to any other remedies which
it may have by law, terminate this Agreement by sending notice of
termination in writing to the other party.  Such termination shall be
effective as of the date of the receipt of such notice.

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"
#include "mymalloc.h"
#include "read_mol2.h"
#include "adj_list.h"
#include "check_connectivity.h"
#include "demo.h"

#ifdef DEMO_VERSION

#include "check_key.h"

#endif

int  main ( int  argc, char  **argv )
{
  molecule_t       molecule;
  char             slide_dir[256];
  int              i;

  if ( argc != 2 )
    {
      sprintf ( slide_dir,
		"%s <mol2-file>\n", 
		argv[0] );
      err_usage ( slide_dir );
    }
  strcpy ( slide_dir, getenv ( "SLIDE_DIR" ) );
  if ( slide_dir == NULL )
    err_panic ( "main",
                "SLIDE_DIR environment variable not set" );  

#ifdef DEMO_VERSION

	check_key(slide_dir);

#endif

  molecule.atoms = 
    (atom_pt) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (atom_t) );
  molecule.bonds = 
    (bond_pt) mymalloc ( MAX_NUMBER_OF_MOL2_BONDS * sizeof (bond_t) );
  molecule.neighbors = 
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    molecule.neighbors[i] =
      (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
  molecule.bond_to_neighbor = 
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    molecule.bond_to_neighbor[i] =
      (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
  molecule.number_of_neighbors =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );
  read_mol2 ( argv[1],
	      &molecule );
  construct_adjacency_list ( &molecule );
  if ( check_connectivity ( &molecule ) == FAILURE )
    {
      fprintf ( stderr,
		"ERROR: unconnected atoms in file %s\n",
		argv[1] );
      exit ( 42 );
    }
  exit ( EXIT_SUCCESS );
}

