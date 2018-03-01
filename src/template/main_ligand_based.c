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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "types.h"
#include <basics.h>
#include "read_mol2.h"
#include "read_flex_defn.h"
#include "read_hyd_defn.h"
#include "adj_list.h"
#include "find_hyd_atoms.h"
#include "assign_hydrogens.h"
#include "bitstrings.h"
#include "find_carbon_ring_centers.h"
#include "complete_link_clustering.h"
#include "demo.h"
#include <mymalloc.h>
#include <slide_itimer.h>

#ifdef DEMO_VERSION

#include "check_key.h"

#endif

int  main ( int  argc, char  **argv )
{
  molecule_t       molecule;
  flex_bond_defn_t flex_bond_rules[MAX_NUMBER_OF_FLEX_BOND_RULES];
  hyd_defn_t       hyd_atom_rules;
  point_t          acceptor_points[MAX_NUMBER_OF_TEMPLATE_POINTS],
                   donor_points[MAX_NUMBER_OF_TEMPLATE_POINTS],
                   doneptor_points[MAX_NUMBER_OF_TEMPLATE_POINTS],
                   hphob_points[MAX_NUMBER_OF_TEMPLATE_POINTS],
                   cluster_points[MAX_NUMBER_OF_TEMPLATE_POINTS];
  FILE             *fp_pdb,
                   *fp_tmp,
                   *fp_template,
                   *fp_log;
  char             filename[FILENAME_MAX];
  char             slide_dir[FILENAME_MAX];
  char             slide_data_dir[FILENAME_MAX];
  const char       *tmp;
  float            threshold;
  int              number_of_acceptor_points,
                   number_of_donor_points,
                   number_of_doneptor_points,
                   number_of_hphob_points,
                   number_of_cluster_points,
                   index;
  int              i, j;
  char template_dir[FILENAME_MAX]; 
  char time_buf[256];

  int mol_start_index = 4;
  const char *target_str = argv[1];
  const char *template_str = argv[2];

  if(argc < 5){
    sprintf(filename, "%s <target> <template> <clustering_threshold> "
            "<mol2-file> [<mol2-file>...]\n", argv[0]);
    err_usage(filename);
  }
  threshold = atof(argv[3]);

  tmp = getenv("SLIDE_DIR");
  if(!tmp) err_panic("main", "SLIDE_DIR environment variable not set");
  strncpy(slide_dir, tmp, FILENAME_MAX);
  slide_dir[FILENAME_MAX - 1] = '\0';
  tmp = getenv("SLIDE_DATA_DIR");
  if(!tmp) err_panic("main", "SLIDE_DATA_DIR environment variable not set");
  strncpy(slide_data_dir, tmp, FILENAME_MAX);
  slide_data_dir[FILENAME_MAX - 1] = '\0';

#ifdef DEMO_VERSION

  check_key(slide_dir);

#endif

  snprintf(template_dir, FILENAME_MAX, "%s/%s/%s/in", slide_data_dir, 
           target_str, template_str);
  snprintf(filename, FILENAME_MAX, "%s/template.pdb", template_dir); 
  fp_pdb = slide_fopen(filename, "w" );
  if(fp_pdb == NULL) err_panic ("main", "Unable to open template.pdb file");
  snprintf(filename, FILENAME_MAX, "%s/template", template_dir); 
  fp_template = slide_fopen(filename, "w" );
  if (fp_template == NULL) err_panic ("main", "Unable to open 'template' file");
  snprintf(filename, FILENAME_MAX, "%s/template.log", template_dir); 
  fp_log = slide_fopen(filename, "w" );
  if (fp_log == NULL) err_panic ("main", "Unable to open log file");

  if ( fp_pdb == NULL || fp_template == NULL || fp_log == NULL )
    /* will never get here unless there is a cosmic event */
    err_panic ( "main", "trouble opening output file" );

  fprintf(fp_template, "# SLIDE template; version: %s\n", VERSION);
  fprintf(fp_template, "# Command-Line Args:\n");
  fprintf(fp_template, "# Target :                             %s\n", 
          target_str);
  fprintf(fp_template, "# Template type :                      ligand based "
          "(biased) template\n");
  fprintf(fp_template, "# Hydrophobic clustering threshold :   %.2f (A)\n#\n",
          threshold);


  printf("\nSLIDE template; version: %s", VERSION);
  fprintf(fp_log, "\nSLIDE template; version: %s", VERSION);
  printf("\nCommand-Line Args: \n");
  fprintf(fp_log, "\nCommand-Line Args: \n");
  for ( i = 0; i < argc; i++ ){
    printf("\t%s\n", argv[i]);
    fprintf(fp_log, "\t%s\n", argv[i]);
  }

  /* Add the timing if desired at a different date */
  /* if(slide_itimer_start(10)) return -1; */

  if(slide_get_local_time(time_buf, 256)){
    printf("Local start time: %s\n", time_buf );
    fprintf(fp_log, "Local start time: %s\n", time_buf );
    fprintf(fp_template, "# Local start time: %s\n#\n", time_buf );
  }

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
      molecule.atoms[i].fragments = 
	bitstring_create ( MAX_NUMBER_OF_FLEXIBLE_BONDS );
      bitstring_clear_all ( molecule.atoms[i].fragments );
    }
  molecule.bond_to_neighbor = 
    (int **) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int *) );
  for ( i = 0; i < MAX_NUMBER_OF_MOL2_ATOMS; i++ )
    molecule.bond_to_neighbor[i] =
      (int *) mymalloc ( MAX_NEIGHBORS * sizeof (int) );
  molecule.number_of_neighbors =
    (int *) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (int) );

  sprintf(filename, "%s/params/hbond.defn", slide_dir );
  read_hyd_defn(filename, &hyd_atom_rules );
  sprintf(filename, "%s/params/flex.defn", slide_dir );
  read_flex_defn(filename, &flex_bond_rules[0] );

  number_of_acceptor_points = 0;
  number_of_donor_points = 0;
  number_of_doneptor_points = 0;
  number_of_hphob_points = 0;
  fprintf(fp_template, "# Ligands:\n");
  for(index = mol_start_index; index < argc; index++){
    /* Note should be a check here in the event that the file name(s) might
     * have been misspelled */

    fprintf(fp_template, "#   %s\n", argv[index]);

    read_mol2(argv[index], &molecule);
    construct_adjacency_list(&molecule);
    find_carbon_ring_centers(&molecule);
    assign_hydrogens(&molecule);  
    find_hyd_atoms(&molecule, &hyd_atom_rules);
      
      for ( i = 0; i < molecule.number_of_atoms; i++ )
	switch ( molecule.atoms[i].hyd)
	  {
	  case DONOR:
	    donor_points[number_of_donor_points].type = DONOR;
	    for ( j = 0; j < 3; j++ )
	      donor_points[number_of_donor_points].pos[j] = 
		molecule.atoms[i].pos[j];
	    number_of_donor_points++;
	    break;
	  case ACCEPTOR:
	    acceptor_points[number_of_acceptor_points].type = ACCEPTOR;
	    for ( j = 0; j < 3; j++ )
	      acceptor_points[number_of_acceptor_points].pos[j] = 
		molecule.atoms[i].pos[j];
	    number_of_acceptor_points++;
	    break;
	  case DONEPTOR:
	    doneptor_points[number_of_doneptor_points].type = DONEPTOR;
	    for ( j = 0; j < 3; j++ )
	      doneptor_points[number_of_doneptor_points].pos[j] = 
		molecule.atoms[i].pos[j];
	    number_of_doneptor_points++;
	    break;
	  default:
	    break;
	  }
      for ( i = 0; i < molecule.number_of_carbon_rings; i++ )
	{
	  hphob_points[number_of_hphob_points].type = HPHOB;
	  for ( j = 0; j < 3; j++ )
	    hphob_points[number_of_hphob_points].pos[j] =
	      molecule.carbon_ring_centers[i][j];
	  number_of_hphob_points++;
	}
    }  
  printf ( "%d hydrophobic points\n",
           number_of_hphob_points );
  fprintf ( fp_log,
	    "%d hydrophobic points\n",
	    number_of_hphob_points );
  printf ( "%d initial h-bonding points (%d A, %d D, %d N)\n",
           number_of_acceptor_points
           + number_of_donor_points
           + number_of_doneptor_points,
           number_of_acceptor_points,
           number_of_donor_points,
           number_of_doneptor_points );
  fflush ( stdout );
  fprintf ( fp_log,
	    "%d initial h-bonding points (%d A, %d D, %d N)\n",
	    number_of_acceptor_points
	    + number_of_donor_points
	    + number_of_doneptor_points,
	    number_of_acceptor_points,
	    number_of_donor_points,
	    number_of_doneptor_points );
  number_of_cluster_points =
    complete_link_clustering(hphob_points, number_of_hphob_points,
                             cluster_points, threshold );
  printf("%d hydrophobic cluster points\n", number_of_cluster_points );
  fflush(stdout );
  fprintf(fp_log, "%d hydrophobic cluster points\n", number_of_cluster_points );

  snprintf(filename, FILENAME_MAX, "%s/hphob.pdb", template_dir); 
  fp_tmp = slide_fopen(filename, "w" );
  if(fp_tmp == NULL) err_panic ("main", "Unable to open hphob.pdb file");

  for ( i = 0; i < number_of_cluster_points; i++ )
    {
      fprintf ( fp_pdb, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00100.00\n", 
                i + 1,
                i + 1,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_tmp, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00100.00\n", 
                i + 1,
                i + 1,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_template, 
                "H*  %7.3f %7.3f %7.3f\n", 
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
    }
  fclose ( fp_tmp );
  number_of_cluster_points =
    complete_link_clustering ( acceptor_points,
                               number_of_acceptor_points,
                               cluster_points,
                               threshold );
  printf("%d acceptor cluster points\n", number_of_cluster_points );
  fflush ( stdout );
  fprintf ( fp_log, "%d acceptor cluster points\n", number_of_cluster_points );

  snprintf(filename, FILENAME_MAX, "%s/acceptor.pdb", template_dir); 
  fp_tmp = slide_fopen(filename, "w" );
  if(fp_tmp == NULL) err_panic ("main", "Unable to open acceptor.pdb file");

  for ( i = 0; i < number_of_cluster_points; i++ )
    {
      fprintf ( fp_pdb, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
                1000 + i,
                1000 + i,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_tmp, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
                i + 1,
                i + 1,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_template, 
                "A*  %7.3f %7.3f %7.3f\n", 
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
    }
  fclose ( fp_tmp );
  number_of_cluster_points =
    complete_link_clustering ( donor_points,
                               number_of_donor_points,
                               cluster_points,
                               threshold );
  printf ( "%d donor cluster points\n",
           number_of_cluster_points );
  fflush ( stdout );
  fprintf ( fp_log, "%d donor cluster points\n", number_of_cluster_points );

  snprintf(filename, FILENAME_MAX, "%s/donor.pdb", template_dir); 
  fp_tmp = slide_fopen(filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open donor.pdb file");

  for ( i = 0; i < number_of_cluster_points; i++ )
    {
      fprintf ( fp_pdb, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 50.00\n", 
                2000 + i,
                2000 + i,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_tmp, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 50.00\n", 
                i + 1,
                i + 1,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_template, 
                "D*  %7.3f %7.3f %7.3f\n", 
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
    }
  fclose ( fp_tmp );
  number_of_cluster_points =
    complete_link_clustering ( doneptor_points,
                               number_of_doneptor_points,
                               cluster_points,
                               threshold );
  printf ( "%d doneptor cluster points\n", number_of_cluster_points );
  fflush ( stdout );
  fprintf ( fp_log, "%d doneptor cluster points\n", number_of_cluster_points );

  snprintf(filename, FILENAME_MAX, "%s/doneptor.pdb", template_dir); 
  fp_tmp = slide_fopen(filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open doneptor.pdb file");
  for ( i = 0; i < number_of_cluster_points; i++ )
    {
      fprintf ( fp_pdb, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 25.00\n", 
                3000 + i,
                3000 + i,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_tmp, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 25.00\n", 
                i + 1,
                i + 1,
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
      fprintf ( fp_template, 
                "N*  %7.3f %7.3f %7.3f\n", 
                cluster_points[i].pos[X],
                cluster_points[i].pos[Y],
                cluster_points[i].pos[Z] );
    }
  fclose ( fp_tmp );
  fclose ( fp_pdb );
  fclose ( fp_template );
  fclose ( fp_log );
  exit ( 0 );
}


