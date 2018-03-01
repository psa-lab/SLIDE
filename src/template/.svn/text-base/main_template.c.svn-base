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
#include <ctype.h>
#include "defs.h"
#include "types.h"
#include <basics.h>
#include <mymalloc.h>
#include "read_pdb_file.h"
#include "read_mol2.h"
#include "read_borders_file.h"
#include "nbh.h"
#include "find_hbond_template_points.h"
#include "find_hphob_template_points.h"
#include "assign_hphil.h"
#include "generate_points.h"
#include "filter_points.h"
#include "complete_link_clustering.h"
/* ADDED BY PCS -- 05-Apr-2000 */
#include "math_util.h"
#include "transform_points.h"
/* ADDED BY PCS -- 05-Jun-2000 */
#include "assign_atom_class.h"
#include "demo.h"
#include <slide_itimer.h>

#ifdef DEMO_VERSION

#include "check_key.h"

#endif

void print_template_header(FILE* template, FILE* log, const char* target_str, 
                           float grid_spacing, float hphob_threshold, 
                           const int argc, const char** argv);

#define TRAC
#undef PDB_TEST
/*******************************/

int  main (const int  argc, const char  **argv )
{
  FILE        *fp,
              *fp_tmp,
              *fp_log,
              *fp_pdb,
              *fp_template,
/* ADDED BY PCS -- 05-Apr-2000 */
              *fp_test;
/*******************************/
  pdb_atom_pt atoms;
  point_pt    surface_points,
              center_points,
              ut_center_points,
              hphob_template_points,
              metal_template_points,
              acceptor_template_points,
              donor_template_points,
              doneptor_template_points,
              hphob_cluster_points,
              hbond_cluster_points;
  point_t     additional_donors[64],
              additional_acceptors[64];
  molecule_t  molecule;
  char        filename[FILENAME_MAX];
  char        pdbfilename[FILENAME_MAX];
  char        slide_dir[FILENAME_MAX];
  double      threshold, dist,
/* ADDED BY PCS -- 05-Apr-2000 */
              hphob_threshold,
              cutoff,
              grid_spacing,
              origin_point[3];
  short       hbond_pt_density;
/*******************************/
  int         number_of_atoms,
              number_of_surface_points,
              number_of_center_points,
              number_of_hphob_template_points,
              number_of_metal_template_points,
              number_of_acceptor_template_points,
              number_of_donor_template_points,
              number_of_doneptor_template_points,
              number_of_hphob_cluster_points,
              number_of_hbond_cluster_points,
              number_of_additional_donors,
              number_of_additional_acceptors,
              number_of_deleted_points;
  int         i, j,
/* ADDED BY PCS -- 05-Apr-2000 */
              boxflag = FALSE;
/*******************************/

/* ADDED BY PCS -- 05-Apr-2000 */
  double      transform_vectors[3][3];
/*******************************/
  const char* tmp;
  char template_dir[FILENAME_MAX];
  char slide_data_dir[FILENAME_MAX];
  const char *target_str = argv[1];
  const char *template_str = argv[2];

  number_of_surface_points = 0;
  if ( argc != 6 && argc != 8 ){
    printf("Usage: %s <target> <template> <dense|sparse|minimal> <grid_spacing>"
           " <hphobic clustering threshold> [<mol2-reference> <dist>]\n", 
           argv[0] );
    exit ( 1 );
  }

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

  /* Parse CL Args */ 
  grid_spacing = atof ( argv[SPACING_CLARG_POS] ); 
  hphob_threshold = atof ( argv[HPHOBTHRESH_CLARG_POS] );
  /** threshold is used for bump-checking hbond and hphob points
	in case of a bump, hbond point prevails **/
  threshold = hphob_threshold/2;

  if(toupper(argv[HBONDTHRESH_CLARG_POS][0])=='D') hbond_pt_density=DENSE;	
  else if(toupper(argv[HBONDTHRESH_CLARG_POS][0])=='S') hbond_pt_density=SPARSE;
  else if(toupper(argv[HBONDTHRESH_CLARG_POS][0])=='M') 
    hbond_pt_density=MINIMAL;	
  else{
    fprintf(stderr,"ERROR: Unknown HBond density type: %s\n", 
            argv[HBONDTHRESH_CLARG_POS]); 
    exit(1); 
  }

  /* open files and print header information */
  printf("opening files\n");
  snprintf(template_dir, FILENAME_MAX, "%s/%s/%s/in", slide_data_dir,
           target_str, template_str);
  snprintf(filename, FILENAME_MAX, "%s/template.log", template_dir);
  fp_log = slide_fopen(filename, "w" );
  if(fp_log == NULL) err_panic ("main", "Unable to open log file") ;
  snprintf(filename, FILENAME_MAX, "%s/template", template_dir);
  fp_template = slide_fopen ( filename, "w" );
  if (fp == NULL) err_panic ( "main", "Unable to open 'template' file");

  print_template_header(fp_template, fp_log, target_str, grid_spacing,
                        hphob_threshold, argc, argv);

  atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
  if(argc == 8){
    molecule.atoms = 
      (atom_pt) mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (atom_t) );
    molecule.bonds = 
      (bond_pt) mymalloc ( MAX_NUMBER_OF_MOL2_BONDS * sizeof (bond_t) );
  }
  snprintf(pdbfilename, FILENAME_MAX, "%s/%s/%s.pdb", slide_data_dir,
           target_str, target_str);
  number_of_atoms = read_pdb_file ( pdbfilename, atoms );
  assign_hphil ( atoms, number_of_atoms );    
  assign_atom_class ( atoms, number_of_atoms );  /* PCS-06Jun00 */
  center_points = (point_pt) mymalloc ( MAX_NUMBER_OF_CENTER_POINTS
					* sizeof (point_t) );
  ut_center_points = (point_pt) mymalloc ( MAX_NUMBER_OF_CENTER_POINTS
					* sizeof (point_t) );

  snprintf(filename, FILENAME_MAX, "%s/borders.xyz", template_dir);
  number_of_center_points = read_borders_file(filename, center_points );
  for(i=0;i<number_of_center_points;i++)
  { 
    ut_center_points[i].pos[X] = center_points[i].pos[X];
    ut_center_points[i].pos[Y] = center_points[i].pos[Y];
    ut_center_points[i].pos[Z] = center_points[i].pos[Z];
  }

/* ADDED BY PCS -- 05-Apr-2000 */
  if (number_of_center_points == 8) {
    printf("8 border points found. Running in 'box' mode.\n");
    boxflag = TRUE;
#ifdef PDB_TEST
    fp_test = fopen ( "test.pdb", "w");
    if (fp_test == NULL) err_panic ("main", "Unable to open test.pdb file");
    for (i=0;i<8;i++) {
      fprintf(fp_test,"HETATM%5d BOX  BOX %4d     %8.3f%8.3f%8.3f\n",
	      i+1,i+1,center_points[i].pos[X],
	      center_points[i].pos[Y],
	      center_points[i].pos[Z]);
    }
#endif  

    transform_center_pts (center_points, number_of_center_points,
			  origin_point, transform_vectors);

#ifdef TRACE
    /* Check transformation vectors */
    printf("Vec1 (P2P1): ");
    for (i=0;i<=Z;i++)
      printf("+ %7.3f%c ",transform_vectors[0][i],i+120);
    printf("\nVec2 (P3P1): ");
    for (i=0;i<=Z;i++)
      printf("+ %7.3f%c ",transform_vectors[1][i],i+120);
    printf("\nVec3 (P5P1): ");
    for (i=0;i<=Z;i++)
      printf("+ %7.3f%c ",transform_vectors[2][i],i+120);
    printf("\n");
    /* check new cornerpts */
    for(i=0;i<number_of_center_points;i++) {
      printf("C Pt %2d: ",i+1);
      for (j=0;j<=Z;j++) {
	printf ("%8.3f ",center_points[i].pos[j]);
      }
      printf("\n");
    }
#endif
  }
/*******************************/

  if ( argc == 6 )
    fprintf(fp_log, "%s %s %s %s %s %s\n", argv[0], argv[TARGET_CLARG_POS], 
            argv[TEMPLATE_CLARG_POS], argv[HBONDTHRESH_CLARG_POS], 
	    argv[SPACING_CLARG_POS], argv[HPHOBTHRESH_CLARG_POS] );
  else
    fprintf(fp_log, "%s %s %s %s %s %s %s %s\n", argv[0], 
            argv[TARGET_CLARG_POS], argv[TEMPLATE_CLARG_POS], 
            argv[HBONDTHRESH_CLARG_POS], argv[SPACING_CLARG_POS],
	    argv[HPHOBTHRESH_CLARG_POS], argv[MOL2REF_CLARG_POS], 
	    argv[DIST_CLARG_POS] ); 

  surface_points = 
    generate_grid_points(&number_of_surface_points, center_points, 
                         number_of_center_points, grid_spacing, fp_log);

/* ADDED BY PCS -- 05-Apr-2000 */
  if (boxflag) {
    fprintf ( fp_log, "Running in box mode.\n");
/* Transform generated points away from axes */ 
    transform_surface_pts (surface_points, number_of_surface_points,
			   origin_point, transform_vectors);

#ifdef PDB_TEST
    for (i=0;i<number_of_surface_points;i++) {
      fprintf(fp_test,"HETATM%5d PTS  PTS %4d     %8.3f%8.3f%8.3f\n",
	      i+9,i+9,surface_points[i].pos[X],
	      surface_points[i].pos[Y],
	      surface_points[i].pos[Z]);
    }
    fprintf(fp_test,"CONECT    1    2\n");
    fprintf(fp_test,"CONECT    1    3\n");
    fprintf(fp_test,"CONECT    1    4\n");
    fprintf(fp_test,"CONECT    1    5\n");
    fprintf(fp_test,"CONECT    1    6\n");
    fprintf(fp_test,"CONECT    1    7\n");
    fprintf(fp_test,"CONECT    2    3\n");
    fprintf(fp_test,"CONECT    2    4\n");
    fprintf(fp_test,"CONECT    2    5\n");
    fprintf(fp_test,"CONECT    2    6\n");
    fprintf(fp_test,"CONECT    2    8\n");
    fprintf(fp_test,"CONECT    3    4\n");
    fprintf(fp_test,"CONECT    3    5\n");
    fprintf(fp_test,"CONECT    3    7\n");
    fprintf(fp_test,"CONECT    3    8\n");
    fprintf(fp_test,"CONECT    4    6\n");
    fprintf(fp_test,"CONECT    4    7\n");
    fprintf(fp_test,"CONECT    4    8\n");
    fprintf(fp_test,"CONECT    5    6\n");
    fprintf(fp_test,"CONECT    5    7\n");
    fprintf(fp_test,"CONECT    5    8\n");
    fprintf(fp_test,"CONECT    6    7\n");
    fprintf(fp_test,"CONECT    6    8\n");
    fprintf(fp_test,"CONECT    7    8\n");
#endif
  }
/*******************************/  

  filter_points(surface_points, number_of_surface_points, atoms, 
                number_of_atoms);

  if( argc == 8){
    read_mol2(argv[MOL2REF_CLARG_POS], &molecule );
    cutoff = atof ( argv[DIST_CLARG_POS] );
    filter_ligand_points(surface_points, number_of_surface_points, &molecule, 
                         cutoff);
  }
  snprintf(filename, FILENAME_MAX, "%s/points.usr", template_dir);
  fp = slide_fopen(filename, "w");
  if(fp == NULL) err_panic ( "main", "Unable to open usr file");

  fprintf ( fp, "DOTS\n" );
  j = 0;
  for ( i = 0; i < number_of_surface_points; i++ )
    if ( surface_points[i].type == OK){ 
      fprintf(fp, "%8.3f %8.3f %8.3f\n", surface_points[i].pos[X], 
              surface_points[i].pos[Y], surface_points[i].pos[Z] );
      j++; 
    }

  fclose ( fp );
  fprintf(fp_log, "%d close surface points\n", j );
  printf("%d close surface points\n", j );
  if(j <= 0)
    err_panic("main", "No protein pdb atoms found inside borders.xyz box");
  fflush ( stdout );

  metal_template_points =
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  acceptor_template_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  donor_template_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  doneptor_template_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  hphob_template_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  number_of_hphob_template_points =
    find_hphob_template_points(atoms, number_of_atoms, surface_points, 
                               number_of_surface_points, hphob_template_points,
                               fp_log );
  printf ( "%d hydrophobic surface points\n", number_of_hphob_template_points );
  fflush ( stdout );
  fprintf(fp_log, "%d hydrophobic surface points\n", 
          number_of_hphob_template_points );

  /*  find_metal_template_points ( atoms,
			       number_of_atoms,
			       surface_points,
			       number_of_surface_points,
			       metal_template_points,
			       number_of_metal_template_points );*/

  find_hbond_template_points ( pdbfilename,
			       surface_points,
			       number_of_surface_points,
			       ut_center_points,
			       number_of_center_points,
			       acceptor_template_points,
			       &number_of_acceptor_template_points,
			       donor_template_points,
			       &number_of_donor_template_points, 
			       doneptor_template_points,
			       &number_of_doneptor_template_points,
                   hbond_pt_density);
  printf ( "%d H-bonding surface points (%d A, %d D, %d N)\n",
	   number_of_acceptor_template_points
	   + number_of_donor_template_points
	   + number_of_doneptor_template_points,
	   number_of_acceptor_template_points,
	   number_of_donor_template_points,
	   number_of_doneptor_template_points );
  fflush ( stdout );
  fprintf ( fp_log,
	    "%d H-bonding surface points (%d A, %d D, %d N)\n",
	    number_of_acceptor_template_points
	    + number_of_donor_template_points
	    + number_of_doneptor_template_points,
	    number_of_acceptor_template_points,
	    number_of_donor_template_points,
	    number_of_doneptor_template_points );

  hbond_cluster_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );
  hphob_cluster_points = 
    (point_pt ) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS
			   * sizeof (point_t) );

  /** Cluster HPhobic points **/
  number_of_hphob_cluster_points =
    complete_link_clustering ( hphob_template_points,
			       number_of_hphob_template_points,
			       hphob_cluster_points,
			       hphob_threshold );
  printf("%d hydrophobic cluster points\n", number_of_hphob_cluster_points);
  fprintf(fp_log, "%d hydrophobic cluster points\n", 
          number_of_hphob_cluster_points );
  for(i=0; i<number_of_hphob_cluster_points; i++)
    hphob_cluster_points[i].type=HPHOB;

  snprintf(filename, FILENAME_MAX, "%s/template.pdb", template_dir);
  fp_pdb = slide_fopen(filename, "w");
  if (fp == NULL) err_panic ( "main", "Unable to open 'template.pdb' file"); 

  number_of_hbond_cluster_points = 0;
  number_of_deleted_points = 0;

  /** Acceptor Output **/
  printf("%d acceptor cluster points\n", number_of_acceptor_template_points);
  fprintf(fp_log, "%d acceptor cluster points\n", 
          number_of_acceptor_template_points);

  snprintf(filename, FILENAME_MAX, "%s/acceptor.pdb", template_dir);
  fp_tmp = slide_fopen ( filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open tmp file");
  snprintf(filename, FILENAME_MAX, "%s/acceptor.usr", template_dir);
  fp = slide_fopen ( filename, "w" );
  if (fp == NULL) err_panic ( "main", "Unable to open usr file");
  fprintf ( fp, "DOTS\n" );

  for ( i = 0; i < number_of_acceptor_template_points; i++ )
  {
     j=number_of_hbond_cluster_points;
     hbond_cluster_points[j].pos[X] = acceptor_template_points[i].pos[X];
     hbond_cluster_points[j].pos[Y] = acceptor_template_points[i].pos[Y];
     hbond_cluster_points[j].pos[Z] = acceptor_template_points[i].pos[Z];
     hbond_cluster_points[j].type = ACCEPTOR; 
     number_of_hbond_cluster_points++;

    fprintf ( fp,
	      "%8.3f %8.3f %8.3f\n",
	      acceptor_template_points[i].pos[X],
	      acceptor_template_points[i].pos[Y],
	      acceptor_template_points[i].pos[Z] );

	fprintf ( fp_tmp, 
		  "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
		  i + 1,
		  i + 1,
		  acceptor_template_points[i].pos[X],
		  acceptor_template_points[i].pos[Y],
		  acceptor_template_points[i].pos[Z] );
	fprintf ( fp_pdb, 
		  "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  0.00\n", 
		  1000 + i,
		  1000 + i,
		  acceptor_template_points[i].pos[X],
		  acceptor_template_points[i].pos[Y],
		  acceptor_template_points[i].pos[Z] );
	fprintf ( fp_template, 
		  "A*  %7.3f %7.3f %7.3f\n", 
		  acceptor_template_points[i].pos[X],
		  acceptor_template_points[i].pos[Y],
		  acceptor_template_points[i].pos[Z] );

  	for ( j = 0; j < number_of_hphob_cluster_points; j++ )
    {
      if(hphob_cluster_points[j].type==HPHOB)
      {
        dist = find_dist(acceptor_template_points[i].pos[X], 
                       acceptor_template_points[i].pos[Y], 
                       acceptor_template_points[i].pos[Z], 
                       hphob_cluster_points[j].pos[X], 
                       hphob_cluster_points[j].pos[Y], 
                       hphob_cluster_points[j].pos[Z]);
        if(dist<threshold) 
        { 
            hphob_cluster_points[j].type=NOTHING;
            number_of_deleted_points++; 
            break; 
        }
      }
	}
  }
  fclose ( fp_tmp );
  fclose ( fp );

  /** Donor Output **/
  printf ( "%d donor cluster points\n", number_of_donor_template_points);
  fprintf(fp_log, "%d donor cluster points\n", number_of_donor_template_points);
  snprintf(filename, FILENAME_MAX, "%s/donor.pdb", template_dir);
  fp_tmp = slide_fopen ( filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open tmp file");
  snprintf(filename, FILENAME_MAX, "%s/donor.usr", template_dir);
  fp = slide_fopen ( filename, "w" );
  if (fp == NULL) err_panic ( "main", "Unable to open usr file");
  fprintf ( fp, "DOTS\n" );

  for ( i = 0; i < number_of_donor_template_points; i++ )
  {
     j=number_of_hbond_cluster_points;
     hbond_cluster_points[j].pos[X] = donor_template_points[i].pos[X];
     hbond_cluster_points[j].pos[Y] = donor_template_points[i].pos[Y];
     hbond_cluster_points[j].pos[Z] = donor_template_points[i].pos[Z];
     hbond_cluster_points[j].type = DONOR; 
     number_of_hbond_cluster_points++;

    fprintf ( fp,
	      "%8.3f %8.3f %8.3f\n",
	      donor_template_points[i].pos[X],
	      donor_template_points[i].pos[Y],
	      donor_template_points[i].pos[Z] );

	fprintf ( fp_tmp, 
		  "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 50.00\n", 
		  i + 1,
		  i + 1,
		  donor_template_points[i].pos[X],
		  donor_template_points[i].pos[Y],
		  donor_template_points[i].pos[Z] );
	fprintf ( fp_pdb, 
		  "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 50.00\n", 
		  2000 + i,
		  2000 + i,
		  donor_template_points[i].pos[X],
		  donor_template_points[i].pos[Y],
		  donor_template_points[i].pos[Z] );
	fprintf ( fp_template, 
		  "D*  %7.3f %7.3f %7.3f\n", 
		  donor_template_points[i].pos[X],
		  donor_template_points[i].pos[Y],
		  donor_template_points[i].pos[Z] );

  	for ( j = 0; j < number_of_hphob_cluster_points; j++ )
    {
      if(hphob_cluster_points[j].type==HPHOB)
      {
        dist = find_dist(donor_template_points[i].pos[X], 
                       donor_template_points[i].pos[Y], 
                       donor_template_points[i].pos[Z], 
                       hphob_cluster_points[j].pos[X], 
                       hphob_cluster_points[j].pos[Y], 
                       hphob_cluster_points[j].pos[Z]);
        if(dist<threshold) 
        { 
            hphob_cluster_points[j].type=NOTHING;
            number_of_deleted_points++; 
            break; 
        }
      }
    }
  }

  fclose ( fp_tmp );
  fclose ( fp );

  /** Doneptor output **/
  printf ( "%d doneptor cluster points\n", number_of_doneptor_template_points);
  fprintf(fp_log, "%d doneptor cluster points\n", 
          number_of_doneptor_template_points);
  snprintf(filename, FILENAME_MAX, "%s/doneptor.pdb", template_dir);
  fp_tmp = slide_fopen ( filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open tmp file");
  snprintf(filename, FILENAME_MAX, "%s/doneptor.usr", template_dir);
  fp = slide_fopen ( filename, "w" );
  if (fp == NULL) err_panic ( "main", "Unable to open usr file");
  fprintf ( fp, "DOTS\n" );

  for ( i = 0; i < number_of_doneptor_template_points; i++ )
  {
     j=number_of_hbond_cluster_points;
     hbond_cluster_points[j].pos[X] = doneptor_template_points[i].pos[X];
     hbond_cluster_points[j].pos[Y] = doneptor_template_points[i].pos[Y];
     hbond_cluster_points[j].pos[Z] = doneptor_template_points[i].pos[Z];
     hbond_cluster_points[j].type = DONEPTOR; 
     number_of_hbond_cluster_points++;

    fprintf ( fp,
       "%8.3f %8.3f %8.3f\n",
        doneptor_template_points[i].pos[X],
        doneptor_template_points[i].pos[Y],
        doneptor_template_points[i].pos[Z] );

    fprintf ( fp_tmp, 
        "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 25.00\n", 
        i + 1,
        i + 1,
        doneptor_template_points[i].pos[X],
        doneptor_template_points[i].pos[Y],
        doneptor_template_points[i].pos[Z] );
    fprintf ( fp_pdb, 
	    "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00 25.00\n", 
	    3000 + i,
	    3000 + i,
	    doneptor_template_points[i].pos[X],
	    doneptor_template_points[i].pos[Y],
	    doneptor_template_points[i].pos[Z] );
    fprintf ( fp_template, 
	    "N*  %7.3f %7.3f %7.3f\n", 
	    doneptor_template_points[i].pos[X],
	    doneptor_template_points[i].pos[Y],
	    doneptor_template_points[i].pos[Z] );

  	for ( j = 0; j < number_of_hphob_cluster_points; j++ )
    {
      if(hphob_cluster_points[j].type==HPHOB)
      {
        dist = find_dist(doneptor_template_points[i].pos[X], 
                       doneptor_template_points[i].pos[Y], 
                       doneptor_template_points[i].pos[Z], 
                       hphob_cluster_points[j].pos[X], 
                       hphob_cluster_points[j].pos[Y], 
                       hphob_cluster_points[j].pos[Z]);
        if(dist<threshold) 
        { 
            hphob_cluster_points[j].type=NOTHING;
            number_of_deleted_points++; 
            break; 
        }
      }
    }
  }
  fclose ( fp_tmp );
  fclose ( fp );

  printf ( "%d deleted points\n", number_of_deleted_points );
  number_of_hphob_cluster_points-=number_of_deleted_points;
  printf ( "%d unbumped hydrophobic cluster points\n", number_of_hphob_cluster_points );
  fflush ( stdout );
  fprintf ( fp_log, "%d unbumped hydrophobic cluster points\n", number_of_hphob_cluster_points );

  snprintf(filename, FILENAME_MAX, "%s/hphob.pdb", template_dir);
  fp_tmp = slide_fopen ( filename, "w" );
  if (fp_tmp == NULL) err_panic ("main", "Unable to open tmp file");
  snprintf(filename, FILENAME_MAX, "%s/hphob.usr", template_dir);
  fp = slide_fopen ( filename, "w" );
  if (fp == NULL) err_panic ( "main", "Unable to open usr file");
  fprintf ( fp, "DOTS\n" );

  for ( i = 0; i < (number_of_hphob_cluster_points+number_of_deleted_points); i++ )
  {
    if( hphob_cluster_points[i].type==HPHOB)
    {
      fprintf ( fp, 
	      "%8.3f %8.3f %8.3f\n",
	      hphob_cluster_points[i].pos[X],
	      hphob_cluster_points[i].pos[Y],
	      hphob_cluster_points[i].pos[Z] );

      fprintf ( fp_tmp, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00100.00\n", 
		i + 1,
		i + 1,
		hphob_cluster_points[i].pos[X],
		hphob_cluster_points[i].pos[Y],
		hphob_cluster_points[i].pos[Z] );
      fprintf ( fp_pdb, 
                "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00100.00\n", 
		i + 1,
		i + 1,
		hphob_cluster_points[i].pos[X],
		hphob_cluster_points[i].pos[Y],
		hphob_cluster_points[i].pos[Z] );
      fprintf ( fp_template, 
		"H   %7.3f %7.3f %7.3f\n", /* toneroma 06FEB07 - turned HPHOB key points off */
		hphob_cluster_points[i].pos[X],
		hphob_cluster_points[i].pos[Y],
		hphob_cluster_points[i].pos[Z] );
    }
  }
  fclose ( fp );
  fclose ( fp_tmp );
  fclose ( fp_pdb );
  fclose ( fp_template );
  fclose ( fp_log );
  exit ( 0 );
}

void print_template_header(FILE* template, FILE* log, const char* target_str, 
                           float grid_spacing, float hphob_threshold, 
                           const int argc, const char** argv)
{
  char time_buf[256];

  /* This is a mess no matter the method I try -- this is why people tend to
   * use C++ strings, python, etc for I/O rather than C 
   *
   * A waste of time to try to do it elegantly.  To keep things eaiser for 
   * those unacustomed to pointers, just copy the print statements.
   */
  printf("# SLIDE template; version: %s\n", VERSION);
  fprintf(template, "# SLIDE template; version: %s\n", VERSION);
  fprintf(log, "# SLIDE template; version: %s\n", VERSION);
  printf("# Command-Line Args:\n");
  fprintf(template, "# Command-Line Args:\n");
  fprintf(log, "# Command-Line Args:\n");
  printf("#   Target :                             %s\n", target_str);
  fprintf(template, "#   Target :                             %s\n", 
          target_str);
  fprintf(log, "#   Target :                             %s\n", target_str);
  printf("#   Template type :                      unbiased template\n");
  fprintf(template, "#   Template type :                      unbiased template"
          "\n");
  fprintf(log, "#   Template type :                      unbiased template\n");
  printf("#   Hydrogen bond points density :       %s\n", argv[3]);
  fprintf(template, "#   Hydrogen bond points density :       %s\n", 
          argv[3]);
  fprintf(log, "#   Hydrogen bond points density :       %s\n", argv[3]);
  printf("#   Grid spacing :                       %.2f (A)\n", grid_spacing);
  fprintf(template, "#   Grid spacing :                       %.2f (A)\n", 
          grid_spacing);
  fprintf(log, "#   Grid spacing :                       %.2f (A)\n", 
          grid_spacing);
  printf("#   Hydrophobic clustering threshold :   %.2f (A)\n#\n", 
         hphob_threshold);
  fprintf(template, "#   Hydrophobic clustering threshold :   %.2f (A)\n#\n", 
          hphob_threshold);
  fprintf(log, "#   Hydrophobic clustering threshold :   %.2f (A)\n#\n", 
          hphob_threshold);
  if(argc == 8){
    printf("#   Reference ligand file :              %s\n", argv[6]);
    fprintf(template, "#   Reference ligand file :              %s\n", argv[6]);
    fprintf(log, "#   Reference ligand file :              %s\n", argv[6]);
    printf("#   Distance from reference ligand :     %s (A)\n", argv[7]);
    fprintf(template, "#   Distance from reference ligand :     %s (A)\n", 
            argv[7]);
    fprintf(log, "#   Distance from reference ligand :     %s (A)\n", argv[7]);
  }
  if(slide_get_local_time(time_buf, 256)){
    printf("#\n# Local start time: %s\n", time_buf);
    fprintf(template, "#\n# Local start time: %s\n", time_buf);
    fprintf(log, "#\n# Local start time: %s\n", time_buf);
  }

}

#if 0
/* There has to be a way to get the number of strings without explicitly passing
 * the number -- e.g. fprintf()
 * Oh -- they use the format string of course -- well there is nothing
 * to be gained from recoding the parsing of the format string -- probalby
 * better to use lex/yacc anyhow.
 */
void print_template_and_log(FILE *template, FILE *log, const char *format, 
                            int num_strs, ...)
{
  va_list al;
  char** strings;
  int i;

  va_start(al, num_strs);
  strings = (char**) mymalloc(num_strs * sizeof(char*) );
  for(i = 0; i < num_strs);
  
   
  fprintf(fp_template, "%s", header);
  printf("%s", header);
  fprintf(fp_log, "%s", header);
}
#endif
