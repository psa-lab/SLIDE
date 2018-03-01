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
3.0.1 of the SLIDE software for screening molecular databases for
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

#define _MAIN_
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include <sys/param.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <netdb.h>
#include <errno.h>
#include <basics.h>
#include "defs.h"
#include "types.h"
#include "malloc.h"
#include "read_pdb.h"
#include "err_handle.h"
#include "initialize.h"
#include "read_template.h"
#include "read_pts_file.h"
#include "hashing.h"
#include "match_triangles.h"
#include "read_flex_defn.h"
#include "neighbors.h"
#include "read_waters.h"
#include "intra_hbonds.h"
#include "read_parameter_file.h"
#include "write_log_file.h"
#include "screen_single_compounds.h"
#include "demo.h"
#include "slide_itimer.h"
#include "mymalloc.h"
#include <distance_matrices.h>

/*#define TRACE*/

#ifdef DEMO_VERSION
#include "check_key.h"
#endif


int  main ( int  argc, char  **argv )
{
  struct stat file_stat_buf;
  FILE            *fp;
  global_data_pt  global;
  char            file[FILENAME_MAX + 1],
    pts_file[FILENAME_MAX + 1],
    hname[MAXHOSTNAMELEN],
    line[FILENAME_MAX],
    time_buf[256];
  int             number_of_template_points;
  int i;
  int dbase_entries;
  double elapsed_real_time;
  double elapsed_virt_time;
  double elapsed_prof_time;
  char  err_msg[FILENAME_MAX];
  char  database_num_inputs_string[FILENAME_MAX];
  int errsv;  /* Need to save errno since printfs can overwrite errno */

  if(gethostname(hname, MAXHOSTNAMELEN)<0)
    err_warning("main", "Cannot determine Hostname");

  if ( argc < 4 ) {
    sprintf ( line, "%s <protein> <template> <database> [<mol2-file>..]\n",
	      argv[0] );
    err_usage ( line );
  }
  printf ( "\n" );
  printf ( "Command-Line Args: ");
  for ( i = 0; i < argc; i++ )
    printf ( "%s ", argv[i] );
  printf ( "\n\nSLIDE Version: %s", VERSION);
  if(slide_itimer_start(10)) return -1;
  if(slide_get_local_time(time_buf, 256)){
    printf("\n\nstarted on %s at %s\n", hname, time_buf );
    fprintf(stderr, "\n\nrun started on %s at %s\n\n", hname, time_buf );
    sprintf(err_msg, "\n\n*****************************************\n\n"
            "run started on %s at %s\n\n", hname, time_buf );
  }
  global = initialize_global_data_structure ( );
  
  strncpy ( global->slide_dir, getenv ( "SLIDE_DIR" ), FILENAME_MAX );
  if ( *global->slide_dir == '\0' )
    err_panic( "main", "SLIDE_DIR environment variable not set" );
  strncpy ( global->data_root, getenv ( "SLIDE_DATA_DIR" ), FILENAME_MAX );
  if ( *global->data_root == '\0' )
    err_panic( "main", "SLIDE_DATA_DIR environment variable not set" );
  printf ( "SLIDE_DIR = %s\n", global->slide_dir );
  printf ( "SLIDE_DATA_DIR = %s\n\n", global->data_root );

#ifdef DEMO_VERSION

  check_key(global->slide_dir);

#endif

  strcpy ( global->protein, argv[1] );
  strcpy ( global->template, argv[2] );
  strcpy ( global->database, argv[3] );

  /* Initialize the static log file name stored in err_handle.c */
  sprintf(file, "%s/%s/%s/log/%s.err", global->data_root, global->protein,
          global->template, global->database );
  set_err_filename(file);
  err_print(err_msg);

  /* first read the global parameter file */
  printf( "\nReading the default parameters set:");
  snprintf(file, sizeof(file), "%s/params/slide.parameters", global->slide_dir);
  if(read_parameter_file ( global, file ) == FAILURE)
    err_panic2("main", "global parameter file not found or has a bad "
              "parameter name");

  /* now check, if there is a local parameter file, all settings in that
     file overrule the global settings */
  sprintf(file, "%s/%s/%s/in/slide.parameters", global->data_root,
	  global->protein, global->template );
  if (read_parameter_file(global, file) == SUCCESS )
    printf ("These new parameters will overwrite the default parameters for this run.\n");
  
  /* read the compile-time parameters file and print out the user adjustable 
   * parameters */
  printf("\nCompile-time parameters as defined in: %s/src/slide/inc/params.h\n",
         global->slide_dir );
  printf("> TRIANGLE_MIN_PERIMETER (A):          %5.2f\n", 
         TRIANGLE_MIN_PERIMETER);
  printf("> TRIANGLE_MAX_PERIMETER (A):          %5.2f\n", 
         TRIANGLE_MAX_PERIMETER);
  printf("> TRIANGLE_MIN_SHORTEST_SIDE (A):      %5.2f\n", 
         TRIANGLE_MIN_SHORTEST_SIDE);
  printf("> TRIANGLE_MAX_SHORTEST_SIDE (A):      %5.2f\n", 
         TRIANGLE_MAX_SHORTEST_SIDE);
  printf("> TRIANGLE_MAX_LONGEST_SIDE (A):       %5.2f\n", 
         TRIANGLE_MAX_LONGEST_SIDE);
  printf("> SMALL_TRIANGLE_MIN_LONGEST_SIDE (A): %5.2f\n", 
         SMALL_TRIANGLE_MIN_LONGEST_SIDE);
  printf("> LARGE_TRIANGLE_MIN_LONGEST_SIDE (A): %5.2f\n", 
         LARGE_TRIANGLE_MIN_LONGEST_SIDE);
  printf("> SMALL_LIGAND_DIAMETER (A):           %5.2f\n", 
         SMALL_LIGAND_DIAMETER);

  /* now see, if any environment variables are set for the parameters */
  check_environment_parameter_variables ( global );
  sprintf ( file, "%s/params/flex.defn", global->slide_dir );
  printf ( "\nreading flex-bond definition file %s\n", file );
  global->number_of_flex_bond_rules =
    read_flex_defn ( file, global->flex_bond_rules);

  sprintf ( file, "%s/%s/%s/in/%s.rad", global->data_root, global->protein,
	    global->template, global->protein );
  printf ( "reading target-file %s\n", file );

  read_pdb(file, global->target_atoms, global->target_residues, ALSO_HETERO,
           &global->number_of_target_atoms, &global->number_of_target_residues);

  /* make backup of target atom positions, since the positions in
     'global->target_atoms' are modified when target side chains
     are rotated during bump resolvement */
  global->target_atoms->neighbors = 
    (atom_pt *) mymalloc(global->number_of_target_atoms * MAX_NEIGHBOR_ATOMS * 
                         sizeof(atom_pt));
  global->target_atoms->neighbor_dist = 
    (float *) mymalloc(global->number_of_target_atoms * MAX_NEIGHBOR_ATOMS *
                       sizeof(float));
  for(i = 0; i < global->number_of_target_atoms; i++){
    global->orig_target_atom_act[i] = global->target_atoms[i].act;
    global->target_atoms[i].neighbors = 
      &(global->target_atoms->neighbors[i*MAX_NEIGHBOR_ATOMS]);
    global->target_atoms[i].neighbor_dist = 
      &(global->target_atoms->neighbor_dist[i*MAX_NEIGHBOR_ATOMS]);
    global->target_atoms[i].num_nbrs = 0;
    global->target_atoms[i].num_added_nbrs = 0;
  }
  memcpy(global->orig_target_atom_positions, global->target_atom_positions,
         3*global->number_of_target_atoms *
         sizeof(*global->target_atom_positions));
  memset(global->target_rotations, 0, global->number_of_target_residues * 
         sizeof(*global->target_rotations));
  memset(global->target_intra_overlap, 0, global->number_of_target_residues * 
         sizeof(*global->target_intra_overlap));

  /* this is the array for the lookup table of inter-atomic distances,
     before doing the very first bump-check after transforming a ligand,
     this array is filled with the distances between all pairs of ligand
     and target atoms, so that we will avoid most of the calls of
     'dist_fun()' during the modeling of the induced complementarity - Volker*/
  /* Sameer 20 Nov 04 - Integrating scoring function (45) code from Litian.
   * Below two tables are allocated before the call to identify_atom_neighbors
   * in main_score.c by Litian He.
   */

  global->target_ligand_distances = 
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS * 
                        global->number_of_target_atoms * sizeof (float) );
  global->target_ligand_sq_dists = 
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS * 
                        global->number_of_target_atoms * sizeof (float) );
  distance_array(&global->target_dists_array, global->target_atoms, 
                 global->number_of_target_atoms, 4.0, 5.0);
  init_target_nbr_arrays(global->target_atoms, global->number_of_target_atoms);

  /* Allocate & initialize memory for flag arrays used in the scoring function 
   */
  global->target_flag =
    (short *) mymalloc ( global->number_of_target_atoms * sizeof(short) );

#ifdef NON_METALBONDED_REPULSIVE
  /* Build a list of metal atom indices */
  global->number_of_metals = 0;
  global->metal_atom_indices = (int*) mymalloc (MAX_TARGET_METAL_ATOMS * sizeof(int) );

  for(i = 0; i < global->number_of_target_atoms; i++)
    if( global->target_atoms[i].act == METAL_1 || global->target_atoms[i].act == METAL_2)
      global->metal_atom_indices[global->number_of_metals++] = i;

  if(MAX_TARGET_METAL_ATOMS <= global->number_of_metals)
    err_panic2 ( "main", "Number of metals exceed MAX_TARGET_METAL_ATOMS");
#endif

  /****************SA memory allocation modifications end********************/
  sprintf(file, "%s/databases/%s/number_of_mol2_entries", global->data_root, 
          global->database);
  if(stat(file, &file_stat_buf)){
    printf("\nUnable to open the database file %s\n%s\n", file, 
           strerror(errno));
    global->database_num_exist = FALSE;
  }else if(file_stat_buf.st_size == 0){
    printf("Database file %s is empty\n", file);
    global->database_num_exist = FALSE;
  }
  fp = slide_fopen(file, "r");
  if(fp == NULL) printf("Unable to open database file: %s\n", file);
  else{ 
    dbase_entries = 0;
    while( fgets(line, FILENAME_MAX, fp))
      if ( line[0] != ' ' && line[0] != '\0' && line[0] != '\n' ) {
	sscanf ( line, "%s",database_num_inputs_string );
	dbase_entries++;
	global->database_num_conformers = atoi(database_num_inputs_string);
	global->database_num_exist = TRUE;
      }
    if(dbase_entries == 0) global->database_num_exist = FALSE;
  }
  
  
  sprintf(file, "%s/%s/%s/in/template",
	  global->data_root, global->protein, global->template );
  printf("reading template-file %s\n",   file );
  number_of_template_points = read_template ( file, global );
  global->number_of_waters = 0;
  printf("Water handling disabled in current version of SLIDE => #waters = 0"
         "\n\n");
#if 0
  /* toneroma 06FEB07 - handling of waters removed from SLIDE */
  sprintf(file, "%s/%s/%s/in/waters.pdb", global->data_root, global->protein,
	  global->template );
  printf ( "reading water-file %s\n", file );
      if ( read_waters ( global, file ) == SUCCESS ) {
      if(global->number_of_waters > 0 ) {
      global->target_water_distances = (float *) 
      mymalloc(global->number_of_target_atoms * global->number_of_waters * 
      sizeof (float ) );
      global->ligand_water_distances = (float *) 
      mymalloc(MAX_NUMBER_OF_MOL2_ATOMS * global->number_of_waters * 
      sizeof (float ) );
      global->intra_water_distances = (float *)
      mymalloc(global->number_of_waters * global->number_of_waters * 
      sizeof (float ) );
      }
      printf ( "water-file being used\n\n", file );
      } 
      else
      fprintf(stderr, "no water-file or not able to open water-file %s => "
      "#waters = 0\n\n", file );
  for ( i = 0; i < global->number_of_waters; i++ )
    for ( j = 0; j < 3; j++ )
      global->orig_water_positions[i][j] = global->waters[i].pos[j];
  if ( global->number_of_waters > 0 )
    identify_water_neighbors ( global->waters, global->target_atoms,
			       global->number_of_waters,
			       global->number_of_target_atoms,
			       global);
#endif

  global->number_of_orig_intra_target_hbonds
    = intra_target_hbonds ( global ) + target_to_water_hbonds ( global );

  printf ( "target: %d residues, %d atoms, %d waters\n",
	   global->number_of_target_residues,
	   global->number_of_target_atoms, global->number_of_waters );
  printf ( "template: %d points, %d key points, %d secondary key points\n",
	   number_of_template_points, global->template_key_points, global->template_key2_points );
  if(global->database_num_exist == TRUE)
    printf("database: %d total input molecules or conformers\n",
	   global->database_num_conformers );
  printf("internal H-bonds: %d intra-target, %d to water\n\n",
	 intra_target_hbonds ( global ), target_to_water_hbonds ( global ) );
  if(create_hash_table(global) == FALSE){
    return -1;
  }

  /***********************CODEMM2*****************************************/
  if(argc == 4){
    /* we want to screen a database */
    sprintf(file, "%s/databases/%s/%s.db", global->data_root, global->database,
            global->database );
    if(stat(file, &file_stat_buf)){
      errsv = errno;
      sprintf(err_msg, "\nFATAL ERROR: Unable to open the file %s\n%s\n\n", 
              file, strerror(errsv));
      fprintf(stderr, err_msg);
      err_print(err_msg);
      return -1;
    }else if(file_stat_buf.st_size == 0){
      sprintf(err_msg, "ERROR: Database file %s is empty\n", file);
      fprintf(stderr, err_msg);
      err_print(err_msg);
    }
    fp = slide_fopen(file, "r");
    if(fp == NULL){
      sprintf( err_msg, "database: %s\n", file );
      fprintf( stderr, err_msg);
      err_print(err_msg);
      err_panic2( "main", "unable to open database");
    }
    
    printf ( "\nreading database-definition file %s\n\n", file );


    printf("***************************************************************"
           "*****************\n"
           "1)  [Ligand Name]_[Conformer]_[Orientation]\n"
           "2)  Orientation Score\n"
           "3)  Affinity Score / # Heavy Ligand Atoms (Ligand Efficiency)\n"
           "4)  Affinity Score\n"
           "5)  Buried Protein Hydrophobic Term\n"
           "6)  Hydrophobic Complementarity Term\n"
           "7)  Polar Component Term\n"
           "8)  Number of Protein-Ligand Hydrophobic Contacts\n"
           "9)  Number of Protein-Ligand H-bonds\n"
           "10) Number of Protein-Ligand Salt-bridges\n"
           "11) Number of Metal-Ligand Bonds\n"
           "12) Number of Interfacial Unsatisfied Polar Atoms\n"
           "13) Number of Interfacial Unsatisfied Charged Atoms\n"
           "14) Buried Carbons\n"
           "15) Remaining van der Waals Collisions\n"
           "16) Remaining van der Waals Overlap\n"
           "****************************"
           "****************************************************\n\n  1                  2      3       4      5       6   7   8   9  10  11  12    13  14\n");

    dbase_entries = 0;
    while( fgets(line, FILENAME_MAX, fp)){
      if(line[0] == ' ' || line[0] == '\0' || line[0] == '\n') continue;

#ifdef TRACE
      printf("%s\n",line);
#endif
      sscanf(line, "%s %s", global->compound_dir, pts_file);
      strcpy(global->compound_name, pts_file);
      sprintf(file, "%s/databases/%s/%s", global->data_root, global->database,
              pts_file );
      if(read_pts_file(global, file) == SUCCESS){
        match_triangles(global);
        dbase_entries++;
      }else if(global->restart_molecule_check == FALSE){
        sprintf(err_msg, "\nWarning!\npoints file %s could not be read\n",
                file);
        fprintf(stderr, err_msg);
        err_print(err_msg);
      }
    }
    if(dbase_entries == 0){
      sprintf(err_msg, "\nWarning!\ndatabase file %s/databases/%s/%s.db\n"
              "has no valid entries.\n\n", global->data_root, global->database,
              global->database);
      fprintf(stderr, err_msg);
      err_print(err_msg);
    }
    if(global->database_num_exist == TRUE){
      /*      fprintf(stderr, "Screened %d molecules, (total screened conformers: "
	      "%d => %.3f %%)\r", global->total_num_molecules_noconf,
	      global->number_of_screened_compounds,
	      (100*((float)global->number_of_screened_compounds) /
	      (float)global->database_num_conformers));*/
      fprintf(stderr, "Screened:%8d(%6.3f%%)||Top:%-7.3f(AVG:%-7.3f,STD:%6.3f) %s\r", global->total_num_molecules_noconf, 
              (100*((float)global->number_of_screened_compounds-1) /
	       (float)global->database_num_conformers),
	      global->best_affiscore_so_far,
	      global->affiscore_mean,
	      global->affiscore_stdd,
	      global->best_affi_name);
    }else{
      /*      fprintf(stderr, "Screened %d molecules, (total screened conformers: %d)"
              "\r", global->total_num_molecules_noconf, 
              global->number_of_screened_compounds-1);*/
      fprintf(stderr, "Screened:%8d||Top:%-7.3f(AVG:%-7.3f,STD:%6.3f) %s\r", global->total_num_molecules_noconf, 
	      global->best_affiscore_so_far,
	      global->affiscore_mean,
	      global->affiscore_stdd,
	      global->best_affi_name);
    }
    
    
    fclose(fp);
    sprintf(file, "%s/%s/%s/log/%s.log", global->data_root, global->protein,
            global->template, global->database );
    printf("\nwriting log-file %s\n", file);
    write_log_file(global, file);
  }
  /****************************************************************/
  else {
    /* there is a list of mol2-files specified that we are going to
       screen */
    if ( strncmp ( argv[4] + strlen ( argv[4] ) - 5, ".mol2", 5 ) == 0 )
      screen_single_compounds ( global, argv, argc, (FILE *) NULL );
    else {
      fp = slide_fopen(argv[4], "r" );
      if ( fp == NULL ) {
	sprintf ( line, "unable to open file %s", argv[4] );
	err_panic2 ( "main", line);
      }
      fgets ( line, sizeof ( line ), fp );
      if ( strchr ( line, ' ' ) != NULL ){
	/* this line contains two strings, so we assume that it is a
	   database definition file */
	printf ( "\ngetting compounds to screen from " );
	printf ( "local database definition file %s\n\n", argv[4] );
	/* now the screening begins, so redirect stderr to the
	   error-log file */
	do
	  if ( line[0] != ' ' && line[0] != '\0' && line[0] != '\n' ) {
	    sscanf ( line, "%s %s", global->compound_dir, pts_file );
	    /** reading file from the databases directory instead of global compound
		dir which might require root permission to modify
		-- RSK (apr 4, 2002) **/
	    sprintf ( file, "%s/databases/%s/%s",
		      global->data_root, global->database, pts_file );
	    printf ( "processing pts-file %s\n", file );
	    if ( read_pts_file ( global, file ) == SUCCESS )
	      match_triangles ( global );
	  }
	while ( fgets ( line, sizeof ( line ), fp ) );
	sprintf ( file, "%s/%s/%s/log/%s.log", global->data_root,
		  global->protein, global->template, global->database );
	
	printf ( "\nwriting log-file %s\n", file );
	write_log_file ( global, file );
      } else
	/* since there is only one entry in the first line of the file,
	   we assume that it contains a list of mol2 files for screening */
	screen_single_compounds ( global, argv, argc, fp );

      fclose ( fp );
    }
  }
#ifndef OUTPUT_ALL_MATCHES
#ifndef VERBOSE_OUTPUT  
  printf ("\nFor each ligand that passed the testing criteria the top scoring orientation was saved.\n");
#endif
#endif
  
  printf ("\nPlease check the error file %s/%s/%s/log/%s.err to identify any problems with the current run.\n", global->data_root, global->protein, global->template, global->database );

  if(slide_get_local_time(time_buf, 256)){
    printf ( "\nrun finished at %s\n", time_buf);
    sprintf(err_msg, "\nrun finished at %s\n\n********************************"
            "*********\n\n", time_buf);
    fprintf(stderr, err_msg);
    err_print(err_msg);
  }
  if(slide_itimer_get(&elapsed_real_time, &elapsed_virt_time, 
                      &elapsed_prof_time))
    fprintf(stderr, "Unable to get the elapsed time from the profile timer\n");
  else{
    printf("Screening time for %d %s (%d %s)\n\tWall Clock Time:  \t\t%8.2f "
           "sec\n", global->total_num_molecules_noconf, 
           (global->total_num_molecules_noconf > 1) ? "molecules" : 
           "molecule", global->number_of_screened_compounds, 
           (global->number_of_screened_compounds > 1) ? "conformers" : 
           "conformer", elapsed_real_time);
    printf("\tSystem Time:      \t\t%8.2f sec\n", elapsed_prof_time);
    printf("\tUser Space Time:  \t\t%8.2f sec\n", elapsed_virt_time);
    printf("\tKernel Space Time:\t\t%8.2f sec\n", elapsed_prof_time - 
           elapsed_virt_time);
  }

  if(global->target_atoms->neighbors) free(global->target_atoms->neighbors);
  if(global->target_atoms->neighbor_dist) 
    free(global->target_atoms->neighbor_dist);
  global->target_atoms->neighbors = 0;
  global->target_atoms->neighbor_dist = 0;

  return 0;
}
