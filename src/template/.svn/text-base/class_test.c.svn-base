#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "types.h"
#include "basics.h"
#include "read_pdb_file.h"
#include "read_borders_file.h"
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

int main ( int argc, char **argv ) {
  
  FILE *fp_pdb;
  pdb_atom_pt atoms;
  char filename[256],
       classarray[6];
  int number_of_atoms;
  int i,j;


  if ( argc != 2 ) {
    printf ("Usage: %s <pdb file>\n",argv[0]);
    exit (1);
  }  

  sprintf( filename, "%s", argv[1]);
  atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
  number_of_atoms = read_pdb_file ( filename,
				    atoms );
  /*  for (i=0; i<number_of_atoms; i++) {
    printf ("%5s  |%3s|  %3s\n", atoms[i].residue_number, 
                                    atoms[i].residue_name,
                                    atoms[i].atom_name);
  }*/
  assign_hphil ( atoms, number_of_atoms );    
  assign_atom_class ( atoms, number_of_atoms );  /* PCS-06Jun00 */

  /* initialize an array to translate class ints to output chars */
  classarray[CLASS_A] = 'A';
  classarray[CLASS_D] = 'D';
  classarray[CLASS_C] = 'C';
  classarray[CLASS_O] = 'O';
  classarray[CLASS_N] = 'N';
  classarray[CLASS_H] = 'H';

  for (i=0; i<number_of_atoms; i++) {
    printf ("%5s  %3s  %-3s  %1c\n", atoms[i].residue_number, 
                                    atoms[i].residue_name,
                                    atoms[i].atom_name,
                                    classarray[atoms[i].class]);
  }
  
  return OK;
}
