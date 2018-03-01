#include <string.h>
#include "defs.h"
#include "types.h"


/* P. Sanschagrin -- 5-Jun-00 */
/* assign_atom_class.c -- assign the class of the atom, where class is 
 * defined as acceptor, donor, doneptor, hydrophobic (S or C bonded only 
 * to C or H), non-hydrophobic C, or other. This is based solely on atom
 * and residue type and does not check any hydrogen bond angles/distances
 */

int assign_atom_class ( pdb_atom_pt atoms,
			int         number_of_atoms ) {

  int restype, atomtype;
  int i;

  for (i = 0 ; i < number_of_atoms ; i++) {
    switch (atoms[i].atom_name[0]) {
    case 'N': 
      atoms[i].class = CLASS_D; /* Atom is a DONOR */
      break;
    case 'O':
      switch (atoms[i].atom_name[1]) {
      case 'H':
      case 'G':
	atoms[i].class = CLASS_N; /* Atom is a DONEPTOR */
	break;
      default:
	atoms[i].class = CLASS_A; /* Atom is an ACCEPTOR */
	break;
      }
      break;
    case 'S': 
      atoms[i].class = CLASS_H; /* All S's are HYDROPHOBIC */
      break;      
    case 'C': /* Carbon -- must distinguisg b/n hphobic C and non-h C */
      atoms[i].class = CLASS_C; /* Atom is non-H CARBON unless reset below */
      if ( strncmp(atoms[i].residue_name,"CA",2) ==  0){
	atoms[i].class = CLASS_O; /* This CA is Ca2+ -- atom class is OTHER */
	break;
      }
      else if (atoms[i].residue_type == ALA) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == ARG) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == ASN) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == ASP) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == GLN) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == GLU) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == HIS) {
	if (atoms[i].type == CB)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == ILE) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG1) ||
	    (atoms[i].type == CG2) ||
	    (atoms[i].type == CD1))
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == LEU) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG) ||
	    (atoms[i].type == CD1) ||
	    (atoms[i].type == CD2)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == LYS) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG) ||
	    (atoms[i].type == CD)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == MET) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG) ||
	    (atoms[i].type == CE))
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == PHE) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG) ||
	    (atoms[i].type == CD1) ||
	    (atoms[i].type == CD2) ||
	    (atoms[i].type == CE1) ||
	    (atoms[i].type == CE2) ||
	    (atoms[i].type == CZ)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == PRO) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == THR) {
	if (atoms[i].type == CG2)
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == TRP) {
	if ((atoms[i].type == CB)  ||
	    (atoms[i].type == CG)  ||
	    (atoms[i].type == CD2) ||
	    (atoms[i].type == CZ2) ||
	    (atoms[i].type == CH2) ||
	    (atoms[i].type == CE3) ||
	    (atoms[i].type == CZ3)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == TYR) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG) ||
	    (atoms[i].type == CD1) ||
	    (atoms[i].type == CE1) ||
	    (atoms[i].type == CD2) ||
	    (atoms[i].type == CE2))
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      else if (atoms[i].residue_type == VAL) {
	if ((atoms[i].type == CB) ||
	    (atoms[i].type == CG1) ||
	    (atoms[i].type == CG2)) 
	  atoms[i].class = CLASS_H; /* Atom is HYDROPHOBIC */
      }
      break;
     default:
      atoms[i].class = CLASS_O; /* Atom isn't one of above, therefore assigned
			   * as class OTHER */
      break;
    }
  }
  return OK;
}
