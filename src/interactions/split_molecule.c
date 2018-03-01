#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"

int  mark_residues ( molecule_pt  molecule,
		     int          *label,
		     int          atom,
		     int          back,
		     int          res_counter,
		     int          number_of_marked_atoms )
{
  int  *neighbors;
  int  i;

  neighbors = molecule->neighbors[atom];
  if ( label[atom] == UNVISITED )
    {
      /* this is the first time we visit this node */
      label[atom] = VISITED;
      molecule->atoms[atom].hyd = res_counter;
      number_of_marked_atoms++;
      for ( i = 0; i < molecule->number_of_neighbors[atom]; i++ )
	/* visit all neighbors, but don't take direct way back */
	if ( neighbors[i] != back )
	  {
	    number_of_marked_atoms =
	      mark_residues ( molecule,
			      label,
			      neighbors[i],
			      atom,
			      res_counter,
			      number_of_marked_atoms );
	  }
    }
  return number_of_marked_atoms;
}
  
void  split_molecule ( molecule_pt  ligand,
		       char         *ligand_file_name )
{
  FILE    *fp;
  atom_pt atoms;
  bond_pt bonds;
  int     labels[MAX_NUMBER_OF_MOL2_ATOMS];
  char    filename[256];
  int     number_of_marked_atoms,
          number_of_residues,
          res,
          atom_count,
          bond_count;
  int     i;

  atoms = ligand->atoms;
  bonds = ligand->bonds;
  /* cut off the ".mol2" from the original file name */
  ligand_file_name[strlen(ligand_file_name)-5] = '\0';
  for ( i = 0; i < ligand->number_of_atoms; i++ )
    labels[i] = UNVISITED;
  number_of_marked_atoms = mark_residues ( ligand, 
					   labels, 
					   0, 
					   -1, 
					   0, 
					   0 );
  number_of_residues = 1;
  while ( number_of_marked_atoms < ligand->number_of_atoms )
    {
      i = 0; 
      while ( labels[i] == VISITED )
	/* find the first unvisited atom */
	i++;
      number_of_marked_atoms = mark_residues ( ligand, 
					       labels, 
					       i, 
					       -1, 
					       number_of_residues, 
					       number_of_marked_atoms );
      number_of_residues++;
    }
  /* renumber all atoms in the different residues */
  for ( res = 0; res < number_of_residues; res++ )
    {
      atom_count = 1;
      for ( i = 0; i < ligand->number_of_atoms; i++ )
	if ( atoms[i].hyd == res )
	  atoms[i].number = atom_count++;
      bond_count = 1;
      for ( i = 0; i < ligand->number_of_bonds; i++ )
	if ( atoms[bonds[i].atom1].hyd == res )
	  bonds[i].number = bond_count++;
      /* open a file to output this residue */
      sprintf ( filename, 
		"%s_%d.mol2", 
		ligand_file_name,
		res );
      fp = fopen ( filename, "w" );
      if ( fp == NULL )
	err_warning ( "split_molecule", 
		      "file open failed" );
      fprintf ( fp, "@<TRIPOS>MOLECULE\n" );
      fprintf ( fp, "%s\n", ligand->name );
      fprintf ( fp,
		"%5d %5d     1 \n",
		atom_count - 1,
		bond_count - 1 );
      fprintf ( fp, "SMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" );
      for ( i = 0; 
	    i < ligand->number_of_atoms;
	    i++ )
	if ( atoms[i].hyd == res )
	  fprintf ( fp,
		    "  %5d %-8s %9.4f %9.4f %9.4f %-9s %3d %-15s  %7.4f\n",
		    atoms[i].number,
		    atoms[i].name,
		    atoms[i].pos[X],
		    atoms[i].pos[Y],
		    atoms[i].pos[Z],
		    atoms[i].type_str,
		    atoms[i].substrid,
		    atoms[i].substrname,
		    atoms[i].charge );
      fprintf ( fp, "@<TRIPOS>BOND\n" );
      for ( i = 0; 
	    i < ligand->number_of_bonds;
	    i++ )
	if ( atoms[bonds[i].atom1].hyd == res )
	  fprintf ( fp,
		    "%5d %5d %5d %s\n",
		    bonds[i].number,
		    atoms[bonds[i].atom1].number,
		    atoms[bonds[i].atom2].number,
		    bonds[i].type_str );
      fprintf ( fp, "@<TRIPOS>SUBSTRUCTURE\n      1 RES1       1\n" );
      fclose ( fp );
    }
}


