/*
 *
 *    assign_hydrogens.c   Volker Schnecke     Mon Feb 16 21:17:11 EST 1998
 *
 *    function: assign_hydrogens()
 *
 *    assigns hydrogens to those nitrogens and oxygens where they seem to 
 *    be missing in the structure (adds only atoms but no positions)
 *
 */
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "basics.h"

/*
 *  This routine checks all nitrogens and oxygens in 'molecule' by counting 
 *  the number of bonds (considering double bonds as two bonds, triple as 
 *  three) and adding hydrogens, if less than three (two for oxigen) bonds 
 *  are found.
 */
void  assign_hydrogens ( molecule_pt  molecule )
{
  atom_pt  atom;
  bond_pt  bond;
  int      number_of_atoms,
           bond_type,
           bond_counter,
           total_bonds,
           index;
  int      i, j;

  number_of_atoms = molecule->number_of_atoms;
  for ( i = 0; i < number_of_atoms; i++ )
    if ( molecule->atoms[i].type == N || molecule->atoms[i].type == O )
      /* only check nitrogens and oxygens */
      {
	if ( molecule->atoms[i].type == N )
	  /* nitrogen should have three bonds, otherwise it is assumed
	     that hydrogens are not included in the structure */
	  total_bonds = 3;
	else
	  /* oxygens should have two bonds */
	  total_bonds = 2;
	if ( molecule->number_of_neighbors[i] != total_bonds )
	  {
#ifdef TRACE
	    printf ( "atom %d: %s needs %d bonds, has %d neighbors, bond_count = ", 
		     i + 1, molecule->atoms[i].type_str,
		     total_bonds, molecule->number_of_neighbors[i] );
#endif
	    bond_counter = 0;
	    for ( j = 0; j < molecule->number_of_neighbors[i]; j++ )
	      {
		bond_type = 
		  molecule->bonds[molecule->bond_to_neighbor[i][j]].type;
#ifdef TRACE
		printf ( "(%d) ", bond_type );
#endif
		if ( bond_type < 4 )
		  /* count double bonds twice, triple threefold */
		  bond_counter += bond_type;
		else
		  /* aromatic bonds count double, although it would be
		     more correct to count them as 1.5 */
		  if ( bond_type == AROMATIC )
		    bond_counter += 2;
		  else
		    /* all other bonds */
		    bond_counter++;
	      }
#ifdef TRACE
	    printf ( "%d\n", bond_counter );
#endif
	    for ( ; bond_counter < total_bonds; bond_counter++ )
	      /* if bond_counter is less than total_bonds add one hydrogen 
		 for each missing bond */
	      {
		if ( molecule->number_of_atoms > MAX_NUMBER_OF_MOL2_ATOMS )
		  err_panic ( "assign_hydrogens",
			      "more than MAX_NUMBER_OF_MOL2_ATOMS atoms" );
		index = molecule->number_of_atoms;
		/* add new atom to neighbor list */
		molecule->neighbors[i][molecule->number_of_neighbors[i]] 
		  = index;
		/* add the new bond to the list of bonds for atom i */
		molecule->bond_to_neighbor[i]
		  [molecule->number_of_neighbors[i]]
		  = molecule->number_of_bonds;
		molecule->number_of_neighbors[i]++;
		molecule->neighbors[index][0] = i;
		molecule->number_of_neighbors[index] = 1;
		/* add new atom to atom list */
		atom = &molecule->atoms[index];
		atom->type = H;
		atom->hyd = NOTHING;
#ifdef TRACE
		printf ( "added Hydrogen to %d %s\n", 
			 i + 1, 
			 molecule->atoms[i].type_str );
#endif
		/* "*_ADD" is the only marker which signals that this is
		   an artificially added atom, couldn't change the type since
		   for all rule checks this has to appear as an ordinary
		   Hydorgen */
		sprintf ( atom->type_str, "H_ADD" );
		sprintf ( atom->name, "H%d", i + 1 );
		atom->pos[X] = atom->pos[Y] = atom->pos[Z] = 0.0;
		atom->number = molecule->number_of_atoms + 1;
		molecule->number_of_atoms++;
		/* add the bond to the Hydrogen */
		bond = &molecule->bonds[molecule->number_of_bonds];
		bond->number = molecule->number_of_bonds + 1;
		bond->atom1 = i;
		bond->atom2 = index;
		bond->type = SINGLE;
		sprintf ( bond->type_str, "B_ADD" );
		molecule->number_of_bonds++;
	      }
	  }
      }
}
