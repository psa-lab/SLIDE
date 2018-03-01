#include <stdlib.h>
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"
#include "mymalloc.h"

int  check_atoms ( molecule_pt  molecule,
		   int          *label,
		   int          *start_counter,
		   int          bond,
		   int          atom,
		   int          back )
{
  int  *neighbors;
  int  result,
       cycles;
  int  i;

  result = 0;
  neighbors = molecule->neighbors[atom];
  if ( label[atom] == UNVISITED )
    {
      /* this is the first time we visit this node */
      label[atom] = VISITED;
      for ( i = 0; i < molecule->number_of_neighbors[atom]; i++ )
	/* visit all neighbors, but don't take direct way back */
	if ( neighbors[i] != back )
	  {
	    bond = molecule->bond_to_neighbor[atom][i];
	    cycles = check_atoms ( molecule,
				   label,
				   start_counter,
				   molecule->bond_to_neighbor[atom][i],
				   neighbors[i],
				   atom );
	    if ( cycles > 0 )
	      /* the bond to this neighbor is in a cycle */
	      if ( molecule->bonds[bond].type < CYCLE_BOND )
		molecule->bonds[bond].type += CYCLE_BOND;
	    /* add the number of cycles this atom is in due to the 
	       the current neighbor atom to the overall number of
	       cycles for the current atom */
	    result += cycles;
	  }
      if ( result > 0 )
	{
	  /* this atom is in a cycle */
	  if ( start_counter[atom] > 0 )
	    /* if the cycle has started at this atom, then don't return
	       CYCLE any longer */
	    {
	      result--;
	      start_counter[atom]--;
	    }
	  else
	    label[atom] = CYCLE;
	}
    }
  else
    if ( label[atom] == VISITED )
      /* we have been here before, so this is the first atom in the cycle */
      {
	start_counter[atom]++;
	result++;
      }
    else
      if ( label[atom] != CYCLE )
      {	
	fprintf ( stderr, "atom label: %d\n", label[atom] );
	err_panic ( "check_atoms", "unknown label" );
      }
  return result;
}

void  find_cycles ( molecule_pt  molecule )
{
  int     labels[MAX_NUMBER_OF_MOL2_ATOMS],
          start[MAX_NUMBER_OF_MOL2_ATOMS];
  int     i;

  for ( i = 0; 
	i < molecule->number_of_atoms;
	i++ )
    {
      labels[i] = UNVISITED;
      start[i] = 0;
    }
  check_atoms ( molecule, labels, start, -1, 0, -1 );
}
