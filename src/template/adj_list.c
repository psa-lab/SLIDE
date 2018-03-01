#include "types.h"
#include "defs.h"

void  construct_adjacency_list ( molecule_pt molecule )
{
  bond_pt  bonds;
  int      **neighbors;
  int      **bond;
  int      *number_of_neighbors;
  int      i, atom1, atom2;
  
  number_of_neighbors = molecule->number_of_neighbors;
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    number_of_neighbors[i] = 0;
  bonds = molecule->bonds;
  neighbors = molecule->neighbors;
  number_of_neighbors = molecule->number_of_neighbors;
  bond = molecule->bond_to_neighbor;
  for ( i = 0; i < molecule->number_of_bonds; i++ )
    {
      atom1 = bonds[i].atom1;
      atom2 = bonds[i].atom2;
      neighbors[atom1][number_of_neighbors[atom1]] = atom2;
      bond[atom1][number_of_neighbors[atom1]++] = i;
      neighbors[atom2][number_of_neighbors[atom2]] = atom1;
      bond[atom2][number_of_neighbors[atom2]++] = i;
    }
}
