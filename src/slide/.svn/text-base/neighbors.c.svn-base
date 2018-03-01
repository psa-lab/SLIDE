#include <string.h>
#include <types.h>
#include "dist_fun.h"
#include "err_handle.h"

static inline int diff_int(int  a, int  b)
{
  return (a > b ? a - b : b - a);
}
   
#if 0
int non_covalently_bound(atom_pt atom1, atom_pt atom2)
{
  /* both atoms in the same residue, if level difference is 0, 1, 2 or 3
   * then they are covalently bound */
  if(atom1->residue_index == atom2->residue_index){
    if(diff_int(atom1->level, atom2->level) > 3 && atom1->res != HIS && 
       atom1->res != PRO && atom1->res != TYR && atom1->res != PHE && 
       atom1->res != TRP) return TRUE;
    else return FALSE;
  }
  /* atoms are main-chain atoms of neighbored residues, if they are both main 
   * chain atoms, it is assumed that they are covalently bound. this is ok 
   * since we do not change the main-chain conformation, thus there is no need 
   * to check for a bump here */
  if(diff_int ( atom1->residue_index, atom2->residue_index) == 1 && 
     atom1->level == ALPHA && atom2->level == ALPHA) return FALSE;
  /* disulfide bond */
  if(atom1->res == CYS && atom1->type == SG && 
     atom2->res == CYS && atom2->type == SG )return FALSE;
  return TRUE;
}

void identify_atom_neighbors(atom_pt atoms, int number_of_atoms)
{
  int      number_of_neighbors[MAX_PDB_ATOMS];
  float    distance;
  int      i, j;

  memset(number_of_neighbors, 0, number_of_atoms * sizeof(*number_of_neighbors));
  for(i = 0; i < number_of_atoms; i++ )
    for(j = i + 1; j < number_of_atoms; j++){
      distance = dist_fun ( atoms[i].pos, atoms[j].pos );
      if(distance < NEIGHBOR_DISTANCE && 
         non_covalently_bound ( &atoms[i], &atoms[j])){
	    atoms[i].neighbors[number_of_neighbors[i]] = j;
	    atoms[i].neighbor_dist[number_of_neighbors[i]] = distance;
	    etoms[j].neighbors[number_of_neighbors[j]] = i;
	    atoms[j].neighbor_dist[number_of_neighbors[j]] = distance;
	    number_of_neighbors[i]++;
	    number_of_neighbors[j]++;
	    if ( number_of_neighbors[i] >= MAX_NEIGHBOR_ATOMS 
		 || number_of_neighbors[j] >= MAX_NEIGHBOR_ATOMS )
	      err_panic2 ( "identify_atom_neighbors", 
			  "more than MAX_NEIGHBOR_ATOMS neighbors");
	  }
      }
  for ( i = 0; i < number_of_atoms; i++ )
    atoms[i].neighbors[number_of_neighbors[i]] = NO_MORE;
}
#endif

void  identify_water_neighbors ( atom_pt  waters,
				 atom_pt  target_atoms,
				 int      number_of_waters,
				 int      number_of_target_atoms,
				 global_data_pt global)
{
  int      number_of_neighbors[MAX_PDB_ATOMS];
  float    distance;
  int      i, j;

  for ( i = 0; i < number_of_waters; i++ )
    number_of_neighbors[i] = 0;
  for ( i = 0; i < number_of_waters; i++ )
    for ( j = 0; j < number_of_target_atoms; j++ )
      {
	distance = dist_fun ( waters[i].pos, target_atoms[j].pos );
	if ( distance < NEIGHBOR_DISTANCE 
	     && non_covalently_bound ( &waters[i], &target_atoms[j] ) )
	  {
	    waters[i].neighbors[number_of_neighbors[i]] = j;
	    waters[i].neighbor_dist[number_of_neighbors[i]] = distance;
	    number_of_neighbors[i]++;
	    if ( number_of_neighbors[i] >= MAX_NEIGHBOR_ATOMS )
	      err_panic2 ( "identify_water_neighbors", 
			  "more than MAX_NEIGHBOR_WATERS neighbors");
	  }
      }
  for ( i = 0; i < number_of_waters; i++ )
    waters[i].neighbors[number_of_neighbors[i]] = NO_MORE;
}
