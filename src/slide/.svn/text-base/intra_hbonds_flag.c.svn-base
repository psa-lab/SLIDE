/*
 * The functions are roughly same as in intra_hbonds.c file
 * The only difference is that in these functions the atoms that
 * form the hbond are recorded in an array. This was required
 * in the scoring function. -- RSK,  Aug 16, 2002
 *
 * intra_ligand_hbonds_flag() and ligand_to_water_hbonds_flag() are
 * added by Litian,  01/14/2004
 */

#include <math.h>
#include <types.h>
#include <dist_fun.h>
#include <hbond_check.h>
#include <mymalloc.h>
#include <check_complementarity.h>

/*#define PRINT_INTERACTIONS*/

#ifdef PRINT_INTERACTIONS
#include "print_interaction.h"
#endif

/*#define TRACE*/

int intra_target_polar_flag(residue_pt target_residues, atom_pt target_atoms,  
                            const int num_targ_atoms, short *target_flag, 
                            FILE *fout, float *num_of_intra_target_salt_bridges,
                            global_data_pt global)
{
  int j, k;
  float sq_dist;
  float angle;
  float preacc_angle;
  int donor;
  int acceptor;
  int sum;
  float hydrogen_position[3];
  int number_of_hbonds = 0;
  int rv = FAILURE;
  int res_idx;
  atom_pt res_start_atom;
  int num_res_atoms;

  for( j = 0; j < num_targ_atoms; j++ ){
    if(target_flag[j] != INTERFACIAL_ATOM || target_atoms[j].act == NOTHING) 
      continue;

    for(k = 0; k < num_targ_atoms; ++k){
      if(target_atoms[k].act == NOTHING) continue;

      sq_dist = squared_dist(target_atoms[j].pos, target_atoms[k].pos);
      if(is_salt_bridge(target_atoms[j].act, target_atoms[j].charge, 
                        target_atoms[k].act, target_atoms[k].charge, sq_dist)){
        target_flag[j] = INTRA_TARGET_SALT_BRIDGE;
        ++(*num_of_intra_target_salt_bridges);
#ifdef PRINT_INTERACTIONS
        print_interaction(&target_atoms[j] , TARGET, j, &target_atoms[k], 
                          TARGET, k, "intraSB", stdout);
#endif
        break;
      }
    }
    if(INTRA_TARGET_SALT_BRIDGE == target_flag[j]) continue;

    for(k = 0; k < num_targ_atoms; ++k){
      if(target_atoms[k].act == NOTHING) continue;

      sq_dist = squared_dist(target_atoms[j].pos, target_atoms[k].pos);
      if(sq_dist < MIN_HBOND_LENGTH_2 || MAX_HBOND_LENGTH_2 < sq_dist) continue;

      /* Check if this is an intramolecular hydrogen bond in the target protein 
       */
      sum = target_atoms[k].act + target_atoms[j].act;
      if(sum != ACCEPTOR_DONOR && sum != DONOR_DONEPTOR &&
         sum != ACCEPTOR_DONEPTOR && sum != DONEPTOR_DONEPTOR) continue;

      /* now compute the donor-H-acceptor angle */
      if(target_atoms[k].act == DONOR || 
         (target_atoms[k].act == DONEPTOR && target_atoms[j].act != DONOR )){
        donor = k;
        acceptor = j;
      }else{
        donor = j;
        acceptor = k;
      }
      angle = compute_target_hydrogen_angle(target_atoms, target_residues,
                                            donor, target_atoms[acceptor].pos,
                                            hydrogen_position, global);

#ifdef TRACE
      printf("%s %s %s | %s %s %s | %f %f\n", 
           target_atoms[j].residue, target_atoms[j].residue_num, 
           target_atoms[j].name, target_atoms[k].residue, 
           target_atoms[k].residue_num, target_atoms[k].name, 
           sqrtf(sq_dist), angle);
#endif

      /* Check the DHA angle */
      if(angle < MIN_HYDROGEN_ANGLE || MAX_HYDROGEN_ANGLE < angle) continue;

      res_idx = target_atoms[acceptor].residue_index;
      res_start_atom = &(target_atoms[target_residues[res_idx].start_atom]);
      num_res_atoms = target_residues[res_idx].number_of_atoms;
      rv = target_preacc_angle(&target_atoms[acceptor], hydrogen_position,
                               res_start_atom, num_res_atoms, &preacc_angle);

      if(rv != SUCCESS || preacc_angle < MIN_PREACCEPTOR_ANGLE ||
         MAX_PREACCEPTOR_ANGLE < preacc_angle) continue;

      number_of_hbonds++;
      target_flag[donor] = INTRA_TARGET_HBOND;
      target_flag[acceptor] = INTRA_TARGET_HBOND;

#ifdef PRINT_INTERACTIONS
      print_interaction(&target_atoms[k], TARGET, k, &target[j], TARGET, j, 
                        "intraHB", global, fout);
#endif
    }
  }
  return number_of_hbonds;
}

/* the follows are added by litian,  01/14/2004*/
void intra_ligand_hbonds_flag ( global_data_pt  global,
			       int             index,
			       short           ligand_flag[MAX_PDB_ATOMS] )
{
  atom_pt ligand;
  float   angle, sq_dist;
  int     donor,
          acceptor,
          sum,
          number_of_hbonds;
  int     j;

  ligand = global->ligand->atoms;
  number_of_hbonds = 0;

  for(j = 0; j < global->ligand->number_of_atoms; j++){
    if(j == index || ligand[j].act == NOTHING) continue;

    sq_dist = squared_dist ( ligand[index].pos, ligand[j].pos ); 
    /* since we only want to check for hydrogen bonds, only DONOR, 
     * ACCEPTOR, and DONEPTOR atoms are of interest */
    if(MIN_HBOND_LENGTH_2 <= sq_dist  && sq_dist <= MAX_HBOND_LENGTH_2){ 
      sum = ligand[index].act + ligand[j].act; 
      /* this is an intramolecular hydrogen bond in the ligand protein */
      if(sum == ACCEPTOR_DONOR || sum == DONOR_DONEPTOR || 
         sum == ACCEPTOR_DONEPTOR || sum == DONEPTOR_DONEPTOR){ 
        if((ligand[index].act == DONOR  || ligand[index].act == DONEPTOR) && 
           ligand[j].act != DONOR ){ 
          donor = index; 
          acceptor = j; 
        }else{ 
          donor = j; 
          acceptor = index; 
        } 
          
        /* now compute the donor-H-acceptor angle  */ 
        angle = compute_ligand_hydrogen_angle(global, donor, 
                                              ligand[acceptor].pos, &angle );
        if ( MIN_HYDROGEN_ANGLE < angle && angle < MAX_HYDROGEN_ANGLE ){ 
          number_of_hbonds++; 
          ligand_flag[index] = INTRA_LIGAND_HBOND;
        } 
      } 
    }
  }
}


#if 0
int intra_target_polar_flag(residue_pt target_residues, atom_pt target_atoms,  
                            const int num_targ_atoms, octree_p target_octree, 
                            short *target_flag, FILE *fout, 
                            float *num_of_intra_target_salt_bridges,
                            global_data_pt global)
{
  int j;
  float sq_dist;
  float angle;
  atom_pt *close_atoms;
  atom_pt *close_atom_p;
  atom_pt close_atom;
  size_t num_atoms;
  int donor;
  int acceptor;
  int sum;
  float hydrogen_position[3];
  int number_of_hbonds = 0;

  for( j = 0; j < num_targ_atoms; j++ ){
    if(target_flag[j] != INTERFACIAL_ATOM || target_atoms[j].act == NOTHING) 
      continue;

    get_close_atoms(&global->target_octree, target_atoms[j].pos,
                    DONT_CARE_BUMP_DISTANCE, &close_atoms, &num_atoms);
    close_atom_p = close_atoms;
    for( ; close_atom_p < &close_atoms[num_atoms]; ++close_atom_p){
      close_atom = *close_atom_p;
      sq_dist = squared_dist(target_atoms[j].pos, close_atom->pos);
      if(is_salt_bridge(target_atoms[j].act, target_atoms[j].charge, 
                        close_atom->act, close_atom->charge, sq_dist)){
        target_flag[j] = INTRA_TARGET_SALT_BRIDGE;
        ++(*num_of_intra_target_salt_bridges);
#ifdef PRINT_INTERACTIONS
        k = close_atom - target_atoms;
        print_interaction(&target_atoms[j] , TARGET, j, close_atom, TARGET, k, 
                          "intraSB", stdout);
#endif
          break;
      }
    }
    if(INTRA_TARGET_SALT_BRIDGE == target_flag[j]){
      if(close_atoms) my_free(close_atoms);
      continue;
    }

    close_atom_p = close_atoms;
    for( ; close_atom_p < &close_atoms[num_atoms]; ++close_atom_p){
      close_atom = *close_atom_p;
      sq_dist = squared_dist(target_atoms[j].pos, close_atom->pos);
      if(sq_dist < MIN_HBOND_LENGTH_2 || MAX_HBOND_LENGTH_2 < sq_dist){
        continue;
      }

      /* Check if this is an intramolecular hydrogen bond in the target protein 
       */
      sum = close_atom->act + target_atoms[j].act;
      if(sum != ACCEPTOR_DONOR && sum != DONOR_DONEPTOR &&
         sum != ACCEPTOR_DONEPTOR && sum != DONEPTOR_DONEPTOR) continue;

      /* now compute the donor-H-acceptor angle */
      if(close_atom->act == DONOR || 
         (close_atom->act == DONEPTOR && target_atoms[j].act != DONOR )){
        donor = close_atom - target_atoms;
        acceptor = j;
      }else{
        donor = j;
        acceptor = close_atom - target_atoms;
      }
      angle = compute_target_hydrogen_angle(target_atoms, target_residues,
                                            donor, target_atoms[acceptor].pos,
                                            hydrogen_position, global);
#ifdef TRACE
      printf("%s %s %s | %s %s %s | %f %f\n", 
           close_atom->residue, close_atom->residue_num, close_atom->name,
           target_atoms[j].residue, target_atoms[j].residue_num, 
           target_atoms[j].name, sqrtf(sq_dist), angle);
#endif
      if(MIN_HYDROGEN_ANGLE < angle && angle < MAX_HYDROGEN_ANGLE){
        number_of_hbonds++;
        target_flag[donor] = INTRA_TARGET_HBOND;
        target_flag[acceptor] = INTRA_TARGET_HBOND;
#ifdef PRINT_INTERACTIONS
      print_interaction(close_atom, TARGET, close_atom - target_atoms, 
                        &target[j], TARGET, j, "intraHB", global, fout);
#endif
      }
    }

    if(close_atoms) my_free(close_atoms);
  }
  return number_of_hbonds;
}

int intra_target_hbonds_flag(global_data_pt  global, int index,
                             short target_flag[MAX_PDB_ATOMS], FILE *fout)
{
  float   angle, sq_dist;
  int     donor,
          acceptor,
          sum;
  int     j;
  float   hydrogen_position[3];
  atom_pt target = global->target_atoms;
  int number_of_hbonds = 0;
  
  for ( j = 0; j < global->number_of_target_atoms; j++ ){
    /* since we only want to check for hydrogen bonds, only DONOR,
       ACCEPTOR, and DONEPTOR atoms are of interest */
    /* Sameer 11/28/04 - check for min distance too*/
    if(j == index || target[j].act == NOTHING) continue;

    /* This should be precomputed, but needs to be done carefully since we do 
     * allow the target atoms to move.  If nothing else, a smart space 
     * partioning would be useful.  Finally, the whole SLIDE code would need
     * to be revamped and the idea might not be readily apparent, but computing 
     * sqrt() all the time is needless computation since we can just compare the
     * squared distance with respect to square tolerances.  */
    /* besides at this point it would be faster to port the unbumping and 
     * scoring plus miscellaneous code to the SimSite3D code base, rather
     * than muck around in the SLIDE code */
    sq_dist = squared_dist(target[index].pos, target[j].pos);
    if(sq_dist < MIN_HBOND_LENGTH_2 || MAX_HBOND_LENGTH_2 < sq_dist) continue;

    /* Check if this is an intramolecular hydrogen bond in the target protein */
    sum = target[index].act + target[j].act;
    if(sum != ACCEPTOR_DONOR && sum != DONOR_DONEPTOR &&
       sum != ACCEPTOR_DONEPTOR && sum != DONEPTOR_DONEPTOR) continue;

    if(target[index].act == DONOR 
       || (target[index].act == DONEPTOR && target[j].act != DONOR )){
      donor = index;
      acceptor = j;
    }else{
      donor = j;
      acceptor = index;
    }
    /* now compute the donor-H-acceptor angle */
    angle = compute_target_hydrogen_angle(target, global->target_residues,
                                          donor, target[acceptor].pos,
                                          hydrogen_position, global);
#ifdef TRACE
    printf("%s %s %s | %s %s %s | %f %f\n", 
           target[index].residue, target[index].residue_num, target[index].name,
           target[j].residue, target[j].residue_num, target[j].name,
           sqrtf(sq_dist), angle);
#endif
    if(MIN_HYDROGEN_ANGLE < angle && angle < MAX_HYDROGEN_ANGLE){
      number_of_hbonds++;
      target_flag[index] = INTRA_TARGET_HBOND;
      target_flag[j] = INTRA_TARGET_HBOND;
#ifdef PRINT_INTERACTIONS
      print_interaction(&target[index], TARGET, index, &target[j], TARGET, j,
                        "intraHB", global, fout);
#endif
    }
#ifdef TRACE
#endif
  }
  return number_of_hbonds;
}
#endif

void target_to_water_hbonds_flag ( global_data_pt  global,
				  int             index,
				  short           target_flag[MAX_PDB_ATOMS] )
{
    atom_pt target,
	    waters;
    float   angle, sq_dist;
    float   number_of_hbonds;
    int     j;
    float   hydrogen_position[3];

  /* no binding-site waters */
  if(global->number_of_waters == 0) return;

  target = global->target_atoms;
  waters = global->waters;
  number_of_hbonds = 0.0;

  for( j = 0; j < global->number_of_waters; j++ ){
    sq_dist = squared_dist ( target[index].pos, waters[j].pos );

    /* Sameer 11/28/04 - check for min distance too*/
    if(waters[j].state == CONSERVED && 
       MIN_HBOND_LENGTH_2 <= sq_dist && sq_dist <= MAX_HBOND_LENGTH_2 ){ 
      /* compute the donor-H-acceptor angle */ 
      if(target[index].act == DONOR ) { 
        angle = compute_target_hydrogen_angle(target, global->target_residues, 
                                              index, waters[j].pos, 
                                              hydrogen_position, global);
	      if ( MIN_HYDROGEN_ANGLE < angle  && angle < MAX_HYDROGEN_ANGLE )
		{
		  number_of_hbonds += waters[j].prediction;
		  if ( waters[j].prediction > 0.75 )
		    target_flag[index] = TARGET_WATER_HBOND;
		}
#ifdef TRACE
	      printf ("index = %d - %d (%7.3f,%7.3f,%7.3f), angle = %7.2f\n",
		      index, j, waters[j].pos[X], waters[j].pos[Y],
		      waters[j].pos[Z], angle );
#endif
	    }
	  else if (target[index].act == DONEPTOR || target[index].act == ACCEPTOR   )
	    { /* Sameer 28/Nov/04 - What if target is HYDROPHOBIC*/
	      /* target is either DONEPTOR or ACCEPTOR, in both cases
		 we already have a hydrogen bond, since we only consider
		 the distance between this atom and the water */

	      number_of_hbonds += waters[j].prediction;
	      if ( waters[j].prediction > 0.75 )
		target_flag[index] = TARGET_WATER_HBOND;
#ifdef TRACE
	      printf ("index = %d - %d (%7.3f, %7.3f,%7.3f)\n",
		      index, j, waters[j].pos[X], waters[j].pos[Y],
		      waters[j].pos[Z] );
#endif
	    }
	}
      /*return (int) number_of_hbonds;*/
	} /*end of for loop*/
}


void ligand_to_water_hbonds_flag ( global_data_pt  global,
				   int             index,
				   short           ligand_flag[MAX_PDB_ATOMS] )
{
  atom_pt  ligand, waters;
  float   angle, sq_dist;
  float   number_of_hbonds;
  int     j;

  /* no binding-site waters */
  if ( global->number_of_waters == 0 ) return;

  ligand = global->ligand->atoms;
  waters = global->waters;
  number_of_hbonds = 0.0;

  /*for ( i = 0; i < global->number_of_ligand_atoms; i++ )*/
  /*  if ( ligand[i].act != NOTHING )*/

  for ( j = 0; j < global->number_of_waters; j++ )
    if ( waters[j].state == CONSERVED ){
		sq_dist = squared_dist ( ligand[index].pos, waters[j].pos );
	if ( MIN_HBOND_LENGTH_2 <= sq_dist && sq_dist <= MAX_HBOND_LENGTH_2 )
	  {
	    if ( ligand[index].act == DONOR )
	      /* compute the donor-H-acceptor angle */
	      {
		angle =
		  compute_ligand_hydrogen_angle ( global,
						  index,
						  waters[j].pos,
						  &angle );
		if ( MIN_HYDROGEN_ANGLE < angle && angle < MAX_HYDROGEN_ANGLE )
		  {
		    number_of_hbonds += waters[j].prediction;
		    if ( waters[j].prediction > 0.75 )
		      ligand_flag[index] = LIGAND_WATER_HBOND;
		  }
	      }
	    else if (ligand[index].act == DONEPTOR || ligand[index].act == ACCEPTOR   )
	      { /* Sameer 28/Nov/04 - What if ligand is HYDROPHOBIC*/
		/* ligand is either DONEPTOR or ACCEPTOR, in both cases
		   we already have a hydrogen bond, since we only consider
		   the distance between this atom and the water */

		number_of_hbonds += waters[j].prediction;
		if ( waters[j].prediction > 0.75 )
		  ligand_flag[index] = LIGAND_WATER_HBOND;
	      }
	  }
      }
  /*return (int) number_of_hbonds;*/
}

