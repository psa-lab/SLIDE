/*
 *  check_complementarity.c
 *                     modified by Litian He         12/01/2004
 *                     modified by Sameer Arora 	 5/17/2005
 *                     Avoid counting atom-pair as repulsive contact if both atoms
 *                     are close to a metal. Also include DONOR too for metal-bond. (no longer count DONOR, since it doesn't make sense - toneroma 06JUN07
 *					   Also since we assume that scoring happens after a bump-free docking has been achieved,
 *					   we can remove the restrictive lower-bound distance check with MIN_HBOND_LENGTH.
 *
 *  Functions in this file provide the complementarity check for
 *  interactions happened between atom pairs
 */
/* NOTE1 : toneroma 30MAY06 - in order to keep distance values and scores consistent between internal and
 *         external scoring (distances are only precise up to the rounded 3rd decimal place), all distance
 *         comparisons will be handled as follows:
 *           >= VALUE changed to >= VALUE - 0.0005
 *           <= VALUE changed to  < VALUE + 0.0005
 *            > VALUE changed to >= VALUE + 0.0005
 *            < VALUE changed to  < VALUE - 0.0005
 *         changes are also made in corresponding score_complex.c functions
 */

#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "dist_fun.h"
#include "err_handle.h"
#include "hbond_check.h"
#include "intra_hbonds.h"
#include <check_complementarity.h>


#define NON_METALBONDED_REPULSIVE /* to not count those contacts as repulsive if both parties are close to metal*/

/* toneroma 30MAY06 - redefined the function between() in order to keep internal and external scoring the same
 */
int between(float dist, float lower, float upper)
{
  return (dist >= (lower - 0.0005)) && (dist < (upper + 0.0005));
}


/* Return TRUE if salt bridge b/w two atoms, else return FALSE */
/* Not used at present */
int is_repulsive_polar( global_data_pt global, 
			atom_t lig, int lig_atom_idx, /* ligand*/
			atom_t tgt, int tgt_atom_idx, /* target*/
			float lig_tgt_dist)
{
  int act_sum = lig.act + tgt.act;
  /* Modification recommended by MIZ - Mar/29/2005. Implementing on May/15/2005
   * If both atoms are withing bonding distance of a metal, they should not be
   * considered mutually repulsive
   * Sameer May/15/2005
   */
#ifdef NON_METALBONDED_REPULSIVE
  float lig_metal_dist, tgt_metal_dist;
  int i, metal_index;
#endif

  if(between(lig_tgt_dist, MIN_HBOND_LENGTH, MAX_HBOND_LENGTH) &&
     (DONOR_DONOR == act_sum || ACCEPTOR_ACCEPTOR == act_sum)){
 
#ifdef NON_METALBONDED_REPULSIVE
    /* These potentially form repulsive contact. However, if both are bonded to
     * metal, this contact should not be counted as repulsive.  
     * SA - 16/May/2005
     */
    for(i = 0; i < global->number_of_metals; i++){
      /* MIZ-Leslie-Sameer discussions - 26/May/2005 - distance should be 
       * within MAX_HBOND_LENGTH for a polar-atom - metal interaction if we 
       * want to NOT count polar-atom pair as repulsive when both are within 
       * MAX_HBOND_LENGTH of a metal.  
       * Also since we assume that scoring happens after a bump-free docking 
       * has been achieved, we can remove the restrictive lower-bound distance 
       * check with MIN_HBOND_LENGTH.
       */
      metal_index = global->metal_atom_indices[i];
      if(global->target_atoms[metal_index].act == METAL_1 ||
	 global->target_atoms[metal_index].act == METAL_2){
        lig_metal_dist = global->target_ligand_distances[metal_index * MAX_NUMBER_OF_MOL2_ATOMS + lig_atom_idx];
        tgt_metal_dist = global->intra_target_distances[metal_index * global->number_of_target_atoms + tgt_atom_idx];
        /* NOTE1 */
        if(lig_metal_dist < MAX_HBOND_LENGTH + 0.0005 && 
           tgt_metal_dist < MAX_HBOND_LENGTH + 0.0005){  
          /* Both a and b are close to metal. So not really repulsive
             fprintf(stderr, "Metal 1: Skipping metal-induced Polar repulsive\n");*/
          return FALSE;
        } 
      }
    }
#endif
      return TRUE;
  } else return FALSE;
}

/* Not used at present */
int is_repulsive_charged( global_data_pt global, 
			  atom_t lig, int lig_atom_idx, /* ligand*/
			  atom_t tgt, int tgt_atom_idx, /* target*/
			  float lig_tgt_dist)
{
  /* here, the first atom has to be the ligand atom and 2nd is the target atom */	
  int act_sum = lig.act + tgt.act;
  
  
#ifdef NON_METALBONDED_REPULSIVE
  float lig_metal_dist, tgt_metal_dist;
  int i, metal_index;
#endif
  
  
  if ( between ( lig_tgt_dist, MIN_SALT_BRIDGE_LENGTH, MAX_SALT_BRIDGE_LENGTH ) &&
       ( ( tgt.charge > 0 && lig.charge_sum > 0 ) ||
	 ( tgt.charge < 0 && lig.charge_sum < 0 )
	 ) &&
       ( DONOR_DONOR == act_sum   ||
	 ACCEPTOR_ACCEPTOR == act_sum
	 )
       )
    {
#ifdef NON_METALBONDED_REPULSIVE
      /* These potentially form repulsive contact. However, if both are bonded to metal, 
       * this contact should not be counted as repulsive.
       * SA - 16/May/2005
       */
      
      for(i = 0; i < global->number_of_metals; i++)
	{
	  metal_index = global->metal_atom_indices[i];
	  lig_metal_dist = global->target_ligand_distances[ metal_index * MAX_NUMBER_OF_MOL2_ATOMS + lig_atom_idx];
	  tgt_metal_dist = global->intra_target_distances[ metal_index * global->number_of_target_atoms + tgt_atom_idx];
	  
	  if( global->target_atoms[metal_index].act == METAL_1 )
	    /*if ( lig_metal_dist <= MAX_METAL_1_HBOND_LENGTH && tgt_metal_dist <= MAX_METAL_1_HBOND_LENGTH )*/
	    if ( lig_metal_dist < MAX_SALT_BRIDGE_LENGTH + 0.0005 && tgt_metal_dist < MAX_SALT_BRIDGE_LENGTH + 0.0005) /* NOTE1*/
	      { /* Both a and b are close to metal. So not really repulsive
		fprintf(stderr, "Metal 1: Skipping metal-induced Charged repulsive\n");*/
		return FALSE;
	      } 
	  
	  if( global->target_atoms[metal_index].act == METAL_2 )
	    if ( lig_metal_dist < MAX_SALT_BRIDGE_LENGTH + 0.0005 && tgt_metal_dist < MAX_SALT_BRIDGE_LENGTH + 0.0005) /* NOTE1 */
	      { /* Both a and b are close to metal. So not really repulsive
		fprintf(stderr, "Metal 2: Skipping metal-induced Charged repulsive\n"); */
		return FALSE;
	      } 
	  
	}
#endif
      
      return TRUE;
    }
  else
    return FALSE;
  
}

int is_hbond(atom_t a, atom_t b, int i, int j, float sq_dist, float *angle,  
             global_data_pt global)
{
  /* Note: the <= and < are to keep the code inline with SLIDE 3.0.2 */
  if(MIN_HBOND_LENGTH_MINUS_TOL_2 <= sq_dist && 
     sq_dist < MAX_HBOND_LENGTH_PLUS_TOL_2 &&
     is_hbond_interaction(a.act, b.act) &&
     check_hydrogen_angle(global, j, i, angle)) return TRUE;
  return FALSE;
}

int is_metal_salt_bridge(atom_t lig, atom_t metal, float sq_dist)
{
  /* If a metal and a charged ligand atom are less than MAX_SALT_BRIDGE_LENGTH 
   * but ore than  MAX_METAL_?_HBOND_LENGTHs distance then count the 
   * interaction as a salt-bridge 
   * - Discussions b/w Maria, Leslie and Sameer - 26/May/05
   */
  /* NOTE1*/
  if(sq_dist < MAX_SALT_BRIDGE_LENGTH_PLUS_TOL_2 && lig.charge_sum != 0  &&
     (lig.act == ACCEPTOR || lig.act == DONEPTOR) &&
     ((metal.act == METAL_1 && sq_dist >= MAX_METAL_1_HBOND_LENGTH_PLUS_TOL_2)||
      (metal.act == METAL_2 && sq_dist >= MAX_METAL_2_HBOND_LENGTH_PLUS_TOL_2)))
    return TRUE;
       /*       ( lig.act == ACCEPTOR || lig.act == DONEPTOR || lig.act == DONOR ) removed DONOR type, since metals shouldn't interact with hydrogens.  Left in DONEPTOR for the case when it is an ACCEPTOR - toneroma 06JUN07 */
  return FALSE;
}

int is_metal_hbond(atom_t polar_atom, atom_t metal, float sq_dist)
{
  /* MIZ -recommended modification - include all polar atoms,
   * not just acceptor or doneptor for metal-bond
   * SA - 17/May/2005
   */
  /* if ( ( a.act == ACCEPTOR || a.act == DONEPTOR || a.act == DONOR ) 
   * removed DONOR type, since metals shouldn't interact with hydrogens.  
   * Left in DONEPTOR for the case when it is an ACCEPTOR 
   * toneroma 06JUN07 */  
  /* NOTE1*/
  if((polar_atom.act == ACCEPTOR || polar_atom.act == DONEPTOR ) &&
     ((metal.act == METAL_1 && sq_dist < MAX_METAL_1_HBOND_LENGTH_PLUS_TOL_2) ||
      (metal.act == METAL_2 && sq_dist < MAX_METAL_2_HBOND_LENGTH_PLUS_TOL_2)))
    return TRUE;
  return FALSE;
}

/* need to pass in charges as arguments since target atoms have a .charge
 * value and ligands use the .charge_sum variable */
int is_salt_bridge(int a_act, float a_charge, int b_act, float b_charge, 
                   float sq_dist)
{
  /* Note: the <= and < are to keep the code inline with SLIDE 3.0.2 */
  if(MIN_SALT_BRIDGE_LENGTH_MINUS_TOL_2 <= sq_dist && 
     sq_dist < MAX_SALT_BRIDGE_LENGTH_PLUS_TOL_2 &&
     is_hbond_interaction(a_act, b_act) &&
     ((b_charge > 0.0 && a_charge < 0.0 ) || 
      ( b_charge < 0.0 && a_charge > 0.0)))
    return TRUE;
  else return FALSE;
}

/* Not used at present */
int is_hphobic_complementary( atom_t a , atom_t b, float dist)
{
  float a_hydro = a.hydro - HYDROPHOB_VALUE_SHIFT;
  float b_hydro = b.hydro - HYDROPHOB_VALUE_SHIFT;
  
  /* here the assumption is that a and b do not bump */
  /* NOTE1*/
  if( dist >= HYDRO_DIST + 0.0005) return FALSE;

  /* Assuming for now that hphobic interactions are only b/w carbons - 
   * Sameer 09/Mar/05 */
  if(a.type != C || b.type != C) return FALSE;
  
  /* both a and b are hphob  */
  if ( a_hydro < 0.0 && b_hydro < 0.0 ){ 
    printf("a_hydro %f, b_hydro %f ", a_hydro, b_hydro);
    return TRUE;
  }else return FALSE;
}
