#include <math.h>
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "dist_fun.h"
#include "unbump_water.h"
#include <check_complementarity.h>
#include <find_all_bumps.h>

/*#define  TRACE*/

/* Seems to tally total_overlap */
int find_all_bumps(dock_feats_pt features, global_data_pt global)
{
  atom_pt atom;
  atom_pt water;
  atom_pt *neighbors;
  float distance;
  float sq_dist;
  float overlap_tolerance = global->side_chain_overlap;
  float overlap;
  float tmp;
  int number_of_bumps = 0,
           sum,
           dist_int;
  int      i, j, k;

  atom_pt target = global->target_atoms;
  atom_pt ligand = global->ligand->atoms;
  int *target_bumps = global->target_bumps;
  int *ligand_bumps = global->ligand_bumps;
  const int num_targ_atoms = global->number_of_target_atoms;
  features->total_overlap = 0.0;

#ifdef TRACE
  printf("\n");
#endif

  for ( i = 0; i < num_targ_atoms; i++ ){
    for ( j = 0; j < global->ligand->number_of_atoms; j++ ){
      /* Ligand hydrogen atoms are considered if there were not added by SLIDE 
       * */
      if(ligand[j].orbit == ADDED) continue;

      sq_dist = global->target_ligand_sq_dists[j*num_targ_atoms + i];
      if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

      if(atoms_overlap2(&target[i], &ligand[j], sq_dist, overlap_tolerance, 
                        &overlap, &distance)){
#ifdef TRACE
        printf("sum=%d, distance=%f, ligand[j].rad=%f, target[i].rad=%f, "
               "tolerated overlap=%f, remaining overlap=%f\n", 
               target[i].act + ligand[j].act, distance, ligand[j].rad, 
               target[i].rad, overlap_tolerance, overlap);
        printf("cumulative overlap = %f\n", features->total_overlap);
        printf("BUMP no. %2d: %3s %s %3s (%7.3f %7.3f %7.3f) <--> %2d %4s %2d "
               "(%7.3f %7.3f %7.3f) : %7.3f (%7.3f)\n", number_of_bumps,
               target[i].name, target[i].residue, target[i].residue_num,
               target[i].pos[X], target[i].pos[Y], target[i].pos[Z],
               j + 1, ligand[j].type_str, ligand[j].fragment,
               ligand[j].pos[X], ligand[j].pos[Y], ligand[j].pos[Z],
               distance, overlap);
#endif
        /* the vector `bump_point' is already full, thus there is no hope for 
         * this ligand ==> bump-check failed */
        if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS ) {
          global->number_of_bumps = number_of_bumps;
          return FAILURE;
        }
        target_bumps[number_of_bumps] = i;
        ligand_bumps[number_of_bumps] = j;
        features->total_overlap += overlap;
        number_of_bumps++;
      }
    }
  }
#ifdef TRACE
  printf("cumulative overlap = %f\n", features->total_overlap);
#endif

#if 0
  /* This portion of the code could suffer from bit rot -- water handling
   * has not been used for some time */

  for(i = 0; i < global->number_of_waters; i++){
    water = &global->waters[i];
    if(water->state != CONSERVED) continue;

    /* find all bumps of this water with any ligand atom */
    for(j = 0; j < global->ligand->number_of_atoms; j++){
      if(ligand[j].orbit == ADDED) continue;

      distance = dist_fun(water->pos, ligand[j].pos);
      if(distance > DONT_CARE_BUMP_DISTANCE) continue;

      water_atom_overlap(&ligand[j], distance, overlap_tolerance, water);
      if(water->state == WATER_OVERLAP_NEEDS_TO_BE_RESOLVED){
        water->state = CONSERVED;
        if(unbump_water(global, water) == FAILURE){	      
#ifdef TRACE
          printf("BUMP no. %2d: water %d <--> %2d %4s %2d : %7.3f (%7.3f)\n",
          number_of_bumps, water->number, j + 1,
                 ligand[j].type_str, ligand[j].fragment, distance,
                 WATER_RAD + ligand[j].rad - overlap_tolerance - distance);
#endif
          /* the vector `bump_point' is already full, thus there is no hope for
           * this ligand ==> bump-check failed */
          if(number_of_bumps >= MAX_SIDE_CHAIN_BUMPS){
            global->number_of_bumps = number_of_bumps;
            return FAILURE;
          }
          target_bumps[number_of_bumps] = (-1) * ( i + 1 );
          ligand_bumps[number_of_bumps] = j;
          number_of_bumps++;
        }
      }
    }

    /* find all bumps of this water with any target atom */
    if(water->state != CONSERVED) continue;
    for(j = 0; j < global->number_of_target_atoms; j++){
      distance = dist_fun(water->pos, target[j].pos);
      if(distance > DONT_CARE_BUMP_DISTANCE) continue;

      water_atom_overlap(&target[j], distance, overlap_tolerance, water);
      if(water->state == WATER_OVERLAP_NEEDS_TO_BE_RESOLVED){
        if(unbump_water(global, water) == FAILURE) water->state = DISPLACED;
        else water->state = CONSERVED;
      }
    }
  }
#endif

  overlap_tolerance = global->intra_overlap;
  for(i = 0; i < global->number_of_target_residues; i++){
    if(global->target_intra_overlap[i] != YES) continue;
     
    for(j = 0; j < global->target_residues[i].number_of_atoms; j++){
      atom = &target[global->target_residues[i].start_atom+j];
      neighbors = atom->neighbors;

      for(k = 0; k < atom->num_nbrs + atom->num_added_nbrs; ++k){

        sq_dist = squared_dist(atom->pos, neighbors[k]->pos);
        if(sq_dist > DONT_CARE_BUMP_DISTANCE_SQ) continue;

        tmp = atom->neighbor_dist[k] - overlap_tolerance;
        /* Added the first inequality to keep this consistent with previous
         * slide code */
        if(sq_dist < tmp * tmp &&
           atoms_overlap2(atom, neighbors[k], sq_dist, overlap_tolerance, 
                          &overlap, &distance)){
          /* the vector `bump_point' is already full, thus there is no hope for 
           * this ligand ==> bump-check failed */
          if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS ) {
            global->number_of_bumps = number_of_bumps;
            return FAILURE;
          }
          features->total_overlap += overlap;
          target_bumps[number_of_bumps] = 
            global->target_residues[i].start_atom+j;
          ligand_bumps[number_of_bumps] = (-1) * ((neighbors[k] - target) + 1);
          number_of_bumps++;
#ifdef TRACE
          printf("intra BUMP %s %s %s (%7.3f,%d) with %s %s %s (%7.3f,%d): %f "
                 "(%f), orig = %f \n", atom->name, atom->residue,
                 atom->residue_num, atom->rad, atom->act, neighbors[k]->name, 
                 neighbors[k]->residue, neighbors[k]->residue_num, 
                 neighbors[k]->rad, neighbors[k]->act, distance,
                 atom->rad + neighbors[k]->rad - overlap_tolerance - distance,
                 atom->neighbor_dist[k]);
          printf("cumulative overlap = %f\n", features->total_overlap);
#endif
        }
      }
    }
  }

  global->number_of_bumps = number_of_bumps;
  if(number_of_bumps == 0) return NO_BUMP;
  return BUMP;
}

int 
atoms_overlap(const atom_pt a, const atom_pt b, const float distance, 
              const float overlap_tolerance, float *overlap)
{
  *overlap = 0.0;
  if(is_hbond_interaction(a->act, b->act)){ 
    if(distance < MIN_HBOND_LENGTH) *overlap = MIN_HBOND_LENGTH - distance;
  }else if((a->act + b->act == ACCEPTOR_DONOR_HYDROGEN ||
          a->act + b->act == DONEPTOR_DONOR_HYDROGEN)){
    if(distance < MIN_HBOND_LENGTH_HYDROGEN) 
      *overlap = MIN_HBOND_LENGTH_HYDROGEN - distance;
  }else if(a->act == METAL_1 || b->act == METAL_1){
    if(distance < MIN_METAL_1_HBOND_LENGTH)
      *overlap = MIN_METAL_1_HBOND_LENGTH - distance;
  }else if(a->act == METAL_2 || b->act == METAL_2){
    if(distance < MIN_METAL_2_HBOND_LENGTH)
      *overlap = MIN_METAL_2_HBOND_LENGTH - distance;
  }else if(distance < a->rad + b->rad - overlap_tolerance) 
    *overlap = a->rad + b->rad - overlap_tolerance - distance;

  if(*overlap > 0.0) return 1;
  return 0;
} 

int
atoms_overlap2(const atom_pt a, const atom_pt b, const float sq_dist, 
               const float overlap_tolerance, float *overlap, float *dist)
{
  float tmp;
  *overlap = 0.0;

  if(is_hbond_interaction(a->act, b->act)){ 
    if(sq_dist < MIN_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_HBOND_LENGTH - *dist;
    }
  }else if((a->act + b->act == ACCEPTOR_DONOR_HYDROGEN ||
          a->act + b->act == DONEPTOR_DONOR_HYDROGEN)){
    if(sq_dist < MIN_HBOND_LENGTH_HYDROGEN_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_HBOND_LENGTH_HYDROGEN - *dist;
    }
  }else if(a->act == METAL_1 || b->act == METAL_1){
    if(sq_dist < MIN_METAL_1_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_METAL_1_HBOND_LENGTH - *dist;
    }
  }else if(a->act == METAL_2 || b->act == METAL_2){
    if(sq_dist < MIN_METAL_2_HBOND_LENGTH_2){
      *dist = sqrt(sq_dist);
      *overlap = MIN_METAL_2_HBOND_LENGTH - *dist;
    }
  }else{
    tmp = a->rad + b->rad - overlap_tolerance;
    if(sq_dist < tmp*tmp){
      *dist = sqrt(sq_dist);
      *overlap = tmp - *dist;
    }
  }

  if(*overlap > 0.0) return 1;
  return 0;
}

int 
water_atom_overlap(const atom_pt a, const float distance, 
                   const float overlap_tolerance, atom_pt w)
{
  /* there is a bump of this water with a polar ligand atom, so
 can be displaced at not cost */
  if((a->act == ACCEPTOR || a->act == DONEPTOR || a->act == DONOR) && 
     distance < MIN_HBOND_LENGTH)
    w->state = POLAR_DISPLACED;
  /* the ligand atom is pretty close to the water, so we just
     assume that the water will be displaced */
  else if(a->act == NOTHING && 
          distance < WATER_RAD / 2.0 + a->rad - overlap_tolerance)
    w->state = DISPLACED;
  /* there is only a relatively small overlap, so rather than assuming we
   * can displace the water, we attempt to "unbump" it */
  else if(a->act == NOTHING && 
          distance < WATER_RAD + a->rad - overlap_tolerance)
    w->state = WATER_OVERLAP_NEEDS_TO_BE_RESOLVED;

  return SUCCESS;
}

#if 0 
int  find_all_bumps_output ( global_data_pt  global )
{
  atom_pt  target,
           ligand,
           atom,
           water;
  int      *ligand_bumps,
           *target_bumps,
           *neighbors;
  float    distance,
           overlap;
  int      number_of_bumps,
           sum;
  int      i, j, k;
/***** Added by PCS -- 22-Mar-00 *****/
  char     filename[FILENAME_MAX];
/*************************************/

  target = global->target_atoms;
  ligand = global->ligand->atoms;
  target_bumps = global->target_bumps;
  ligand_bumps = global->ligand_bumps;
  overlap = global->side_chain_overlap;
  global->total_overlap = 0.0;
  number_of_bumps = 0;

  for ( i = 0; i < global->number_of_target_atoms; i++ )
    for ( j = 0; j < global->ligand->number_of_atoms; j++ )
      if ( ligand[j].orbit != ADDED )
	{
	  distance 
	    = global->target_ligand_distances[i*MAX_NUMBER_OF_MOL2_ATOMS+j];
	  if ( distance > DONT_CARE_BUMP_DISTANCE )
	    continue;
	  sum = target[i].act + ligand[j].act;
	  if ( ( sum != ACCEPTOR_DONOR
		 && sum != ACCEPTOR_DONEPTOR
		 && sum != DONOR_DONEPTOR
		 && sum != DONEPTOR_DONEPTOR
		 && sum != ACCEPTOR_DONOR_HYDROGEN
		 && sum != DONEPTOR_DONOR_HYDROGEN
		 && target[i].act != METAL_1
		 && target[i].act != METAL_2
		 && distance < ligand[j].rad 
		 + target[i].rad 
		 - overlap  )
	       || ( ( sum == ACCEPTOR_DONOR
		      || sum == DONOR_DONEPTOR
		      || sum == ACCEPTOR_DONEPTOR
		      || sum == DONEPTOR_DONEPTOR )
		    && distance < MIN_HBOND_LENGTH ) 
	       || ( ( sum == ACCEPTOR_DONOR_HYDROGEN
		      || sum == DONEPTOR_DONOR_HYDROGEN ) 
		    && distance < MIN_HBOND_LENGTH_HYDROGEN ) 
	       || ( ( target[i].act == METAL_1 )
		    && ( distance < ligand[j].rad
			 + target[i].rad
			 - overlap )
		    && ( distance < MIN_METAL_1_HBOND_LENGTH ) )
	       || ( ( target[i].act == METAL_2 )
		    && ( distance < ligand[j].rad
			 + target[i].rad
			 - overlap )
		    && ( distance < MIN_METAL_2_HBOND_LENGTH ) )
	       )
	    {
	      if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS )
		/* the vector `bump_point' is already full, thus there
		   is no hope for this ligand ==> bump-check failed */
		{
		  global->number_of_bumps = number_of_bumps;
/***** Added by PCS -- 22-Mar-00 *****/
#ifndef OUTPUT_ALL_BUMPS
		  return FAILURE;
#endif
/*************************************/
		}
	      target_bumps[number_of_bumps] = i;
	      ligand_bumps[number_of_bumps] = j;
	      if ( sum != ACCEPTOR_DONOR
		   && sum != ACCEPTOR_DONEPTOR
		   && sum != DONOR_DONEPTOR
		   && sum != DONEPTOR_DONEPTOR
		   && sum != ACCEPTOR_DONOR_HYDROGEN
		   && sum != DONEPTOR_DONOR_HYDROGEN
		   && target[i].act != METAL_1
		   && target[i].act != METAL_2
		   && distance < ligand[j].rad 
		   + target[i].rad 
		   - overlap  )
		global->total_overlap +=
		  ligand[j].rad + target[i].rad - overlap - distance;
	      if ( ( sum == ACCEPTOR_DONOR
		     || sum == DONOR_DONEPTOR
		     || sum == ACCEPTOR_DONEPTOR
		     || sum == DONEPTOR_DONEPTOR )
		   && distance < MIN_HBOND_LENGTH ) 
		global->total_overlap += 
		  MIN_HBOND_LENGTH - distance;
	      if ( ( sum == ACCEPTOR_DONOR_HYDROGEN
		     || sum == DONEPTOR_DONOR_HYDROGEN ) 
		   && distance < MIN_HBOND_LENGTH_HYDROGEN )
		global->total_overlap += 
		  MIN_HBOND_LENGTH_HYDROGEN - distance;
	      if ( ( target[i].act == METAL_1 )
		   && ( distance < ligand[j].rad
			+ target[i].rad
			- overlap )
		   && ( distance < MIN_METAL_1_HBOND_LENGTH ) )
		{
		  global->total_overlap +=
		    MIN_METAL_1_HBOND_LENGTH - distance;
		}
	      if ( ( target[i].act == METAL_2 )
		   && ( distance < ligand[j].rad
			+ target[i].rad
			- overlap )
		   && ( distance < MIN_METAL_2_HBOND_LENGTH ) )
		{
		  global->total_overlap +=
		    MIN_METAL_2_HBOND_LENGTH - distance;
		}
	      number_of_bumps++;
	      printf ( "BUMP no. %2d: %3s %s %3s (%7.3f %7.3f %7.3f) <--> %2d %6s (%7.3f %7.3f %7.3f) : %7.3f (%7.3f)\n",
		       number_of_bumps-1,
		       target[i].name,
		       target[i].residue,
		       target[i].residue_num,
		       target[i].pos[X],
		       target[i].pos[Y],
		       target[i].pos[Z],
		       j + 1,
		       ligand[j].type_str,
		       ligand[j].pos[X],
		       ligand[j].pos[Y],
		       ligand[j].pos[Z],
		       distance,
		       target[i].rad + ligand[j].rad - overlap - distance );
	    }
	}

  for ( i = 0; i < global->number_of_waters; i++ )
    {
      water = &global->waters[i];
      if ( water->state == CONSERVED )
	/* find all bumps of this water with any ligand atom */
	for ( j = 0; j < global->ligand->number_of_atoms; j++ )
	  if ( ligand[j].orbit != ADDED )
	    {
	      distance = dist_fun ( water->pos, 
				    ligand[j].pos );
	      if ( distance > DONT_CARE_BUMP_DISTANCE )
		continue;
	      if ( ( ligand[j].act == ACCEPTOR
		     || ligand[j].act == DONEPTOR
		     || ligand[j].act == DONOR )
		   && distance < MIN_HBOND_LENGTH )
	      /* there is a bump of this water with a polar ligand atom, so
		 it can be displaced at not cost */
		water->state = POLAR_DISPLACED;
	      else
		if ( ligand[j].act == NOTHING
		     && distance < WATER_RAD / 2.0
		     + ligand[j].rad 
		     - overlap  )
		  /* the ligand atom is pretty close to the water, so we just
		     assume that the water will be displaced */
		  water->state = DISPLACED;
		else
		  if ( ligand[j].act == NOTHING
		       && distance < WATER_RAD
		       + ligand[j].rad 
		       - overlap  )	      
		    /* there is only a relatively small overlap, so don't
		       replace the water */
		    if ( unbump_water ( global, 
					&global->waters[i] ) == FAILURE )
		      {	      
			if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS )
			  /* the vector `bump_point' is already full, thus 
			     there is no hope for this ligand ==> bump-check 
			     failed */
			  {
			    global->number_of_bumps = number_of_bumps;
/***** Added by PCS -- 22-Mar-00 *****/
#ifndef OUTPUT_ALL_BUMPS
			    return FAILURE;
#endif
/*************************************/
			  }
			printf ( "BUMP no. %2d: water %d <--> %2d %4s %2d : %7.3f (%7.3f)\n",
				 number_of_bumps,
				 global->waters[i].number,
				 j + 1,
				 ligand[j].type_str,
				 ligand[j].fragment,
				 distance,
				 WATER_RAD + ligand[j].rad - overlap 
				 - distance );
			target_bumps[number_of_bumps] = (-1) * ( i + 1 );
			ligand_bumps[number_of_bumps] = j;
			number_of_bumps++;
		      }
	    }
      if ( water->state == CONSERVED )
	/* find all bumps of this water with any target atom */
	for ( j = 0; j < global->number_of_target_atoms; j++ )
	  {
	    distance = dist_fun ( water->pos, 
				  target[j].pos );
	    if ( distance > DONT_CARE_BUMP_DISTANCE )
	      continue;
	    if ( ( target[j].act == ACCEPTOR
		   || target[j].act == DONEPTOR
		   || target[j].act == DONOR )
		 && distance < MIN_HBOND_LENGTH )
	      /* there is a bump of this water with a polar target atom, so
		 it can be displaced at not cost */
	      water->state = POLAR_DISPLACED;
	    else
	      if ( target[j].act == NOTHING
		   && distance < WATER_RAD / 2.0
		   + target[j].rad 
		   - overlap  )
		/* the ligand atom is pretty close to the water, so we just
		   assume that the water will be displaced */
		water->state = DISPLACED;
	      else
		if ( target[j].act == NOTHING
		     && distance < WATER_RAD
		     + target[j].rad 
		     - overlap  )	      
		  /* there is only a relatively small overlap, so don't
		     replace the water, try to resolve the collision */
		  if ( unbump_water ( global, 
				      &global->waters[i] ) == FAILURE )
		    water->state = DISPLACED;
	  }
    }
  overlap = global->intra_overlap;

  for ( i = 0; i < global->number_of_target_residues; i++ )
    if ( global->target_intra_overlap[i] == YES )
      {
	for ( j = 0; j < global->target_residues[i].number_of_atoms; j++ )
	  {
	    atom = &target[global->target_residues[i].start_atom+j];
	    neighbors = atom->neighbors;
	    k = 0;
	    while ( neighbors[k] != NO_MORE )
	      {
		distance = dist_fun ( atom->pos, target[neighbors[k]].pos );
		if ( distance > DONT_CARE_BUMP_DISTANCE )
		  {
		    k++;
		    continue;
		  }
		sum = atom->act + target[neighbors[k]].act;
		if ( distance < atom->neighbor_dist[k] - overlap
		     && ( ( sum != ACCEPTOR_DONOR
			    && sum != ACCEPTOR_DONEPTOR
			    && sum != DONOR_DONEPTOR
			    && sum != DONEPTOR_DONEPTOR
			    && sum != ACCEPTOR_DONOR_HYDROGEN
			    && sum != DONEPTOR_DONOR_HYDROGEN
			    && atom->act != METAL_1
			    && atom->act != METAL_2
			    && target[neighbors[k]].act != METAL_1
			    && target[neighbors[k]].act != METAL_2
			    && distance < atom->rad 
			    + target[neighbors[k]].rad 
			    - overlap )
			  || ( ( sum == ACCEPTOR_DONOR
				 || sum == DONOR_DONEPTOR
				 || sum == ACCEPTOR_DONEPTOR
				 || sum == DONEPTOR_DONEPTOR )
			       && distance < MIN_HBOND_LENGTH ) 
			  || ( ( sum == ACCEPTOR_DONOR_HYDROGEN  
				 || sum == DONEPTOR_DONOR_HYDROGEN ) 
			       && distance < MIN_HBOND_LENGTH_HYDROGEN ) 
			  || ( ( atom->act == METAL_1 
				 || target[neighbors[k]].act == METAL_1 )
			       && ( distance < target[neighbors[k]].rad
				    + atom->rad
				    - overlap )
			       && ( distance < MIN_METAL_1_HBOND_LENGTH ) )
			  || ( ( atom->act == METAL_2 
				 || target[neighbors[k]].act == METAL_2 )
			       && ( distance < target[neighbors[k]].rad
				    + atom->rad
				    - overlap )
			       && ( distance < MIN_METAL_2_HBOND_LENGTH ) )
			  ))
		  {
		    printf ( "intra BUMP %s %s %s with %s %s %s: %f (%f), orig = %f \n",
			     atom->name,
			     atom->residue,
			     atom->residue_num,
			     target[neighbors[k]].name,
			     target[neighbors[k]].residue,
			     target[neighbors[k]].residue_num,
			     dist_fun ( atom->pos, target[neighbors[k]].pos ),
			     atom->rad + target[neighbors[k]].rad 
			     - overlap -
			     dist_fun ( atom->pos, target[neighbors[k]].pos ),
			     atom->neighbor_dist[k] );
		    if ( number_of_bumps >= MAX_SIDE_CHAIN_BUMPS )
		      /* the vector `bump_point' is already full, thus there
			 is no hope for this ligand ==> bump-check failed */
		      {
			global->number_of_bumps = number_of_bumps;
			return FAILURE;
		      }
		    if ( sum != ACCEPTOR_DONOR
			 && sum != ACCEPTOR_DONEPTOR
			 && sum != DONOR_DONEPTOR
			 && sum != DONEPTOR_DONEPTOR
			 && sum != ACCEPTOR_DONOR_HYDROGEN
			 && sum != DONEPTOR_DONOR_HYDROGEN
			 && atom->act != METAL_1
			 && atom->act != METAL_2
			 && target[neighbors[k]].act != METAL_1
			 && target[neighbors[k]].act != METAL_2
			 && distance < atom->rad 
			 + target[neighbors[k]].rad 
			 - overlap )
		      global->total_overlap +=
			atom->rad + target[neighbors[k]].rad 
			- overlap - distance;
		    if ( ( sum == ACCEPTOR_DONOR
			   || sum == DONOR_DONEPTOR
			   || sum == ACCEPTOR_DONEPTOR
			   || sum == DONEPTOR_DONEPTOR )
			 && distance < MIN_HBOND_LENGTH )
		      global->total_overlap +=
			MIN_HBOND_LENGTH - distance;
		    if ( ( sum == ACCEPTOR_DONOR_HYDROGEN  
			   || sum == DONEPTOR_DONOR_HYDROGEN ) 
			 && distance < MIN_HBOND_LENGTH_HYDROGEN ) 
		      {
			global->total_overlap +=
			  MIN_HBOND_LENGTH_HYDROGEN - distance;
		      }
		    if ( ( atom->act == METAL_1 
			   || target[neighbors[k]].act == METAL_1 )
			 && ( distance < target[neighbors[k]].rad
			      + atom->rad
			      - overlap )
			 && ( distance < MIN_METAL_1_HBOND_LENGTH ) )
		      {
			global->total_overlap +=
			  MIN_METAL_1_HBOND_LENGTH - distance;
		      }
		    if ( ( atom->act == METAL_2 
			   || target[neighbors[k]].act == METAL_2 )
			 && ( distance < target[neighbors[k]].rad
			      + atom->rad
			      - overlap )
			 && ( distance < MIN_METAL_2_HBOND_LENGTH ) )
		      {
			global->total_overlap +=
			  MIN_METAL_2_HBOND_LENGTH - distance;
		      }
		    
		    target_bumps[number_of_bumps] = 
		      global->target_residues[i].start_atom+j;
		    ligand_bumps[number_of_bumps] = (-1) * ( neighbors[k] + 1 );
		    number_of_bumps++;
		  }
		k++;
	      }
	  }
    }
  global->number_of_bumps = number_of_bumps;
  if ( number_of_bumps == 0 )
    return NO_BUMP;

/***** Added by PCS -- 22-Mar-00 *****/
  /* Added to be able to track bumps in the output and to be able to
     correlate with the output ligand mol2 and target pdb files */
  printf ("Above bumps for bumped match %d for ligand %s\n", 
	  global->number_of_bumped_ligands, global->ligand_file_name);

#ifdef BUMP_FILES
  sprintf ( filename, "%s/%s/%s/%s_ligands/bump_%s_ligand_%04d.mol2",
	    global->data_root, global->protein, global->template,
	    global->database, global->ligand_file_name, 
	    global->number_of_bumped_ligands );
  write_ligand_mol2 ( global, filename );
  sprintf ( filename, "%s/%s/%s/%s_targets/bump_%s_target_%04d.pdb",
	    global->data_root, global->protein, global->template,
	    global->database, global->ligand_file_name,
	    global->number_of_bumped_ligands);
  write_target_pdb ( global, filename );
#endif

  global->number_of_bumped_ligands ++;
/*************************************/
 return BUMP;
}
#endif

