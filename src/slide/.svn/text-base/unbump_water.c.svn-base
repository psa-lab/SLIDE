#include "types.h"
#include "defs.h"
#include "basics.h"
#include "dist_fun.h"
#include "unbump_translate.h"

/*
 *  The name of this function might be confusing: it computes the minimal
 *  translation along the global vector that is needed to get rid of
 *  all bumps, and so the maximal factor of all bumps has to be found.
 */
int  compute_minimal_water_translation ( vector_pt vectors,
					 float     *radii,
					 int       number_of_bumps,
					 vector_t  global_vector,
					 float     overlap,
					 double    *maximal_translation )
{
  int    i;
  double translation;
  
  *maximal_translation = 0.0;
  for ( i = 0; i < number_of_bumps; i++ )
    {
      if ( compute_translation_coefficient ( vectors[i],
					     global_vector,
					     (double) ( radii[i]
							+ WATER_RAD
							- overlap
							+ 0.00001 ),
					     &translation ) == FAILURE )
	return FAILURE;
      if ( compare_double(translation, *maximal_translation) == 1 )
	  *maximal_translation = translation;
    }
  return SUCCESS;
}

int  check_atom_for_bump ( atom_pt   water,
			   atom_pt   atom,
			   vector_pt vectors,
			   float     *radii,
			   int       number_of_bumps,
			   float     overlap )
{
  float  distance;
  int    i;

  distance = dist_fun ( water->pos,
			atom->pos );
  if ( ( atom->act != NOTHING 
	 && distance < MIN_HBOND_LENGTH )
       || ( atom->act == NOTHING
	    && distance < WATER_RAD + atom->rad + overlap ) )
    {
      for ( i = 0; i < 3; i++ )
	vectors[number_of_bumps][i] = water->pos[i] - atom->pos[i];
      radii[number_of_bumps] = atom->rad;
      number_of_bumps++;
    }
  return number_of_bumps;
}

int  unbump_water ( global_data_pt  global,
		    atom_pt         water )
{
  atom_pt   atom;
  vector_t  vectors[MAX_NEIGHBOR_ATOMS];
  vector_t  global_vector;
  float     radii[MAX_NEIGHBOR_ATOMS];
  double    factor;
  int       number_of_bumps,
            number_of_iterations;
  int       i;

  number_of_iterations = 0;
  number_of_bumps = 0;
  do 
    {
      if ( number_of_iterations > 0 )
	/* in the first iteration we do not know the bumps at this point,
	   so skip this part */
	{
	  /* compute the average translation direction to resolve the bumps */
	  compute_global_vector ( vectors,
				  number_of_bumps,
				  global_vector );	  
	  /* compute the actual amount to translate the water into the
	     average direction to get rid of all bumps at once */
	  if ( compute_minimal_water_translation ( vectors,
						   radii,
						   number_of_bumps,
						   global_vector,
						   global->side_chain_overlap,
						   &factor ) == FAILURE )
	    return FAILURE;
	  /* actually translate the water molecule */
	  for ( i = 0; i < 3; i++ )
	    water->pos[i] += ( global_vector[i] * factor );
	  global->number_of_water_translations++;
	}
      number_of_bumps = 0;
      /* find all atoms in the ligand that collide with this water molecule */
      for ( i = 0; i < global->ligand->number_of_atoms; i++ )
	{
	  atom = &global->ligand->atoms[i];
	  if ( atom->orbit != ADDED )
	    /* added hydrogens don't have a position (yet) */
	    number_of_bumps = 
	      check_atom_for_bump ( water,
				    atom,
				    vectors,
				    radii,
				    number_of_bumps,
				    global->side_chain_overlap );
	}
      /* find all atoms in the target that collide with this water molecule */
      i = 0; 
      while ( water->neighbors[i] != NO_MORE )
	{
	  atom = &global->target_atoms[water->neighbors[i]]; 
	  number_of_bumps = 
	    check_atom_for_bump ( water,
				  atom,
				  vectors,
				  radii,
				  number_of_bumps,
				  global->side_chain_overlap );
	  i++;
	}
      /* find all waters that collide with this water molecule */
      for ( i = 0; i < global->number_of_waters; i++ )
	{
	  atom = &global->waters[i];
	  if ( atom != water )
	    /* of course, don't check for collisions with the water itself */
	    number_of_bumps = 
	      check_atom_for_bump ( water,
				    atom,
				    vectors,
				    radii,
				    number_of_bumps,
				    global->side_chain_overlap );
	}
      number_of_iterations++;
    }
  while ( number_of_bumps > 0 
	  && number_of_iterations < MAX_NUMBER_WATER_UNBUMP_ITERATIONS );
  if ( number_of_bumps > 0 )
    return SUCCESS;
  return FAILURE;
}
