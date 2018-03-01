/*
 *  $Source: /psa/share/repository/slide/src/slide/inc/distance_matrices.h,v $
 *  $Revision: 1.2 $
 *  $Author: vanvoor4 $
 *  $Date: 2009/03/09 20:25:34 $
 *  
 *  $Log: distance_matrices.h,v $
 *  Revision 1.2  2009/03/09 20:25:34  vanvoor4
 *  changed from octree to distance_array
 *
 *  Revision 1.1  2009/02/25 19:55:50  vanvoor4
 *  Removed the extern from the front of the function declarations.
 *  Ideally it shouldn't make a difference, but we have found with gcc that
 *  having the extern in front of the declarations allows one to look like they
 *  are using the function without the need to include the appropriate header
 *  file.  Unfortunately, the results are then undefined and to top it off gcc does
 *  not give any warnings.
 *
 *  some of the more simple functions were retooled so that passing in global
 *  is not required.  This is getting to be a big mess with sticking everything
 *  in global and passing it to all the functions.
 *
 *
 */

#ifndef DISTANCE_MATRICES_HEADER_FILE_INCLUDED
#define DISTANCE_MATRICES_HEADER_FILE_INCLUDED

#include <octree_array.h>

static inline int diff_int(int  a, int  b)
{
  return (a > b ? a - b : b - a);
}

void
init_target_nbr_arrays(atom_pt target_atoms, const int num_atoms);

void
initialize_inter_dist_matrix(const float* targ_positions, 
                             const size_t num_targ_atoms, 
                             const float* lig_positions, 
                             const size_t num_lig_atoms,
                             float *M, const float max_dist);

#endif
