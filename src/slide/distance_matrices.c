/*
 *  $Source: /psa/share/repository/slide/src/slide/distance_matrices.c,v $
 *  $Revision: 1.2 $
 *  $Author: vanvoor4 $
 *  $Date: 2009/03/09 20:25:36 $
 *  
 *  $Log: distance_matrices.c,v $
 *  Revision 1.2  2009/03/09 20:25:36  vanvoor4
 *  changed from octree to distance_array
 *
 *  Revision 1.1  2009/02/25 20:08:21  vanvoor4
 *  initial checkin
 *
 *
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <mymalloc.h>
#include <err_handle.h>
#include <dist_fun.h>
#include <distance_matrices.h>

static const float MAX_DIST = 8.0;
static const float MAX_SQRD_DIST = 64.0;
static const float SQRD_NBR_DIST = 64.0;

/*
void fill_row(const octree_node_p cnode, const octree_data_p data_el,
              const octree_data_p data_start, int nbr_check, float *row);
*/

int non_covalently_bound2(atom_pt atom1, atom_pt atom2);

void
init_target_nbr_arrays(atom_pt target_atoms, const int num_atoms)
{
  int i, j; 
  float sq_dist;

  for(i = 0; i < num_atoms; ++i){
    for(j = 0; j < i; ++j){
      sq_dist = squared_dist(target_atoms[i].pos, target_atoms[j].pos);
      if(sq_dist <= SQRD_NBR_DIST && 
         non_covalently_bound2(&target_atoms[i], &target_atoms[j])){
        target_atoms[i].neighbors[target_atoms[i].num_nbrs] = &target_atoms[j];
        target_atoms[j].neighbors[target_atoms[j].num_nbrs] = &target_atoms[i];
        target_atoms[i].neighbor_dist[target_atoms[i].num_nbrs] = 
          target_atoms[j].neighbor_dist[target_atoms[j].num_nbrs] = 
          sqrt(sq_dist);
        ++(target_atoms[i].num_nbrs);
        ++(target_atoms[j].num_nbrs);
      }
    }
  }
}

int non_covalently_bound2(atom_pt atom1, atom_pt atom2)
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

void
initialize_inter_dist_matrix(const float* targ_positions, 
                             const size_t num_targ_atoms,
                             const float* lig_positions, 
                             const size_t num_lig_atoms,
                             float *M, const float max_dist)
{
  size_t i, j;
  float *m = M;
  float sq_dist;
  const float* l_pos;
  const float* t_pos;
  const float max_sqrd_dist = max_dist * max_dist;
  const float* targ_pos_end = targ_positions + 3*num_targ_atoms;
  const float* lig_pos_end = lig_positions + 3*num_lig_atoms;

  for(l_pos = lig_positions; l_pos < lig_pos_end; l_pos += 3)
    for(t_pos = targ_positions; t_pos < targ_pos_end; t_pos += 3, ++m){
      sq_dist = squared_dist(l_pos, t_pos);
      *m = (sq_dist <= max_sqrd_dist ? sq_dist : FLT_MAX);
    }
}

#if 0
void
initialize_intra_dist_matrix(octree_p octree, float *M, int nbr_check)
{
  size_t i,j;
  size_t num_el = octree->data_end - octree->data; 
  float *m;

  /* Fill lower triangle + diagonal */
  for(i = 0; i < num_el; ++i){
    m = &M[i*num_el];
    for(j = 0; j < i; ++j, ++m) *m = FLT_MAX;
    fill_row(&octree->root_node, &octree->data[i], octree->data, nbr_check,
             &M[i*num_el]);
    *m = 0.0;
  }

  /* Make this a square matrix by filling upper triangle */
  for(i = 0; i < num_el; ++i)
    for(j = 0; j < i; ++j) M[j*num_el + i] = M[i*num_el + j];
}

/* This method should be made as efficient as possible since it must be called
 * for each docking.
 */
void
initialize_inter_dist_matrix(const octree_p octree, 
                             const float* octree_positions,
                             const float *positions, const size_t num_pos, 
                             float *M, const size_t stride, 
                             const float max_dist)
{
  size_t npts;
  float *r;
  float sq_dist;
  const float *pos;
  float **close_pts = 0;
  float **close_pt = 0;
  float *row = M;
  const float max_sqrd_dist = max_dist * max_dist;

  for(r = M; r < &M[num_pos * stride]; ++r) *r = FLT_MAX;

  pos = positions;
  for( ; pos < &positions[3*num_pos]; pos += 3, row += stride){
    get_close_points(octree, pos, max_dist, &close_pts, &npts);
    for(close_pt = close_pts ; close_pt < &close_pts[npts]; ++close_pt){
      sq_dist = squared_dist(pos, *close_pt); 
      if(sq_dist <= max_sqrd_dist)
        row[(*close_pt - octree_positions)/3] = sq_dist;
    }
    if(close_pts) free(close_pts);
    close_pts = 0;
    close_pt = 0;
    npts = 0;
  }
}
#endif


#if 0
void
fill_row(const octree_node_p cnode, const octree_data_p data_el, 
         const octree_data_p data_start, int nbr_check, float *row)
{
  octree_data_p *nbr_p;
  octree_data_p nbr;
  float sq_dist;
  size_t i;
  /*
  size_t nbr_idx;
  size_t cur_idx = data_el - data_start;
  */
  int rv;

  if(cnode->num_el){
    nbr_p = cnode->data;
    for(; nbr_p < &(cnode->data[cnode->num_el]); ++nbr_p){
      nbr = *nbr_p;
      sq_dist = squared_dist(data_el->position, nbr->position);
      /*nbr_idx = nbr - data_start;*/
      /*printf("(%d,%d) = %f \n", cur_idx, nbr_idx, sq_dist);*/

      if(*nbr_p >= data_el) continue;
      nbr = *nbr_p;

      sq_dist = squared_dist(data_el->position, nbr->position);
      /*nbr_idx = nbr - data_start;*/
      /*
      if(sq_dist <= MAX_SQRD_DIST) row[nbr_idx] = sq_dist;
      */

      /* Check for neighbors of the current atom */
      if(nbr_check && sq_dist <= SQRD_NBR_DIST &&
         non_covalently_bound2(data_el->atom, nbr->atom)){
        data_el->atom->neighbors[data_el->atom->num_nbrs] = nbr->atom;
        nbr->atom->neighbors[nbr->atom->num_nbrs] = data_el->atom;
        /* Use distance here (rather than squared distance) since this is 
         * computed once per screening run
         */
        data_el->atom->neighbor_dist[data_el->atom->num_nbrs] = 
          nbr->atom->neighbor_dist[nbr->atom->num_nbrs] = 
          sqrt(sq_dist);

        data_el->atom->num_nbrs++;
        nbr->atom->num_nbrs++;
        if(nbr->atom->num_nbrs >= MAX_NEIGHBOR_ATOMS ||
           data_el->atom->num_nbrs >= MAX_NEIGHBOR_ATOMS)
          err_panic2("fill_row", 
                     "An atom has more than MAX_NEIGHBOR_ATOMS neighbors");
      }
    }
    return;
  }

  for(i = 0; i < 8; ++i){
    if(cnode->children[i] == 0) continue;

    rv = pt_in_cube(cnode->children[i]->centroid, 
                    cnode->children[i]->half_width, data_el->position, 
                    MAX_DIST);
    if(rv == CENTER_IN_CUBE ||
       (rv == POSSIBLE_INTERSECTION && 
        cube_and_sphere_intersect(cnode->children[i]->centroid, 
                                  cnode->children[i]->half_width, 
                                  data_el->position, MAX_DIST)))
      fill_row(cnode->children[i], data_el, data_start, nbr_check, row);
  }
}
#endif

