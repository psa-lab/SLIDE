/******************************************************************************
 * Copyright (c) 2006,2007, Michigan State University (MSU) Board of Trustees.
 *   All rights reserved.
 *
 * This file is part of the ASCbase Software project.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Authors: Jeffrey Van Voorst, vanvoor4@msu.edu
 *          Leslie Kuhn, Ph.D., KuhnL@msu.edu 
 *****************************************************************************/
  
/*
 * $Source: /psa/share/repository/slide/src/slide/octree_array.c,v $
 * $Revision: 1.2 $
 * $Author: vanvoor4 $
 * $Date: 2009/03/09 20:43:54 $
 * 
 * $Log: octree_array.c,v $
 * Revision 1.2  2009/03/09 20:43:54  vanvoor4
 * checked in, in the event that someone later would find is useful.
 * Now that is updated, feel free to remove it -- CVS will still keep a copy
 * of course
 *
 * Revision 1.1  2009/02/25 20:08:21  vanvoor4
 * initial checkin
 *
 *
 * 
 */ 

#include <stdio.h>
#include <math.h>
#include <mymalloc.h>
#include <octree_array.h>

#define SLIDE_DEBUG

int line_in_sphere(const float *centroid, const float half_width,
                   const float *center, const float radius, const int idx,
                   const float *sign);

void
init_octree_node(octree_node_p node)
{
  node->parent = 0;
  memset(node->children, 0, 8 * sizeof(*node->children));
  memset(node->centroid, 0, 3 * sizeof(*node->centroid));
  node->half_width = 0;
  node->data = 0;
  node->num_el = 0;
  node->added_data = 0;
  node->num_added_el = 0;
  node->added_size = 0;
  node->leaf_node = 0;
}

void 
build_octree(octree_p tree, atom_pt atoms, float *positions, size_t num_atoms, 
             float min_half_width, size_t max_el_in_bin)
{
  float max_pt[3], min_pt[3];
  const float *p;
  size_t i;
  float longest_side;
  float tmp;
  float *pos = 0;
  octree_data_p data_ptr = 0;
  octree_data_p *root_data_ptr = 0;
  octree_node_p root_node = &tree->root_node;
  float *positions_end = &positions[3*num_atoms];
  atom_pt atom = 0;
  float cube_width;

  tree->min_half_width = min_half_width;
  tree->max_numel = max_el_in_bin;
  init_octree_node(root_node);
  

  /* Spin through the positions to find the min and max values */
  memcpy(max_pt, positions, 3*sizeof(*positions));
  memcpy(min_pt, positions, 3*sizeof(*positions));
  for(p = positions; p < positions_end; p += 3)
    for(i = 0; i < 3; ++i){
      if(max_pt[i] < p[i]) max_pt[i] = p[i];
      if(min_pt[i] > p[i]) min_pt[i] = p[i];
    }

  /*
   * Use a square at the cost of more storage (larger tree) and possibly
   * a longer search time.  We do want a reasonable tolerance on the added
   * to each side so that points do not easily move outside
   */

  /* Compute the centroid and halfwidths */
  longest_side = 0;
  for(i = 0; i < 3; ++i){
    root_node->centroid[i] = 0.5 * (max_pt[i] + min_pt[i]);
    tmp = max_pt[i] - min_pt[i];
    longest_side = (longest_side < tmp ? tmp : longest_side);
  }

  /*
   * Should be an input variable at some point in time -- a Lys or Arg
   * sidechain can easily move more than 10.0 (A)
   */
  longest_side += 40.0;

  /* Tune size to be divisible by min half width otherwise we get 
     a less tight set of boxes -- need to test on a number of different
     point sets in order to get an idea of its impact*/
  cube_width = min_half_width;
  for( ; longest_side > min_half_width; cube_width *= 2.0, longest_side *= 0.5);
  root_node->half_width = 0.5 * cube_width;

  /* Fill the data pointer array */
  tree->data = (octree_data_p) mymalloc(num_atoms * sizeof(octree_data_t));
  tree->data_end = &tree->data[num_atoms];
  atom = atoms;
  pos = positions;
  for(data_ptr = tree->data; data_ptr < tree->data_end; ++data_ptr, ++atom){
    data_ptr->position = pos;
    data_ptr->atom = atom;
    data_ptr->orig_node = 0;
    data_ptr->cnode = 0;
    pos += 3;
  }

  /* Fill the array for the head node */
  root_node->data = (octree_data_p*)mymalloc(num_atoms * sizeof(octree_data_p));
  root_node->num_el = num_atoms; 
  root_data_ptr = root_node->data;
  for(data_ptr = tree->data; data_ptr < tree->data_end; ++data_ptr){
    *root_data_ptr = data_ptr;
    ++root_data_ptr;
  }

  subdivide(tree, root_node, 0);
}

void
subdivide(octree_p tree, octree_node_p cnode, int level)
{
  size_t idx;
  size_t j;
  size_t dim_mask;
  int i;
  const float *centroid = cnode->centroid;
  octree_data_p *data;
  octree_data_p *data_end = &(cnode->data[cnode->num_el]);
  const float child_h_width = 0.5 * cnode->half_width;
  octree_node_p child = 0;
  int octant_count[8];
  int *partition;
  int *octant;

/*
  printf("subdivision level %d\n", level);
  printf("half width %f\n", cnode->half_width);
  printf("centroid: %f %f %f\n", cnode->centroid[0], cnode->centroid[1], 
         cnode->centroid[2]);
  printf("number of positions: %d\n", cnode->num_el);
*/

  if(cnode->num_el <= tree->max_numel || child_h_width < tree->min_half_width)
  {
    for(data = cnode->data; data < data_end; ++data){
      (*data)->orig_node = cnode;
      (*data)->cnode = cnode;
    }
    cnode->leaf_node = 1;
    return;
  }

  memset(octant_count, 0, 8 * sizeof(*octant_count));
  partition = (int*) mymalloc(cnode->num_el * sizeof(int));
  octant = partition;

  /* partition based on the position of each point */
  for(data = cnode->data; data < data_end; ++data, ++octant){
    idx = 0;
    for(i = 2; i > -1; --i){
      idx <<= 1;
      if((*data)->position[i] >= centroid[i]) ++idx;
    }
    *octant = idx;
    ++octant_count[idx];
  }

  /* Allocate the children for the octants with data and the childrens' 
   * data arrays */
  for(idx = 0; idx < 8; ++idx){
    if(octant_count[idx] == 0) continue;

    child = (octree_node_p) mymalloc(sizeof(octree_node_t));
    init_octree_node(child);
    child->parent = cnode;
    child->half_width = child_h_width;
    dim_mask = 1; 
    for(j = 0; j < 3; ++j, dim_mask <<= 1){
      if(idx & dim_mask) child->centroid[j] = centroid[j] + child_h_width;
      else child->centroid[j] = centroid[j] - child_h_width;
    }
    child->data = 
      (octree_data_p*) mymalloc(octant_count[idx] * sizeof(octree_data_p));
    cnode->children[idx] = child;
  }

  /* Partition the data based on the partion array */
  octant = partition;
  for(data = cnode->data; data < data_end; ++data, ++octant){
    child = *(cnode->children + *octant);
    child->data[child->num_el] = *data;
    ++child->num_el;
  }
  free(partition);
  partition = 0;
  octant = 0;
  free(cnode->data);
  cnode->data = 0;
  cnode->num_el = 0;

  /* Recursively subdivide children */
  for(j = 0; j < 8; ++j)
    if(octant_count[j]) subdivide(tree, cnode->children[j], level + 1);
}

void free_octree(octree_p tree)
{
  free_octree_node(&tree->root_node);
}

void free_octree_node(octree_node_p cnode)
{
  int leaf_node = 0;
  int i;

  if(cnode->num_el){
    if(cnode->data) free(cnode->data);
    cnode->data = 0;
    cnode->num_el = 0;
    leaf_node = 1;
  }
  if(cnode->num_added_el){
    if(cnode->added_data) free(cnode->added_data);
    cnode->added_data = 0;
    cnode->num_added_el = 0;
    cnode->added_size = 0;
    leaf_node = 1;
  }
  if(leaf_node) return;

  for(i = 0; i < 8; ++i)
    if(cnode->children[i]){
      free_octree_node(cnode->children[i]);
      free(cnode->children[i]);
      cnode->children[i] = 0;
    }
}

void 
find_close_atoms(const octree_node_p cnode, const float *pt, 
                 const float pt_tol, atom_pt **close_atoms_array_ptr, 
                 size_t *num_atoms, int level)
{
  int i;
  size_t new_num_atoms;
  octree_data_p *data_p;
  atom_pt *atom_p;
  int rv;

  /*printf("\nposition to find: %f %f %f\n", pt[0], pt[1], pt[2]);*/
  /*
  printf("\nlevel: %d\n", level + 1);
  printf("size: %d\n", cnode->num_el);
  printf("h_width: %f\n", cnode->half_width);
  printf("centroid: %f %f %f\n", cnode->centroid[0], cnode->centroid[1], 
         cnode->centroid[2]);
         */
  /* 
   * We could check distance and only push back those meeting the 
   * distance criterion here but that might cause the distance 
   * calculations for the points close to pt to be computed at least twice
   */
  if(cnode->num_el + cnode->num_added_el){
    new_num_atoms = *num_atoms + cnode->num_el + cnode->num_added_el;
    *close_atoms_array_ptr = 
      (atom_pt*) myrealloc(*close_atoms_array_ptr, 
                           new_num_atoms * sizeof(atom_pt));

    if(cnode->num_el){
      atom_p = &((*close_atoms_array_ptr)[*num_atoms]);
      data_p = cnode->data;
      for(; data_p < &(cnode->data[cnode->num_el]); ++data_p, ++atom_p)
        *atom_p = (*data_p)->atom; 
      *num_atoms += cnode->num_el;
    }
    if(cnode->num_added_el){
      atom_p = &((*close_atoms_array_ptr)[*num_atoms]);
      data_p = cnode->added_data;
      for(;data_p < &(cnode->added_data[cnode->num_added_el]); ++data_p){
        *atom_p = (*data_p)->atom; 
	++atom_p;
      }
      *num_atoms += cnode->num_added_el;
    }
    return;
  }

  for(i = 0; i < 8; ++i){
    if(cnode->children[i] == 0) continue;

    rv = pt_in_cube(cnode->children[i]->centroid, 
                    cnode->children[i]->half_width, pt, pt_tol);
    if(rv == CENTER_IN_CUBE ||
       (rv == POSSIBLE_INTERSECTION && 
        cube_and_sphere_intersect(cnode->children[i]->centroid, 
                                  cnode->children[i]->half_width, pt, pt_tol))){
      find_close_atoms(cnode->children[i], pt, pt_tol, close_atoms_array_ptr, 
                       num_atoms, level + 1);
    }
  }
}

void 
find_close_points(const octree_node_p cnode, const float *pt, 
                 const float pt_tol, float ***close_points_p, size_t *npts, 
                 int level)
{
  int i;
  size_t new_npts;
  octree_data_p *data_p;
  float **pos_p;
  int rv;

  /*printf("\nposition to find: %f %f %f\n", pt[0], pt[1], pt[2]);*/
  /*
  printf("\nlevel: %d\n", level + 1);
  printf("size: %d\n", cnode->num_el);
  printf("h_width: %f\n", cnode->half_width);
  printf("cendroid: %f %f %f\n", cnode->centroid[0], cnode->centroid[1], 
         cnode->centroid[2]);
         */

  if(cnode->num_el + cnode->num_added_el){
    new_npts = *npts + cnode->num_el + cnode->num_added_el;
    *close_points_p = 
      (float**) myrealloc(*close_points_p, new_npts * sizeof(float*));

    if(cnode->num_el){
      pos_p = &((*close_points_p)[*npts]);
      data_p = cnode->data;
      for(; data_p < &(cnode->data[cnode->num_el]); ++data_p, ++pos_p)
        *pos_p = (*data_p)->position;
      *npts += cnode->num_el;
    }
    if(cnode->num_added_el){
      pos_p = &((*close_points_p)[*npts]);
      data_p = cnode->added_data;
      for(;data_p < &(cnode->added_data[cnode->num_added_el]); ++data_p){
        *pos_p = (*data_p)->position;
        ++pos_p;
      }
      *npts += cnode->num_added_el;
    }
    return;
  }

  for(i = 0; i < 8; ++i){
    if(cnode->children[i] == 0) continue;

    rv = pt_in_cube(cnode->children[i]->centroid, 
                    cnode->children[i]->half_width, pt, pt_tol);
    if(rv == CENTER_IN_CUBE ||
       (rv == POSSIBLE_INTERSECTION && 
        cube_and_sphere_intersect(cnode->children[i]->centroid, 
                                  cnode->children[i]->half_width, pt, pt_tol))){
      find_close_points(cnode->children[i], pt, pt_tol, close_points_p, npts, 
                        level + 1);
    }
  }
}


int
point_in_cube(const float *centroid, const float half_width, const float *pt)
{
  float dist[3];
  int i;

  for(i = 0; i < 3; ++i){
    dist[i] = centroid[i] - pt[i];     
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }

  /* If pt in box, return true */
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width) 
    return 1;

  return 0;
}

int 
cube_and_sphere_intersect(const float *centroid, const float half_width,
                          const float *center, const float radius)
{
  float dist[3];
  float tmp;
  float sign[3];
  int i;

  sign[0] = sign[1] = sign[2] = 0;

/*  
  printf("\nCentroid: %f %f %f\n", centroid[0], centroid[1], centroid[2]);
  printf("half_width: %f\n", half_width);
  printf("center: %f %f %f\n", center[0], center[1], center[2]);
  printf("radius: %f\n", radius);
  */

  for(i = 0; i < 3; ++i){
    tmp = centroid[i] - half_width;
    if(tmp - radius <= center[i] && center[i] <= tmp) sign[i] = -1;
    
    tmp = centroid[i] + half_width;
    if(tmp <= center[i] && center[i] <= tmp + radius) sign[i] = 1;
  }

  /*
  printf("Sign: %2f %2f %2f\n", sign[0], sign[1], sign[2]);
  */

  /* Center is in the rectangular solid defined by a face and ?? and is closer 
   * than radius to that face? */
  if((sign[0] && !sign[1] && !sign[2]) || (!sign[0] && sign[1] && !sign[2]) ||
     (!sign[0] && !sign[1] && sign[2])){
    /*
    printf("Center of sphere is adjacent to a face\n");
    */
    return CUBE_AND_SPHERE_INTERSECT;
  }

  /* Sphere contains a cube corner? */
  if(sign[0] && sign[1] && sign[2]){
    for(i = 0; i < 3; ++i) 
      dist[i] = centroid[i] + sign[i] * half_width - center[i];
    /*
    printf("dist[]: %f %f %f\n", dist[0], dist[1], dist[2]);
    */
    if(radius*radius >= dist[0]*dist[0] + dist[1]*dist[1] + dist[2]*dist[2]){
      /*
      printf("Sphere contains a cube corner\n");
      */
      return CUBE_AND_SPHERE_INTERSECT;
    }else{
      /*
      printf("Sphere does not contain a cube corner\n");
      */
      return NO_INTERSECTION;
    }
  }

  if(!sign[0]){
    if(line_in_sphere(centroid, half_width, center, radius, 0, sign))
      return CUBE_AND_SPHERE_INTERSECT;
  }else if(!sign[1]){
    if(line_in_sphere(centroid, half_width, center, radius, 1, sign))
      return CUBE_AND_SPHERE_INTERSECT;
  }else if(!sign[2]){
    if(line_in_sphere(centroid, half_width, center, radius, 2, sign))
      return CUBE_AND_SPHERE_INTERSECT;
  }

  /*printf("No intersection!\n");*/
  return NO_INTERSECTION;
}

void
print_bin_trace(const octree_p otree, const octree_node_p cnode, int level)
{
  size_t i;
  if(cnode->num_el){
    printf("\nBin at level %d", level);
    printf("\nCentroid: %f %f %f", cnode->centroid[0], cnode->centroid[1], 
           cnode->centroid[2]);
    printf("\nHalf Width: %f", cnode->half_width);
    printf("\nPositions:\n");
    
    for(i = 0; i < cnode->num_el; ++i){
      printf("\t%4d: %f %f %f\n", cnode->data[i] - otree->data, 
             cnode->data[i]->position[0],
             cnode->data[i]->position[1], cnode->data[i]->position[2]);
    }
    return;
  }

  for(i = 0; i < 8; ++i)
    if(cnode->children[i]) 
      print_bin_trace(otree, cnode->children[i], level + 1);
}

int
line_in_sphere(const float *centroid, const float half_width,
               const float *center, const float radius, const int idx,
               const float *sign)
{
  int i;
  float pt[3];
  float RHS = radius*radius;
  float S1, S2;

  for(i = 0; i < 3; ++i)
    if(i != idx){
      pt[i] = centroid[i] + sign[i]*half_width;
      RHS -= (pt[i] - center[i]) * (pt[i] - center[i]);
    }
  /*
  printf("pt: %f %f %f\n", pt[0], pt[1], pt[2]);
  */

  if(RHS < 0) return 0;
  RHS = sqrt(RHS);
  /*
  printf("RHS: %f\n", RHS);
  */

  S1 = center[idx] - RHS;
  S2 = center[idx] + RHS;
  if((S1 < centroid[idx] - half_width && S2 < centroid[idx] - half_width) ||
     (S1 > centroid[idx] + half_width && S2 > centroid[idx] + half_width))
    return 0;
  return 1;
}

void check_position(octree_p tree, const atom_pt atom)
{
  int i;
  float dist[3];
  octree_data_p data_ptr = tree->data + (atom - tree->data->atom);
  int rv = 0;
  octree_node_p cnode = data_ptr->cnode;

  rv = pt_in_cube(cnode->centroid, cnode->half_width, data_ptr->position, 4.5);
  if(rv == CENTER_IN_CUBE) return;

  fprintf(stderr, "\nMISSED ME!");
  fprintf(stderr, "\nCentroid: %f %f %f", cnode->centroid[0], 
          cnode->centroid[1], cnode->centroid[2]);
  fprintf(stderr, "\nHalf Width: %f", cnode->half_width);
  fprintf(stderr, "\natom position: %f %f %f\n", atom->pos[0], atom->pos[1],
          atom->pos[2]);

#if 0
  for(i = 0; i < 3; ++i){
    dist[i] = cnode->centroid[i] - atom->pos[i];
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }
  fprintf(stderr, "\ndistances: %f %f %f\n", dist

  /* If center in box, return true */
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width)
#endif


}

#if 0
void 
update_bin(octree_p tree, const atom_pt atom)
{
  octree_data_p data_ptr = tree->data + (atom - tree->data.atom);
  _update_bin(data_ptr->cnode, atom)
}
#endif

void
_update_bin(octree_node_p cnode, octree_data_p atom_data)
{
  octree_node_p child;
  float child_h_width;
  int dim_mask;
  int i;
  int idx;
  int rv = 0;
  const float my_tol = 4.5;
  atom_pt atom = atom_data->atom;
  rv = pt_in_cube(cnode->centroid, cnode->half_width, atom->pos, my_tol);
#if 0

  // here for debugging
  if(rv == CENTER_IN_CUBE && cnode->leaf_node && cnode == atom_data->cnode) return;

  fprintf(stderr, "\nCentroid: %f %f %f", cnode->centroid[0], 
          cnode->centroid[1], cnode->centroid[2]);
  fprintf(stderr, "\nHalf Width: %f", cnode->half_width);
#endif

  // point is in current bin
  if(rv == CENTER_IN_CUBE){
    // Is current bin a leaf node?
    if(cnode->leaf_node){
      // Point did not move out of its current bin
      if(cnode == atom_data->cnode) return;
      // Point has moved back to its original bin
      else if(cnode == atom_data->orig_node){
        atom_data->cnode = atom_data->orig_node;
        return; 
      // Point is in a different bin than previously
      }else{
        if(cnode->added_size == 0 || cnode->added_size >= cnode->num_added_el){
          cnode->added_size += 10;
          cnode->added_data = 
            (octree_data_p*) myrealloc(cnode->added_data, cnode->added_size * 
                                       sizeof(cnode->added_data[0]));
        } 
        atom_data->cnode = cnode;
        cnode->added_data[cnode->num_added_el] = atom_data;
        ++cnode->num_added_el;
        return;
      }
    // not a leaf node -- go down the tree
    }else{
      
      idx = 0;
      for(i = 2; i > -1; --i){
        idx <<= 1;
        if(atom->pos[i] >= cnode->centroid[i]) ++idx;
      }
      if(cnode->children[idx]) _update_bin(cnode->children[idx], atom_data);
      // Need to create a bin
      else{
        child = (octree_node_p) mymalloc(sizeof(octree_node_t));
        init_octree_node(child);
        child->parent = cnode;
        child_h_width = cnode->half_width / 2.0;
        dim_mask = 1; 
        for(i = 0; i < 3; ++i, dim_mask <<= 1){
          if(idx & dim_mask) 
            child->centroid[i] = cnode->centroid[i] + child_h_width;
          else child->centroid[i] = cnode->centroid[i] - child_h_width;
        }
        child->added_size = 10;
        child->added_data = 
          (octree_data_p*) mymalloc( 10 * sizeof(child->added_data[0]));
        atom_data->cnode = cnode;
        child->added_data[0] = atom_data;
        child->num_added_el = 1;
        child->leaf_node = 1;
        child->half_width = child_h_width;
        cnode->children[idx] = child;
      } 
    } 
  // Point is not in the current bin -- go up the tree
  }else{
    if(cnode->parent) _update_bin(cnode->parent, atom_data);
    else
      fprintf(stderr, "TREE is too small -- need to increase the initial "
              "volume of consideration\n");
  }
}

void
reset_positions(octree_p tree)
{
  octree_data_p ptr;
  for(ptr = tree->data; ptr < tree->data_end; ++ptr)
    ptr->cnode = ptr->orig_node;

  /* this is not necessary for production, but can be very helpful when 
   * debugging */
  zero_out_moved_positions(&tree->root_node);
}

void 
zero_out_moved_positions(octree_node_p cnode)
{
  int i;
  cnode->num_added_el = 0;

#ifdef SLIDE_DEBUG
  if(cnode->added_data){
    memset(cnode->added_data, 0, 
           cnode->added_size * sizeof(cnode->added_data[0]));
    return;
  }
#endif

  for(i = 0; i < 8; ++i)
    if(cnode->children[i]) zero_out_moved_positions(cnode->children[i]); 
}

