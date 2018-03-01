#include <string.h>
#include <math.h>
#include <types.h>
#include <mymalloc.h>
#include <err_handle.h>
#include <dist_fun.h>

#if 0
#include <initialize.h>
#include <read_pdb.h>
#include <read_mol2.h>

int load_ligand(char *filename, global_data_pt global);

void set_global_junk(global_data_pt global);
#endif

void
distance_array(dist_array_pt darray, const atom_pt atoms, int num_atoms, 
               const float box_width, const float max_dist)
{
  int i, j, k;
  float min_pt[3], max_pt[3];
  atom_pt a;
  const atom_pt atoms_end = atoms + num_atoms;
  dist_bin_pt bin_p;

  darray->box_width = box_width;

  /* Spin through the positions to find the min and max values */
  memcpy(min_pt, atoms->pos, 3*sizeof(min_pt[0]));
  memcpy(max_pt, atoms->pos, 3*sizeof(min_pt[0]));
  for(a = atoms; a < atoms_end; ++a){
    for(i = 0; i < 3; ++i){
      if(max_pt[i] < a->pos[i]) max_pt[i] = a->pos[i];
      if(a->pos[i] < min_pt[i]) min_pt[i] = a->pos[i];
    }
  }

  /* Make a box that has a buffer as well as side lengths divisible by 
   * box_width */
  for(i = 0; i < 3; ++i){
    max_pt[i] = box_width * ceil(max_pt[i] / box_width) + 2*box_width;
    min_pt[i] = box_width * floor(min_pt[i] / box_width) - 2*box_width;
    darray->num_bins[i] = (int) ((max_pt[i] - min_pt[i]) / box_width);
  }
  memcpy(darray->min_corner, min_pt, 3 * sizeof(min_pt[0]));
  memcpy(darray->max_corner, max_pt, 3 * sizeof(max_pt[0]));
 
  /* Populate the arrays */ 
  darray->bins = 
    (dist_bin_pt) mymalloc(darray->num_bins[0] * darray->num_bins[1] * 
                           darray->num_bins[2] * sizeof(dist_bin_t));
  bin_p = darray->bins;
  min_pt[0] = darray->min_corner[0];
  max_pt[0] = darray->min_corner[0] + box_width;
  for(i = 0; i < darray->num_bins[0]; ++i){
    min_pt[1] = darray->min_corner[1];
    max_pt[1] = darray->min_corner[1] + box_width;
    for(j = 0; j < darray->num_bins[1]; ++j){
      min_pt[2] = darray->min_corner[2];
      max_pt[2] = darray->min_corner[2] + box_width;
      for(k = 0; k < darray->num_bins[2]; ++k, ++bin_p){
        /* Intialize the bin variables */
        bin_p->atoms = 0;
        bin_p->num_atoms = 0;
        bin_p->size = 0;
        /* Populate the bin */
        for(a = atoms; a < atoms_end; ++a){
          if(min_pt[0] - max_dist <= a->pos[0] &&
             a->pos[0] <= max_pt[0] + max_dist &&
             min_pt[1] - max_dist <= a->pos[1] &&
             a->pos[1] <= max_pt[1] + max_dist && 
             min_pt[2] - max_dist <= a->pos[2] &&
             a->pos[2] <= max_pt[2] + max_dist){
            if(bin_p->size <= bin_p->num_atoms){
              bin_p->size = (bin_p->size ? 2*bin_p->size : 10);
              bin_p->atoms = (atom_pt *) myrealloc(bin_p->atoms, bin_p->size * 
                                                   sizeof(atom_pt));
            }
            bin_p->atoms[bin_p->num_atoms] = a;
            ++(bin_p->num_atoms);
          }
        }
        min_pt[2] += box_width;
        max_pt[2] += box_width;
      }
      min_pt[1] += box_width; 
      max_pt[1] += box_width; 
    }
    min_pt[0] += box_width;
    max_pt[0] += box_width;
  }
}

/*
 Need to fix this as it causes a segfault -- a bit confusing since my_free
 doesn't attempt to free null pointers.
*/
void
free_distance_array(dist_array_pt d)
{
  int i;
  int num_bins = d->num_bins[0] * d->num_bins[1] * d->num_bins[2];
  for(i = 0; i < num_bins; ++i)
    my_free(d->bins[i].atoms);
  my_free(d->bins);
}

dist_bin_pt
get_bin(dist_array_pt darray, float *pt)
{
  int i;
  int idx[3];

  /* Check if point falls outside the range of the bins */
  for(i = 0; i < 3; ++i)
    if(pt[i] < darray->min_corner[i] || darray->max_corner[i] < pt[i]) return 0;

  /* Compute the index */
  for(i = 0; i < 3; ++i)
    idx[i] = (int) floor((pt[i] - darray->min_corner[i]) / darray->box_width);
 
  return darray->bins + (idx[0] * darray->num_bins[1]*darray->num_bins[2] +
                         idx[1] * darray->num_bins[2] + idx[2]); 
}


/* Code to test distance_array */
#if 0
int
main()
{
  int rv;
  int i;
  int j;
  float sq_dist;
  global_data_pt global;

  global = initialize_global_data_structure();
  read_pdb("/home/vanvoor4/rognan_100_slide_test/data/1aaq/unbiased/in/1aaq.rad",
           global->target_atoms, global->target_residues, ALSO_HETERO,
           &global->number_of_target_atoms, &global->number_of_target_residues);
  set_global_junk(global);

  rv = load_ligand("/home/vanvoor4/rognan_100_slide_test/1aaq_ligand/1aaq.mol2", global);
  if(rv != SUCCESS){
    fprintf(stderr, "Reading of the ligand failed\n");
    exit(-1);
  }

  
  dist_array_t dist_array;
  distance_array(&dist_array, global->target_atoms, global->number_of_target_atoms, 4.5, 5.0);

//const dist_bin_pt
//get_bin(dist_array_pt darray, float *pt)
  dist_bin_pt bin = 0;
  for(i = 0; i < global->ligand->number_of_atoms; ++i){
    bin = get_bin(&dist_array, global->ligand->atoms[i].pos);
    for(j = 0; j < bin->num_atoms; ++j){
      sq_dist = squared_dist(global->ligand->atoms[i].pos, bin->atoms[j]->pos);
      if(sq_dist < 5.0*5.0){
        printf("%d\n", bin->atoms[j]->atom_number);
      }
    } 
  }
   

//  free(&dist_array);
}

int load_ligand(char *filename, global_data_pt global)
{
  FILE *MOL2;
  char err_msg[FILENAME_MAX];
  int rv = FATAL_FAILURE;
  
  if((MOL2 = open_mol2(filename)) == NULL) return FATAL_FAILURE;
  if((rv = read_mol2(MOL2, filename, global, NULL)) != SUCCESS){
    fclose(MOL2);
    return rv;
  }
  fclose(MOL2);
#if 0
  /* -- analyze_ligand() depends on interactions loaded from a pts file -- 
   * Do it the "hard" way */
  construct_adjacency_list ( global->ligand );
  sum_charges(global->ligand);
  find_hyd_atoms(global->ligand, &global->hyd_atom_rules);
  find_cycles(global->ligand);
  find_flexible_bonds(global->ligand, global->flex_bond_rules,
                        global->number_of_flex_bond_rules);

  if(global->ligand_flag != 0) free(global->ligand_flag);
  global->ligand_flag =
    (short *) mymalloc(3 * global->ligand->number_of_atoms * sizeof(short) );

#endif
  return SUCCESS;
}
void set_global_junk(global_data_pt global)
{
  int i;
  char file[FILENAME_MAX + 1];

  strncpy(global->slide_dir, getenv("SLIDE_DIR"), FILENAME_MAX );
  if(*global->slide_dir == '\0')
    err_panic("main_score", "SLIDE_DIR environment variable not set");
#if 0
  sprintf(file, "%s/params/flex.defn", global->slide_dir );
  global->number_of_flex_bond_rules =
    read_flex_defn(file, global->flex_bond_rules);
  sprintf(file, "%s/params/hbond.defn", global->slide_dir );
  read_hyd_defn(file, &global->hyd_atom_rules);
#endif

  /* make backup of target atom positions, since the positions in
   * 'global->target_atoms' are modified when target side chains
   * are rotated during bump resolvement */
  for ( i = 0; i < global->number_of_target_atoms; i++ )
    global->orig_target_atom_act[i] = global->target_atoms[i].act;
  memcpy(global->orig_target_atom_positions, global->target_atom_positions,
         3*global->number_of_target_atoms *
         sizeof(*global->target_atom_positions));
  memset(global->target_rotations, 0, global->number_of_target_residues *
         sizeof(*global->target_rotations));

  /* this is the array for the lookup table of inter-atomic distances,
     before doing the very first bump-check after transforming a ligand,
     this array is filled with the distances between all pairs of ligand
     and target atoms, so that we will avoid most of the calls of
     'dist_fun()' during the modeling of the induced complementarity - Volker*/
  global->target_ligand_distances =
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS *
                        global->number_of_target_atoms * sizeof (float) );
  global->target_ligand_sq_dists =
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS *
                        global->number_of_target_atoms * sizeof (float) );

  /* Allocate & initialize memory for flag arrays used in the scoring function 
   */
  global->target_flag =
    (short *) mymalloc(global->number_of_target_atoms * sizeof(short));

#ifdef NON_METALBONDED_REPULSIVE
  /* Build a list of metal atom indices */
  global->number_of_metals = 0;
  global->metal_atom_indices =
    (int*) mymalloc (MAX_TARGET_METAL_ATOMS * sizeof(int));

  for(i = 0; i < global->number_of_target_atoms; i++)
    if(global->target_atoms[i].act == METAL_1 ||
       global->target_atoms[i].act == METAL_2)
      global->metal_atom_indices[global->number_of_metals++] = i;

  if(MAX_TARGET_METAL_ATOMS <= global->number_of_metals)
    err_panic2("main_score", "Number of metals exceed MAX_TARGET_METAL_ATOMS");
#endif

  global->number_of_waters = 0;
  printf("Water handling disabled in current version of SLIDE => #waters = 0"
         "\n\n");
}

#endif
