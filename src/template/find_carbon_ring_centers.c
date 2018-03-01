/*
 *
 *   find_carbon_ring_centers.c  Volker Schnecke  Fri Feb 27 22:06:49 EST 1998
 *
 *   functions: check_for_ring()
 *              find_cycle_points()
 *              find_carbon_ring_centers()
 *
 *   Finds hydrophobic interaction points in a molecule using the new
 *   method.  
 *
 *   This code is derived from the original interaction point code
 *   as written by Volker Schnecke, hence the not quite correct
 *   naming. All "carbon_ring_centers" are derived from the assignment
 *   algorithm and are not truly ring centers. 
 *             -- Paul Sanschagrin, 06-Jun-01
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "types.h"
#include "find_carbon_centers.h"
#include "complete_link_clustering.h"
#include "basics.h"
#include <mymalloc.h>

#define TRACE

/*
 *  This recursive routine does a depth-first search in the subgraph of 
 *  the molecule graph consisting only of carbons to identify cycles with
 *  up to MAX_CARBON_RING_SIZE atom. Once a ring has been found, the
 *  recursion immediately steps back and stores all node indices on the
 *  way back to the first node in the cycle in the array cycle.
 */
int  check_for_ring ( molecule_pt  molecule,
		      int          *mark,
		      int          *cycle,
		      int          index,
		      int          prohibited,
		      int          count )
{
  int  *neighbors;
  int  result,
       old_mark;
  int  i;

  if ( mark[index] < -1 )
    /* no carbon, so don't go any further */
    return NO_CYCLE;
  if ( count <= MAX_CARBON_RING_SIZE )
    /* go on only if we haven't visited MAX_CARBON_RING_SIZE atoms yet */
    {
      /* mark atom 'index' as atom number 'count' in the potential ring,
	 but save original mark */
      old_mark = mark[index];
      mark[index] = count;
      neighbors = molecule->neighbors[index];
      result = NO_CYCLE;
      for ( i = 0; i < molecule->number_of_neighbors[index]; i++ )
	/* check all neigbors */
	if ( mark[neighbors[i]] > -2 
	     && mark[neighbors[i]] <= 1 
	     && neighbors[i] != prohibited )
	  /* consider only carbons, if you bump into an atom of the current
	     ring other than the first one, stop the search, and don't go 
	     directly back to where you came from */
	  {
	    if ( mark[neighbors[i]] == 1 )
	      /* we have found the first node in our cycle */
	      result = CYCLE;
	    else
	      /* check the next neighbor */
	      result = check_for_ring ( molecule,
					mark,
					cycle,
					neighbors[i],
					index,
					count + 1 );
	    if ( result == CYCLE )
	      {
		/* store index of this ring atom in array 'cycle' */
		cycle[count-1] = index;
		/* this is to mark that this atom is in a cycle that has been
		   already found, i.e. no new search can start at this atom
		   to prevent marking the same cycle twice or more */
		mark[index] = -1;
		return CYCLE;
	      }
	  }
      /* this atom is not part of a cycle, so give back it's original mark */
      mark[index] = old_mark;
    }
  return NO_CYCLE;
}

/*  This recursive routine searches around each ring to place a point at 
    approximately every other atom in a ring based on starting from atoms
    in the ring which have substituents.
*/

void find_cycle_points ( molecule_pt molecule, 
			 int* in_current_ring, 
			 int* substituents,
			 int* visited,
			 int step, 
			 int index, 
			 int previndex ) 
{
  atom_pt atoms;
  int     *neighbors;
  int     i, j, k;

  atoms = molecule->atoms;
  neighbors = molecule->neighbors[index];
  visited[index] = VISITED;

  for (i = 0; i < molecule->number_of_neighbors[index]; i++) 
    {
      if ( ( in_current_ring[neighbors[i]] == TRUE ) && 
	   ( neighbors[i] != previndex ) )
	/* The atom is in the current hphobic ring and isn't the one
	   we're coming from. */
	{
#ifdef TRAC
	  printf ("Atom %4s: m=%2d; Neighbor %4s: m=%2d\n",
		  atoms[index].name,
		  in_current_ring[index],
		  atoms[neighbors[i]].name,
		  in_current_ring[neighbors[i]]);
	  
#endif
	  if ( ( substituents[neighbors[i]] == NO_SUBSTITUENTS ) &&
	       ( visited[neighbors[i]] == UNVISITED) &&
	       ( step == 2) ) 
	    /* We only name an hphobic point every ~other ring atom */
	    {
	      atoms[neighbors[i]].hyd = HPHOB;
	      find_cycle_points ( molecule,
				  in_current_ring,
				  substituents,
				  visited,
				  step-1,
				  neighbors[i],
				  index );
	    } 
	  else if ( ( substituents[neighbors[i]] == NO_SUBSTITUENTS ) &&
		    (visited[neighbors[i]] == UNVISITED) &&
		    ( step == 1 ) )
	    {
	      find_cycle_points ( molecule,
				  in_current_ring,
				  substituents,
				  visited,
				  step+1,
				  neighbors[i],
				  index );
	    }
	  
	}
    }
  
}

/*  This routine searches all carbon rings with up to MAX_CARBON_RING_SIZE
 *  atoms. It is not very efficently programmed, but since it is part of
 *  the global preprocessing when initializing the screening database, it 
 *  has to be only done once per database. Markings for all atoms are
 *  stored in an array 'mark', the entries have the following meanings:
 *   -2   no carbon, ignore this atom during the search
 *   -1   carbon already included in a cycle
 *   0    an unvisited carbon or carbon not included in a cycle yet
 *  The reason for distinguishing between the last two cases is that
 *  a search will never start with an atom that is already included in 
 *  another ring to prevent identifying the same ring more than once.
 */
void  find_carbon_ring_centers ( molecule_pt  molecule )
{
  atom_pt  atoms;
  bond_pt  bonds;
  float    avg[3];
  int      mark[MAX_NUMBER_OF_MOL2_ATOMS],
           visited[MAX_NUMBER_OF_MOL2_ATOMS],
           cycle[MAX_CARBON_RING_SIZE],
           in_current_ring[MAX_NUMBER_OF_MOL2_ATOMS],
           substituents[MAX_NUMBER_OF_MOL2_ATOMS],
           inrings_atoms[MAX_NUMBER_OF_MOL2_ATOMS],
           inrings_bonds[MAX_NUMBER_OF_MOL2_BONDS],
           *neighbors;
  int      count, 
           number_of_carbon_rings,
           number_of_clusters;
  int      i, j, k, l;
  point_pt points_to_cluster,
           clusters;
  
  atoms = molecule->atoms;
  bonds = molecule->bonds;
  number_of_carbon_rings = 0;

  /* initialization */  
  for ( i = 0; i < molecule->number_of_bonds; i++ )
      inrings_bonds[i] = 0;
  for ( i = 0; i < molecule->number_of_atoms; i++ ) 
    in_current_ring[i] = FALSE;
  for ( i = 0; i < molecule->number_of_atoms; i++ ) 
    {
      substituents[i] = NO_SUBSTITUENTS;
      visited[i] = UNVISITED;
      inrings_atoms[i] = 0;
      /* mark all carbons OR sulfurs */
      if (( atoms[i].type == C ) || ( atoms[i].type == S ))
	mark[i] = 0;
      else
	mark[i] = -2;
  }

#ifdef TRAC
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    printf ("%4s: %2d (%1d)\n",
	    atoms[i].name,
	    mark[i],
	    substituents[i]);
#endif

  for ( i = 0; i < MAX_CARBON_RING_SIZE; i++ )
    cycle[i] = -1;
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      if ( mark[i] == 0 )
	/* atom i is a carbon that is not included in any of the rings
	   found yet */
	{
	  /* mark it as the first atom in the potential cycle */
	  mark[i] = 1;
	  neighbors = molecule->neighbors[i];
	  for ( j = 0; 
		j < molecule->number_of_neighbors[i]; 
		j++ )
	    if ( mark[neighbors[j]] != -2 )
	      /* neigbor is carbon */
	      if ( check_for_ring ( molecule,
				    mark,
				    cycle,
				    neighbors[j],
				    i,
				    2 ) == CYCLE )
		/* a ring starting at atom 'i' has been found */
		{
		  cycle[0] = i;
		  /* do not check the other neighbors, we already found
		     one ring */
		  break;
		}
	  /* set mark back to non-ring carbon */
	  mark[i] = 0;
	}
      if ( cycle[0] >= 0 )
	/* a cycle has been found */
	{        
#ifdef TRAC
	  printf ( "+++++ A NEW RING +++++\n");
#endif   
	  /* ok, mark it as a ring carbon */
	  mark[i] = -1;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	    if (cycle[j] >= 0)
	      in_current_ring[cycle[j]] = TRUE;
	  count = 0;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++ )
	    if ( cycle[j] >= 0 )
	      {
		inrings_atoms[cycle[j]]++;
		count++;
		neighbors = molecule->neighbors[cycle[j]];
		for ( k = 0; k < molecule->number_of_neighbors[cycle[j]]; 
		      k++ ) {
		  /* Check for nonH ring substrituents at each postion
		     in the ring. Included are atoms in other
		     rings. */
		  if ( ((mark[neighbors[k]] == 0) ||
			((mark[neighbors[k]] == -1) && 
			 (in_current_ring[neighbors[k]] == FALSE))) && 
		       (substituents[cycle[j]] != HPHIL_SUBSTITUENT))
		    substituents[cycle[j]] = HPHOB_SUBSTITUENT;
		  else if ( (mark[neighbors[k]] == -2) && 
			    (atoms[neighbors[k]].type != H) )
		    substituents[cycle[j]] = HPHIL_SUBSTITUENT;
		}
#ifdef TRAC
		printf ( "%2d: %4s %7.3f   %7.3f   %7.3f: ",
			 cycle[j],
			 atoms[cycle[j]].name,
			 atoms[cycle[j]].pos[X],
			 atoms[cycle[j]].pos[Y],
			 atoms[cycle[j]].pos[Z] );
		for (k = 0; k < molecule->number_of_neighbors[cycle[j]]; 
		     k++ ) {
		  printf ("%4s (%2d) ",
			  atoms[neighbors[k]].name,
			  mark[neighbors[k]]);
		}
		printf("%2d  %2d\n",
		       substituents[cycle[j]],
		       inrings_atoms[cycle[j]]);
		
		
#endif
		if (substituents[cycle[j]] != NO_SUBSTITUENTS)
		  atoms[cycle[j]].hyd = HPHOB;
		if ( j == 0 )
		  /* store index of one atom on the ring */
		  molecule->carbon_ring_atom[number_of_carbon_rings] = 
		    cycle[j];
	      }
	  /* Now that we have a complete ring assigned, assign the points
	     which belong at unsubstiuated positions */
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	    if ( ( substituents[cycle[j]] != NO_SUBSTITUENTS ) &&
		 ( atoms[cycle[j]].hyd == HPHOB ) &&
		 ( cycle[j] >= 0))
	      find_cycle_points ( molecule, 
				  in_current_ring, 
				  substituents,
				  visited,
				  1,			       
				  cycle[j],
				  0);
	  
	  /* Now we have to cluster the points which are in unsubstituated
	     positions that are within a bond of each other */
	  count = 0;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	      if ( cycle[j] >= 0 )
		  in_current_ring[cycle[j]] = FALSE;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	    if ( ( substituents[cycle[j]] == NO_SUBSTITUENTS ) &&
		 ( atoms[cycle[j]].hyd == HPHOB ) &&
		 ( cycle[j] >= 0) ) 
		count++;
	  points_to_cluster = (point_pt) 
	    mymalloc ( MAX_CARBON_RING_SIZE * sizeof (point_t) );
	  clusters = (point_pt)
	    mymalloc ( MAX_CARBON_RING_SIZE * sizeof (point_t) );
	  k = 0;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	    if ( ( substituents[cycle[j]] == NO_SUBSTITUENTS ) &&
		 ( atoms[cycle[j]].hyd == HPHOB ) &&
		 ( cycle[j] >= 0) ) 
	      {
		atoms[cycle[j]].hyd = NOTHING;
		points_to_cluster[k].pos[X] = atoms[cycle[j]].pos[X];
		points_to_cluster[k].pos[Y] = atoms[cycle[j]].pos[Y];
		points_to_cluster[k].pos[Z] = atoms[cycle[j]].pos[Z];
		points_to_cluster[k].type = NOTHING;
		k++;
	      }
#ifdef TRAC
	  for (j = 0; j < count; j ++)
	    printf("%3d:  %8.3f %8.3f %8.3f\n",
		   j,
		   points_to_cluster[j].pos[X],
		   points_to_cluster[j].pos[Y],
		   points_to_cluster[j].pos[Z]);
#endif
	  number_of_clusters = 
	    complete_link_clustering ( points_to_cluster,
					  count,
					  clusters,
					  RING_THRESHOLD );
	  count = 0;
	  /* Transfer the cluster coords to carbon "ring" centers */
	  for ( j = 0; j < number_of_clusters; j++ )
	    {
		molecule->carbon_ring_centers[number_of_carbon_rings + j][X] =
		  clusters[j].pos[X];
		molecule->carbon_ring_centers[number_of_carbon_rings + j][Y] =
		  clusters[j].pos[Y];
		molecule->carbon_ring_centers[number_of_carbon_rings + j][Z] =
		  clusters[j].pos[Z];
		count ++;
	    }
	  number_of_carbon_rings += count;
	  free(points_to_cluster);
	  free(clusters);
	  
	  /* Label bonds as being a member of a ring based on their 
	     neighboring atoms being members of a ring. */
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++ ) 
	    if (cycle[j] >= 0)
	      for (l = 0; l < molecule->number_of_bonds; l++)
		  if ( (bonds[l].atom1 == cycle[j]) &&
		       (inrings_atoms[bonds[l].atom2] != 0) )
		    inrings_bonds[l]++;
	  
#ifdef TRAC
    for ( i = 0; i < molecule->number_of_bonds; i++ )
    {
      printf ("%3d (%3d): %4s  %4s  %3d\n",
	      i,
	      bonds[i].number,
	      atoms[bonds[i].atom1].name,
	      atoms[bonds[i].atom2].name,
	      inrings_bonds[i] );
    }
#endif
    
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++ ) 
	      if ( cycle[j] >= 0 )
		  inrings_atoms[cycle[j]] = 0;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++)
	    cycle[j] = -1;
	}
    }
  /* Below is to strip the HPHOBIC label from the ring atoms which
     have hydrophilic substituents */
  for (i = 0; i < molecule->number_of_atoms; i++)
    if ( ( mark[i] == -1 ) &&
	 ( substituents[i] == HPHIL_SUBSTITUENT ))
      atoms[i].hyd = NOTHING;

     
  /* Next we have to cluster points which are on atoms which share an edge with
     other rings. */
  for (i = 0; i < molecule->number_of_bonds; i++)
      {
      if (inrings_bonds[i] > 1) 
	  /* This bond is contained in more than 1 ring */
	  {
	    molecule->carbon_ring_centers[number_of_carbon_rings][X] =
	      (atoms[bonds[i].atom1].pos[X] + atoms[bonds[i].atom2].pos[X]) / 2;
	    molecule->carbon_ring_centers[number_of_carbon_rings][Y] =
	      (atoms[bonds[i].atom1].pos[Y] + atoms[bonds[i].atom2].pos[Y]) / 2;
	    molecule->carbon_ring_centers[number_of_carbon_rings][Z] =
	      (atoms[bonds[i].atom1].pos[Z] + atoms[bonds[i].atom2].pos[Z]) / 2;
	    number_of_carbon_rings++;
	    /* We must also remove the hydrophobic labels from the atoms of
	       this bond */
	    atoms[bonds[i].atom1].hyd = NOTHING;
	    atoms[bonds[i].atom2].hyd = NOTHING;
	  }
      }

  molecule->number_of_carbon_rings = number_of_carbon_rings;
  find_carbon_centers ( molecule, mark );

#define FINAL_CLUSTERING
#ifdef FINAL_CLUSTERING
  /* We must now do a final clustering step to solve some odd
     assignments in edge sharing rings which combine to a size of less
     than MAX_CARBON_RING_SIZE. */
  points_to_cluster = (point_pt) 
    mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (point_t) );
  clusters = (point_pt)
    mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (point_t) );
  count = 0;

  /* First we assign the hphobic atom points. */
  for (i = 0; i < molecule->number_of_atoms; i++)
    if (atoms[i].hyd == HPHOB)
      {
	points_to_cluster[count].pos[X] = atoms[i].pos[X];
	points_to_cluster[count].pos[Y] = atoms[i].pos[Y];
	points_to_cluster[count].pos[Z] = atoms[i].pos[Z];
	atoms[i].hyd = NOTHING;
	count++;
      }
  /* Then we assign the hphobic ring points. */
  for (i = 0; i < molecule->number_of_carbon_rings; i++)
    {
      points_to_cluster[count].pos[X] = 
	molecule->carbon_ring_centers[i][X];
      points_to_cluster[count].pos[Y] = 
	molecule->carbon_ring_centers[i][Y];
      points_to_cluster[count].pos[Z] = 
	molecule->carbon_ring_centers[i][Z];
      count++;
    }
#ifdef TRAC
  for (j = 0; j < count; j ++)
    printf("%3d:  %8.3f %8.3f %8.3f\n",
	   j,
	   points_to_cluster[j].pos[X],
	   points_to_cluster[j].pos[Y],
	   points_to_cluster[j].pos[Z]);
#endif
  number_of_clusters = 
    complete_link_clustering ( points_to_cluster,
			       count,
			       clusters,
			       FINAL_THRESHOLD );
#ifdef TRAC
  printf ("Number of Clusters: %3d\n",
	  number_of_clusters);
#endif
  /* Transfer the cluster coords to carbon "ring" centers */
  for ( j = 0; j < number_of_clusters; j++ )
    {
      molecule->carbon_ring_centers[j][X] =
	clusters[j].pos[X];
      molecule->carbon_ring_centers[j][Y] =
	clusters[j].pos[Y];
      molecule->carbon_ring_centers[j][Z] =
	clusters[j].pos[Z];
    }
  molecule->number_of_carbon_rings = number_of_clusters;
  free(points_to_cluster);
  free(clusters);

#endif  
}
