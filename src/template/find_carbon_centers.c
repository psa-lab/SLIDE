#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "complete_link_clustering.h"
#include <mymalloc.h>
#include "basics.h"

#define TRAC

void  find_carbon_centers ( molecule_pt  molecule,
			    int          *mark )
{
  atom_pt  atoms;
  int      number_of_carbon_rings;
  int      *neighbors;
  int      terminal[MAX_NUMBER_OF_MOL2_ATOMS];
  int      used[MAX_NUMBER_OF_MOL2_ATOMS];
  int      number_of_neighbors;
  int      i, j, k;
  int      count;
  int      countH;
  int      neighborcount;
  int      numberleft;
  int      nonHneighbor;
  int      usedneighbor;
  int      step = 1;
  int      number_of_clusters;
  float    avg[3];
  point_pt points_to_cluster,
           clusters;
  char     linebuffer[MAX_MOL2_LINELENGTH];

  number_of_carbon_rings = molecule->number_of_carbon_rings;
  atoms = molecule->atoms;

  /* mark all atoms that have only hydrophobic neighbors (i.e., only C
     or S) */
  for ( i = 0; i < molecule->number_of_atoms; i++ )
      {
	  terminal[i] = FALSE;
	  used[i] = FALSE;
      }

  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      if ( mark[i] == 0 )
	/* when searching for carbon ring centers, all carbons were
	   originally marked with 0, all other atoms with -2, and
	   after the ring search, all atoms within a ring of
	   MAX_CARBON_RING_SIZE carbons are marked with -1, thus all
	   atoms that are marked with 0 at this point are carbons not
	   included in a ring system, the center of which is already
	   considered as a hydrophobic interaction center */
	{
	  neighbors = molecule->neighbors[i];
	  number_of_neighbors = molecule->number_of_neighbors[i];
	  count = 0; 
	  countH = 0;
	  for (j = 0; j < number_of_neighbors; j++) 
	    {
	      if ( (atoms[neighbors[j]].type == C) ||
		   (atoms[neighbors[j]].type == S) )
		count ++;
	      if ( atoms[neighbors[j]].type == H)
		countH ++;
	      
	      if ( count == (number_of_neighbors - countH) )
		/* all nonH neighbors are either carbons or sulfurs, and
		   this is a carbon not located on a carbon ring -> it
		   is part of a hydrophobic network */
		mark[i] = 2;
	      if ( (number_of_neighbors-countH) == 1 && mark[i] == 2)
		terminal[i] = TRUE;
	    }
	}
      }
  
#ifdef TRACE
  for (i = 0; i < molecule->number_of_atoms; i++ )
    {
      printf ("%3d (%4s): %2d (m) %1d (t) -- ",
	      i,
	      molecule->atoms[i].name,
	      mark[i],
	      terminal[i] );
      neighbors = molecule->neighbors[i];
      number_of_neighbors = molecule->number_of_neighbors[i];
      for (j = 0; j < number_of_neighbors; j++) 
	printf ("%4s ",atoms[neighbors[j]].name);
      printf ("\n");
    }
#endif  

  for (i = 0; i < molecule->number_of_atoms; i++ )
      if ( (terminal[i] == TRUE) && (used[i] == FALSE) )
	{
#ifdef TRACE
	  printf ("T Atom %3d (%4s): ",
		  i,
		  atoms[i].name );
#endif
	  /* Now, for each of the terminal atoms, we assign a point as
	     the average of the terminal and T-1 atoms */
	  avg[X] = avg[Y] = avg[Z] = 0.0;
	  neighbors = molecule->neighbors[i];
	  count = 0;
	  nonHneighbor = -1;
	  /* First we have to find which neighbor of the terminal atom
	     is the single nonH, C/S atom. */
	  for (j = 0; j < molecule->number_of_neighbors[i]; j++)
	    {
	      if (atoms[neighbors[j]].type != H) 
		{
		  nonHneighbor = neighbors[j];
		  break;
		}	      
	    }
	  if ( nonHneighbor == -1)
	    {
	      sprintf (linebuffer,
		       "NO NON-H NEIGHBORS FOR TERMINAL ATOM: %3d (%4s)\n",
		       i,
		       atoms[i].name);
	      err_panic ("find_carbon_center",linebuffer);
	    }

	  /* Now we must check the number of termianl nonH-neighbors
	     of the T-1 atom */
	  for (j = 0; j < molecule->number_of_neighbors[nonHneighbor]; j++)
	    if ( (atoms[molecule->neighbors[nonHneighbor][j]].type != H) &&
		 (terminal[molecule->neighbors[nonHneighbor][j]] == TRUE) )
	      count++;

	  /* If count == 1, the atom i is a single terminal, but we
	     have to check if the T-1 atom is also in the nphobic
	     network. If so, we still assign a point to the mean T/T-1
	     position */
	  if ( count == 1 )
	    {
	      if ( mark[nonHneighbor] == 2 )
		/* This is a single terminal bonded to an HN atom --
		   assign a point as the average between the two */
		{
		  molecule->carbon_ring_centers[number_of_carbon_rings][X] = 
		    (atoms[i].pos[X] + atoms[nonHneighbor].pos[X])/ 2;
		  molecule->carbon_ring_centers[number_of_carbon_rings][Y] = 
		    (atoms[i].pos[Y] + atoms[nonHneighbor].pos[Y])/ 2;
		  molecule->carbon_ring_centers[number_of_carbon_rings][Z] = 
		    (atoms[i].pos[Z] + atoms[nonHneighbor].pos[Z])/ 2;
		  number_of_carbon_rings++;
		  used[i] = TRUE;
		  used[nonHneighbor] = TRUE;
#ifdef TRACE
		  printf ("Single terminal, assigned average point\n");
#endif
		}
	      else if ( mark[nonHneighbor] == -1) 
		{
		  atoms[i].hyd = HPHOB;
		  used[i] = TRUE;
#ifdef TRACE
		  printf ("Single terminal, assigned atom point (i.e. ring)\n");
#endif
		}
	      else 
		{
		  used[i] = TRUE;
#ifdef TRACE
		  printf (" Single Terminal, no point assigned\n");
#endif
		}
	    }
	  else
	    /* This is not a single terminal, so we assign a point at
	       the average of the terminal and T-1 atoms. */
	    {
	      avg[X] = atoms[nonHneighbor].pos[X];
	      avg[Y] = atoms[nonHneighbor].pos[Y];
	      avg[Z] = atoms[nonHneighbor].pos[Z];
	      for (j = 0; j < molecule->number_of_neighbors[nonHneighbor]; j++)
		if ( (mark[molecule->neighbors[nonHneighbor][j]] == 2) &&
		     (terminal[molecule->neighbors[nonHneighbor][j]] == TRUE))
		  {
#ifdef TRACE
		    printf("%3d (%4s) ",
			   molecule->neighbors[nonHneighbor][j],
			   atoms[molecule->neighbors[nonHneighbor][j]].name);
#endif
		    avg[X] += 
		      atoms[molecule->neighbors[nonHneighbor][j]].pos[X];
		    avg[Y] += 
		      atoms[molecule->neighbors[nonHneighbor][j]].pos[Y];
		    avg[Z] += 
		      atoms[molecule->neighbors[nonHneighbor][j]].pos[Z];
		    used[molecule->neighbors[nonHneighbor][j]] = TRUE;
		  }
	      /* The count+1 below comes from the fact that the count
		 only includes the NonH neighbors of the T-1 atom and
		 not the T-1 atom itself. */
	      molecule->carbon_ring_centers[number_of_carbon_rings][X] = 
		avg[X] / (count + 1);
	      molecule->carbon_ring_centers[number_of_carbon_rings][Y] = 
		avg[Y] / (count + 1);
	      molecule->carbon_ring_centers[number_of_carbon_rings][Z] = 
		avg[Z] / (count + 1);
	      number_of_carbon_rings++;
	      used[nonHneighbor] = TRUE;
#ifdef TRACE
	      printf("\n");
#endif
	      used[nonHneighbor] = TRUE;
	    }
	}

  /* Now that we have assigned all of the terminal Hphobic points, we
     must assign the internal ones. This will use a similar procedure
     as above, with a few changes. */
  
  numberleft = MAX_NUMBER_OF_MOL2_ATOMS; /* Just setting to some
					    number >0 to force one run
					    through the loop. */
  while (numberleft != 0)
    {
      numberleft = 0;
      for (i = 0; i < molecule->number_of_atoms;i++)
	{
	  if ( (mark[i] == 2) && (used[i] != TRUE) )
	    /* We're only interested in the hphobic atoms which have
	     not been used before, but which have a used neighbor */
	    {
#ifdef TRACE
	      printf ("Atom %3d (%4s): ",
		      i,
		      atoms[i].name);
#endif
	      count = 0;
	      neighborcount = 0;
	      usedneighbor = FALSE;
	      neighbors = molecule->neighbors[i];
	      number_of_neighbors = molecule->number_of_neighbors[i];
	      for (j = 0; j < number_of_neighbors; j++)
		{
		  if ( (mark[neighbors[j]] == 2) && 
		       (used[neighbors[j]] != TRUE) )
		    count ++;
		  else if (mark[neighbors[j]] == 0)
		    /* Count any nonHN hphobic neighbors */
		    neighborcount++;
		  if (used[neighbors[j]] == TRUE)
		    usedneighbor = TRUE;
		}
	      /* We must distinguish between cases where all HN
		 neighbors of the atom have been used (and we don't
		 assign a point) and cases where the atom forms a
		 single atom HN (and we do assign a point) */
	      if (count == 0)
		{
		  used[i] = TRUE;
		  if (neighborcount > 0)
		    /* There are some nonHN hphoibic neighbors,
		       therefore this is a single atom HN and gets
		       assigned a point. This assignment could be done
		       by assigning to the atom.hyd property, but is
		       done as a "ring" to facilitate clustering done
		       later. */
		    {
		      molecule->carbon_ring_centers[number_of_carbon_rings][X] =
			atoms[i].pos[X];
		      molecule->carbon_ring_centers[number_of_carbon_rings][Y] =
			atoms[i].pos[Y];
		      molecule->carbon_ring_centers[number_of_carbon_rings][Z] =
			atoms[i].pos[Z];
		      number_of_carbon_rings++;
		    }
#ifdef TRACE
		  printf ("No unused, hphobic neighbors.\n");
#endif
		  continue;
		}
	      else if ( (usedneighbor == TRUE) && (count == 1) )
		/* We have found a new termini */
		{
#ifdef TRACE
		  for (j = 0; j < number_of_neighbors; j++)		
		    {
		      printf("%3d (%4s) ",
			     neighbors[j],
			     atoms[neighbors[j]].name);
		      if (used[neighbors[j]] == TRUE)
			printf ("USED ");
		    }
		  printf ("\n");
#endif
		}
#ifdef TRACE
	      else
		/* This is an internal piece not attached to an
		   external terminal. It gets processed correctly, but
		   TRACE output needs this additional output. */
		printf ("Internal, unattached HN piece\n");
#endif
	      /* Now again, we must find the T-1 atom and count it's
		 neighbors */
	      for (j = 0; j < number_of_neighbors; j++)	
		if ( (mark[neighbors[j]] == 2) &&
		     (used[neighbors[j]] == FALSE) )
		  /* We have found the nonH neighbor of this termini */
		  nonHneighbor = neighbors[j];

	      /* Count the number of nonH Neighbors of the T-1 atom */
	      count = 0;
	      neighbors = molecule->neighbors[nonHneighbor];
	      number_of_neighbors = molecule->number_of_neighbors[nonHneighbor];
	      for (j = 0; j < number_of_neighbors; j++)
		{
		  if ( (mark[neighbors[j]] == 2) )
		    count ++;
		}
	      /* There must be at least 1 hphobic neighbor since we
		 came from one */
	      if (count == 1)
		/* Only the one we came from => there are only 2
		   hphobic atoms in this HN */
		/* We must check which internal step we're on */
		{
		  if (step > 1)
		    /* We don't want a point on anything after the
		       first step */
		    used[i] = TRUE;
		  else
		    /* Assign to average of T and T-1 atoms */
		    {
		      molecule->carbon_ring_centers[number_of_carbon_rings][X] = 
			(atoms[i].pos[X] + atoms[nonHneighbor].pos[X])/ 2;
		      molecule->carbon_ring_centers[number_of_carbon_rings][Y] = 
			(atoms[i].pos[Y] + atoms[nonHneighbor].pos[Y])/ 2;
		      molecule->carbon_ring_centers[number_of_carbon_rings][Z] = 
			(atoms[i].pos[Z] + atoms[nonHneighbor].pos[Z])/ 2;
		      number_of_carbon_rings++;
		      used[i] = TRUE;
		      used[nonHneighbor] = TRUE;		      
		    }
		}
	      else
		/* This is a multiple internal terminal, so assign a
		   point to the mean of each of the atoms neighboring
		   the T-1 atom */
		{
		  neighborcount = 0;
		  avg[X] = atoms[nonHneighbor].pos[X];
		  avg[Y] = atoms[nonHneighbor].pos[Y];
		  avg[Z] = atoms[nonHneighbor].pos[Z];
		  for (j = 0; j < number_of_neighbors; j++)
		    if ( (mark[neighbors[j]] == 2) /*&&
						     (used[neighbors[j]] != TRUE)*/)
		      {
			avg[X] += atoms[neighbors[j]].pos[X];
			avg[Y] += atoms[neighbors[j]].pos[Y];
			avg[Z] += atoms[neighbors[j]].pos[Z];
		      }
		  molecule->carbon_ring_centers[number_of_carbon_rings][X] 
		    = avg[X] / (count + 1);			     
		  molecule->carbon_ring_centers[number_of_carbon_rings][Y] 
		    = avg[Y] / (count + 1);			     
		  molecule->carbon_ring_centers[number_of_carbon_rings][Z] 
		    = avg[Z] / (count + 1);
		  number_of_carbon_rings++;
		  used[j] = TRUE;
		}
	      
	      used[nonHneighbor] = TRUE;
	      numberleft ++;
	    }
	}
      step ++;
    }

  /* Now that we have assigned all of these new points based on the
     hphobic networks, we must cluster them to eliminate points which
     are too close together. */
  count = 0;
  points_to_cluster = (point_pt) 
    mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (point_t) );
  clusters = (point_pt)
    mymalloc ( MAX_NUMBER_OF_MOL2_ATOMS * sizeof (point_t) );

  /* We only want to cluster the new points, not those derived from
     ring structures. This is from
     rings[molecule->number_of_carbon_rings]
     (molecule->number_of_carbon_rings hasn't changed since it gets
     changed only at the end) to rings[total number of rings] */
  /*number_of_new_carbon_rings = number_of_rings - 
    molecule->number_of_carbon_rings;*/
  for ( i = molecule->number_of_carbon_rings; 
	i < number_of_carbon_rings; i++)
    {
      points_to_cluster[i - molecule->number_of_carbon_rings].pos[X] 
	= molecule->carbon_ring_centers[i][X];
      points_to_cluster[i - molecule->number_of_carbon_rings].pos[Y] 
	= molecule->carbon_ring_centers[i][Y];
      points_to_cluster[i - molecule->number_of_carbon_rings].pos[Z] 
	= molecule->carbon_ring_centers[i][Z];
    }
  
#ifdef TRACE
  for ( i = 0; i < (number_of_carbon_rings - 
		    molecule->number_of_carbon_rings); i++)
    printf("%3d:  %8.3f %8.3f %8.3f\n",
	   i,
	   points_to_cluster[i].pos[X],
	   points_to_cluster[i].pos[Y],
	   points_to_cluster[i].pos[Z]);
#endif
  
  number_of_clusters = 
    complete_link_clustering ( points_to_cluster,
			       (number_of_carbon_rings - 
				molecule->number_of_carbon_rings),
			       clusters,
			       OVERALL_THRESHOLD );
#ifdef TRACE
  printf ("Clusters: %3d\n",number_of_clusters);
  for (i = 0; i < number_of_clusters; i++)
    printf("%3d:  %8.3f %8.3f %8.3f\n",
	   i,
	   clusters[i].pos[X],
	   clusters[i].pos[Y],
	   clusters[i].pos[Z]);
#endif

  /* Now we have the new clusters, we must assign a point to each
     center and "delete" the previous points. */
  for ( i = 0; i < number_of_clusters; i++)
    {
      molecule->carbon_ring_centers[i + molecule->number_of_carbon_rings][X] 
	= clusters[i].pos[X];
      molecule->carbon_ring_centers[i + molecule->number_of_carbon_rings][Y] 
	= clusters[i].pos[Y];
      molecule->carbon_ring_centers[i + molecule->number_of_carbon_rings][Z] 
	= clusters[i].pos[Z];
    }

  number_of_carbon_rings = number_of_clusters;
  molecule->number_of_carbon_rings += number_of_carbon_rings;
}
