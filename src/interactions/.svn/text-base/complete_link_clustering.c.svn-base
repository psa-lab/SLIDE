#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "distance.h"
#include "basics.h"
#include "mymalloc.h"
#include "err_handle.h"

int  complete_link_clustering ( point_pt  points,
				int       number_of_points,
				point_pt  cluster_points,
				double    threshold )
{
  double  **proximity_matrix;
  double  dist,
          min_dist;
  int     **cluster;
  int     *number_of_cluster_points,
          *marker;
  int     min_i, 
          min_j,
          number_of_clusters;
  int     i, j, k;

  proximity_matrix = 
    (double **) mymalloc ( number_of_points * sizeof (double *) );
  cluster =
    (int **) mymalloc ( number_of_points * sizeof (int *) );
  number_of_cluster_points =
    (int *) mymalloc ( number_of_points * sizeof (int) );
  marker =
    (int *) mymalloc ( number_of_points * sizeof (int) );
  for ( i = 0; i < number_of_points; i++ )
    {
      proximity_matrix[i] = 
	(double *) mymalloc ( number_of_points * sizeof (double) );
      cluster[i] = 
	(int *) mymalloc ( number_of_points * sizeof (int) );
      cluster[i][0] = i;
      number_of_cluster_points[i] = 1;
      marker[i] = IN;
    }
  min_dist = 999999.9;
  for ( i = 0; i < number_of_points; i++ )
    if ( marker[i] == IN )
      for ( j = i + 1; j < number_of_points; j++ )
	if ( marker[j] == IN )
	  {
	    dist = distance ( points[i].pos, points[j].pos );
	    proximity_matrix[i][j] = dist;
	    if ( compare_double( dist, min_dist ) == -1 )
	      {
		min_dist = dist;
		min_i = i;
		min_j = j;
	      }
	  }
  while ( min_dist < threshold )
    {      
      for ( i = 0; i < number_of_cluster_points[min_j]; i++ )	
	{
	  cluster[min_i][number_of_cluster_points[min_i]] = 
	    cluster[min_j][i];
	  number_of_cluster_points[min_i]++;
	}
      marker[min_j] = OUT;
      for ( i = 0; i < min_i; i++ )
	if ( marker[i] != OUT )
	    if ( compare_double(proximity_matrix[i][min_i], proximity_matrix[i][min_j] ) == -1 )
	    proximity_matrix[i][min_i] = proximity_matrix[i][min_j];

      for ( i = min_i + 1; i < min_j; i++ )
	if ( marker[i] != OUT )
	    if ( compare_double(proximity_matrix[min_i][i], proximity_matrix[i][min_j] ) == -1 )
	    proximity_matrix[min_i][i] = proximity_matrix[i][min_j];

      for ( j = min_j + 1; j < number_of_points; j++ )
	if ( marker[j] != OUT )
	    if ( compare_double(proximity_matrix[min_i][j], proximity_matrix[min_j][j] ) == -1 )
	    proximity_matrix[min_i][j] = proximity_matrix[min_j][j];

      min_dist = 999999.9;

      for ( i = 0; i < number_of_points; i++ )
	if ( marker[i] == IN )
	  for ( j = i + 1; j < number_of_points; j++ )
	    if ( marker[j] == IN )
		if ( compare_double(proximity_matrix[i][j], min_dist ) == -1 )
		{
		  min_dist = proximity_matrix[i][j];
		  min_i = i;
		  min_j = j;
		}
    }
  number_of_clusters = 0;
  for ( i = 0; i < number_of_points; i++ )
    if ( marker[i] != OUT )
      {
	for ( k = 0; k < 3; k++ )
	  cluster_points[number_of_clusters].pos[k] = 0.0;
	for ( j = 0; j < number_of_cluster_points[i]; j++ )
	  for ( k = 0; k < 3; k++ )
	    cluster_points[number_of_clusters].pos[k] +=
	      points[cluster[i][j]].pos[k];
	for ( k = 0; k < 3; k++ )
	  cluster_points[number_of_clusters].pos[k] /= 
	    number_of_cluster_points[i];
	/* Use the int type field to hold the cluster density */
        cluster_points[number_of_clusters].type = 
          number_of_cluster_points[i];

	number_of_clusters++;
	if ( number_of_clusters >= MAX_NUMBER_OF_CLUSTERING_POINTS )
	  err_panic ( "complete_link_clustering",
		      "overflow cluster-point array" );
      }
  for ( i = 0; i < number_of_points; i++ )
    {
      free ( cluster[i] );
      free ( proximity_matrix[i] );
    }
  free ( cluster );
  free ( number_of_cluster_points );
  free ( proximity_matrix );
  return number_of_clusters;
}
