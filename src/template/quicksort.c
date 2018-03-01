/*
 * sort the double array and put sorted order into sort_index
 * The function quicksort takes three arguments : the input array of numbers to be sorted,
 * the starting position of list, the last position of list. 
 */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>

#include "defs.h"
#include "types.h"
#include "basics.h"
#include "quicksort.h"


void swap(struct _cpair in_array[], int i, int j);

void quicksort_cpair(struct _cpair in_array[], int left, int right)
{
	int current=0, last=0;

	if (left >= right)
        	return;
	swap(in_array, left, (left+right)/2);
	last = left;
	for (current = left + 1; current <= right; current++)
	    if ((in_array[current].dist < in_array[left].dist) && ((in_array[current].dist - in_array[left].dist) > MIN_DOUBLE || (in_array[left].dist - in_array[current].dist) > MIN_DOUBLE))
			swap(in_array, ++last, current);
	swap(in_array, left, last);
	quicksort_cpair(in_array, left, last-1);
	quicksort_cpair(in_array, last+1, right);
}

/** The function swap interchanges the values in two positions of an array. **/ 

void swap(struct _cpair in_array[], int i, int j)
{
	double temp=0.0;
	struct _spoint *temp_point=NULL;
	int cnt=0;

	temp = in_array[i].dist;
	in_array[i].dist = in_array[j].dist;
	in_array[j].dist = temp;
	for (cnt=0; cnt<2; cnt++) 
	{
		temp_point = in_array[i].pnts[cnt];
		in_array[i].pnts[cnt] = in_array[j].pnts[cnt];
		in_array[j].pnts[cnt] = temp_point;
	}

}


