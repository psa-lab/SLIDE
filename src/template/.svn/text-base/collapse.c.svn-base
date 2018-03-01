/* Functions for clustering code of the hydrogen    */
/* bonding template points                          */
/*     -- RSK                                       */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

double sqr(double x)
{ 
	return x*x; 
}

/** compare function for qsort **/
int comp_cpair ( const void *a, const void *b)
{
	struct _cpair a1 = *((struct _cpair *)a);
	struct _cpair b1 = *((struct _cpair *)b);

	double diff = a1.dist - b1.dist;
	if (diff>0 && (diff > MIN_DOUBLE || diff < -MIN_DOUBLE)) return 1;
	else if (diff<0 && (diff > MIN_DOUBLE || diff < -MIN_DOUBLE)) return -1;
	else return 0;
}

/** merge two lists into one **/
void cmerge(int a[], int b[], int *na, int *nb)
{
	int i,j;
	int tmpb;

	tmpb=(*nb);

	for(i=0;i<(*na);i++)
	{ 
		b[tmpb]=a[i]; 
		tmpb++; 
	}
	for(i=0;i<(*nb);i++)
	{ 
		a[(*na)]=b[i]; 
		(*na)++; 
	}
	(*nb)=tmpb;

	/** printf("a (%d): ", (*na));
	    for(i=0;i<(*na);i++) printf("%d ", a[i]); printf("\n");
	    printf("b (%d): ", (*nb));
	    for(i=0;i<(*nb);i++) printf("%d ", b[i]); printf("\n");
	**/

}

double find_distance(point_t AA, point_t BB)
{
	return sqrt( sqr(AA.pos[0]-BB.pos[0]) + sqr(AA.pos[1]-BB.pos[1]) + sqr(AA.pos[2]-BB.pos[2]));
}

/** check if two points are same **/
int is_same_pt(point_t a, point_t b)
{
	if( a.pos[0] == b.pos[0] &&
	    a.pos[1] == b.pos[1] &&
	    a.pos[2] == b.pos[2] &&
	    a.type == b.type) return 1;
	else return 0;
}

/** integer compare function for qsort **/
static int comp_int(const void *a, const void *b)
{

	int a1 = *((int *)a);
	int b1 = *((int *)b);

	int diff = a1 - b1;
	if (diff>0) return 1;
	else if (diff<0) return -1;
	else return 0;

}

/** check if both arrays are same **/
int eq_array(int a[], int b[], int na, int nb)
{
	int i;
	int flag;

	if(na != nb) return 0;
	qsort(a, na, sizeof(int), &comp_int);
	qsort(b, nb, sizeof(int), &comp_int);

	for(i=0;i<na;i++)
		if(a[i]!=b[i]) return 0;

	return 1;
}

/* Computes the new point temp that is the mean all points in clusters AA 
and BB */ 

void get_new_point(struct _spoint *AA, struct _spoint *BB, point_t PP[], point_t *temp)
{
	int i;

	temp->pos[0] = temp->pos[1] = temp->pos[2] = 0;

	for(i=0; i<AA->nindx; i++)
	{
		temp->pos[0] += PP[AA->oindx[i]].pos[0];
		temp->pos[1] += PP[AA->oindx[i]].pos[1];
		temp->pos[2] += PP[AA->oindx[i]].pos[2];
	}
	for(i=0; i<BB->nindx; i++)
	{
		temp->pos[0] += PP[BB->oindx[i]].pos[0];
		temp->pos[1] += PP[BB->oindx[i]].pos[1];
		temp->pos[2] += PP[BB->oindx[i]].pos[2];
	}

	temp->pos[0] = temp->pos[0]/(AA->nindx+BB->nindx);
	temp->pos[1] = temp->pos[1]/(AA->nindx+BB->nindx);
	temp->pos[2] = temp->pos[2]/(AA->nindx+BB->nindx);
}


/** check for validity of new point **/
int valid_pt(struct _spoint *AA, struct _spoint *BB, point_t PP[], point_t *temp)
{
	int i, j;

	/* find the mean point */
	get_new_point(AA, BB, PP, temp);

	/** check for drift with points in cluster A **/
	for(i=0; i<AA->nindx; i++)
		if(find_distance((*temp),PP[AA->oindx[i]]) >= (HBOND_THRESHOLD/2))
			return 0;

	/** check for drift with points in cluster B **/
	for(i=0; i<BB->nindx; i++)
		if(find_distance((*temp),PP[BB->oindx[i]]) >= (HBOND_THRESHOLD/2))
			return 0;

	/** check if they belong to the same cluster **/
	/** code can be improved here **/
	/**if (is_same_pt(AA->tpnt, BB->tpnt) &&
			AA->nindx>1 && BB->nindx>1 &&
			eq_array(AA->oindx, BB->oindx, AA->nindx, BB->nindx))
			return 0;**/

	return 1;
}

/* copy the data of cluster to another */
void update_pt(struct _spoint *old, struct _spoint *new)
{
	int i;

	/* copy coordinates */
	old->tpnt.pos[0] = new->tpnt.pos[0];
	old->tpnt.pos[1] = new->tpnt.pos[1];
	old->tpnt.pos[2] = new->tpnt.pos[2];
	/* copy  points in the cluster */
	if(old->nindx < new->nindx)
	{
		for(i=0; i<new->nindx; i++)
			old->oindx[i] = new->oindx[i];
		old->nindx = new->nindx;
	}
	else fprintf(stdout, "ERROR: Merged cluster smaller than parent\n");

}
