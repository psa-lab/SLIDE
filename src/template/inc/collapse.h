#ifndef _COLLAPSE_H
#define _COLLAPSE_H

extern int valid_pt(struct _spoint *AA, struct _spoint *BB, point_t PP[], point_t *temp);
extern double find_distance(point_t AA, point_t BB);
extern int comp_cpair ( const void *a, const void *b);
extern void cmerge(int a[], int b[], int *na, int *nb);

#endif

