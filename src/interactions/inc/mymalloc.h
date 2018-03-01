#ifndef _MYMALLOC_H
#define _MYMALLOC_H

#include <stdlib.h>
#include <malloc.h>

extern void  *mymalloc ( size_t  size );

extern void  *myrealloc ( void    *old,
			  size_t  size );

#endif

