#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "defs.h" 


int compare_double ( double d1, double d2 )
{
    if ( ( d1 > d2 ) && \
	 ( ( d1 - d2 ) > MIN_DOUBLE || \
	   ( d1 - d2 ) < -MIN_DOUBLE ) ) 
	return (1);
    else
	 if ( ( d1 < d2 ) && \
	 ( ( d1 - d2 ) > MIN_DOUBLE || \
	   ( d1 - d2 ) < -MIN_DOUBLE ) ) 
	return (-1);
    else
	return (0);
}
