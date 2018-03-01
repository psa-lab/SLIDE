#include <stdlib.h>
#include <malloc.h>
#include "err_handle.h"

void  *mymalloc ( size_t  size )
{
  void  *adr;

  adr = malloc ( size );
  if ( !adr )
    err_panic ( "", "malloc failed" );
  return ( adr );
}

void  *myrealloc ( void    *old,
		   size_t  size )
{
  void  *adr;

  adr = realloc ( old, size );
  if ( !adr )
    err_panic ( "", "realloc failed" );
  return ( adr );
}

