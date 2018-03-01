#include "mymalloc.h"
#include "err_handle.h"

void  *mymalloc ( size_t  size )
{
  void  *adr;
  adr = malloc ( size );
  if ( !adr ) err_panic( "mymalloc", "malloc failed" );
  return ( adr );
}

void  *myrealloc ( void    *old, size_t  size )
{
  void  *adr;

  adr = realloc ( old, size );
  if ( !adr ) err_panic( "myrealloc", "realloc failed" );
  return ( adr );
}

void my_free(void *ptr)
{
  if(ptr) free(ptr);
  ptr = NULL;
}
