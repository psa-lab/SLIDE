#include <stdio.h>
#include <stdlib.h>

void  err_panic ( char  *function,
		 char  *message )
{
  fprintf ( stderr, 
	    "FATAL ERROR: %s (%s)\n", 
	    message, 
	    function );
  exit ( EXIT_FAILURE );
}
 
void  err_error ( char  *function,
		  char  *message )
{
  fprintf ( stderr, 
	    "ERROR: %s (%s)\n", 
	    message, 
	    function );
}
 
void  err_warning ( char  *function,
		    char  *message )
{
  fprintf ( stderr, 
	    "WARNING: %s (%s)\n", 
	    message, 
	    function );
}
 
void  err_usage ( char  *message )
{
  fprintf ( stderr, 
	    "USAGE: %s\n", 
	    message );
  exit ( EXIT_FAILURE );
}

	    
