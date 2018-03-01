#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"

void  assign_type_and_orbit ( char  *str,
			      int   *type,
			      int   *orbit )
{
  char  *orbit_str;

  /* check if there is a '.' that separates both parts */
  orbit_str = strchr ( str, '.' );
  if ( orbit_str == str || *str == '*' )
    /* no atom type specified */
    *type = ANY;
  else
    /* determine atom type */
    {
      if ( strncmp ( str, "Ag", 2 ) == 0 )
	*type = AG;
      else
      if ( strncmp ( str, "As", 2 ) == 0 )
	*type = AS;
      else
	if ( *str == 'B' )
	  *type = BR;
	else
	  if ( strncmp ( str, "Co", 2 ) == 0 )
	    *type = CO;
	  else
	    if ( strncmp ( str, "Cl", 2 ) == 0 )
	      *type = CL;
	    else
	    if ( strncmp ( str, "Cu", 2 ) == 0 )
	      *type = CU;
	    else
	      if ( *str == 'C' )
		*type = C;
	      else
		if ( strncmp ( str, "Er", 2 ) == 0 )
		  *type = ER;
		else
		if ( strncmp ( str, "Fe", 2 ) == 0 )
		  *type = FE;
		else
		  if ( *str == 'F' )
		    *type = F;
		  else
		  if ( *str == 'H' )
		    *type = H;
		  else
		  if ( *str == 'I' )
		    *type = I;
		  else
		    if ( *str == 'K' )
		      *type = K;
		    else
		      if ( strncmp ( str, "La", 2 ) == 0 )
			*type = LA;
		      else
		      if ( strncmp ( str, "Li", 2 ) == 0 )
			*type = LI;
		      else
		      if ( strncmp ( str, "Mn", 2 ) == 0 )
			*type = MN;
		      else
			if ( strncmp ( str, "Mo", 2 ) == 0 )
			  *type = MO;
			else
			  if ( strncmp ( str, "Ni", 2 ) == 0 )
			    *type = NI;
			  else
			    if ( *str == 'N' )
			      *type = N;
			    else
			      if ( strncmp ( str, "Os", 2 ) == 0 )
				*type = OS;
			      else
				if ( *str == 'O' )
				  *type = O;
				else
				    if ( strncmp ( str, "Pd", 2 ) == 0 )
				      *type = PD;
				    else
				  if ( *str == 'P' )
				    *type = P;
				  else
				    if ( strncmp ( str, "Re", 2 ) == 0 )
				      *type = RE;
				    else
				    if ( strncmp ( str, "Ru", 2 ) == 0 )
				      *type = RU;
				    else
				    if ( strncmp ( str, "Si", 2 ) == 0 )
				      *type = SI;
				    else
				      if ( *str == 'S' )
					*type = S;
				      else
				      if ( *str == 'T' )
					*type = TL;
				      else
					if ( *str == 'U' )
					  *type = U;
					else
                                        if ( strncmp ( str, "Zn", 2 ) == 0 )
                                          *type = ZN;
                                        else
					if ( strncmp ( str, "Du", 2 ) == 0 )
					  *type = DU;
					else
					  {
					    fprintf ( stderr,
						      "atom definition: %s\n",
						      str );
					    err_panic ( "assign_type_and_orbit",
							"unknown atom type" );
					  }
    }
  if ( orbit_str == NULL )
    /* no orbit type specified */
    *orbit = ANY;
  else
    /* determine orbit type */
    {
      orbit_str++;   /* skip '.' */
      if ( *orbit_str == ' ' || *orbit_str == '\0' )
	/* this is for rules like "( . )", i.e. something has to be
	   bound to this, but neither atom type nor orbit are specfied,
	   and for rules like "N.", where any orbit is ok */
	*orbit = ANY;
      else
	if ( *orbit_str == '1' )
	  *orbit = SP1;
	else
	  if ( *orbit_str == '2' )
	    *orbit = SP2;
	  else
	    if ( *orbit_str == '3' )
	      *orbit = SP3;
	    else
	      if ( *orbit_str == '4' )
		*orbit = SP4;
	      else
		if ( strncmp ( orbit_str, "ar", 2 ) == 0 )
		  *orbit = AR;
		else
		  if ( strncmp ( orbit_str, "cat", 3 ) == 0 )
		    *orbit = CAT;
		  else
		    if ( strncmp ( orbit_str, "co2", 3 ) == 0 )
		      *orbit = CO2;
		    else
		      if ( strncmp ( orbit_str, "oh", 2 ) == 0 )
			*orbit = OH;
		      else
		      if ( strncmp ( orbit_str, "o2", 2 ) == 0 )
			*orbit = O2;
		      else
			if ( strncmp ( orbit_str, "am", 2 ) == 0 )
			  *orbit = AM;
			else
			  if ( strncmp ( orbit_str, "pl3", 3 ) == 0 )
			    *orbit = PL3;
			  else
			    if ( strncmp ( orbit_str, "th", 2 ) == 0 )
			      *orbit = TH;
			    else
			      if ( *orbit_str == 'o' )
				*orbit = O;
			      else
				{
				  fprintf ( stderr, 
					    "atom definition: %s\n",
					    str );
				  err_panic ( "assign_type_and_orbit",
					      "unknown orbit type" );
				}
    }
}

