#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "defs.h"
#include "assign_type.h"
#include <mymalloc.h>
#include <err_handle.h>
#include <basics.h>

void  read_mol2 ( char        *filename, molecule_pt molecule )
{
  char  linebuffer[MAX_MOL2_LINELENGTH];
  FILE  *in;
  int   i;

  in = slide_fopen(filename, "r");
  if(in == NULL){
    sprintf ( linebuffer, "unable to open file %s", filename );
    printf("No template was generated!\n");  
    err_panic ( "read_mol2", linebuffer );
  }
  molecule->score = 0.0;
  do
    {
      fgets ( linebuffer, sizeof ( linebuffer ), in );
      if ( strncmp ( linebuffer, "# score", 7 ) == 0 )
	sscanf ( linebuffer, "%*s %*s %*s %f", &molecule->score );
    }
  while ( linebuffer[0] != '@' );
  do
    {
      if ( strncmp ( linebuffer, "@<TRIPOS>MOLECULE", 17 ) == 0 )
	{
	  /* read molecule name and number of atoms and bonds */
	  fgets ( linebuffer, sizeof ( linebuffer ), in );
	  sscanf ( linebuffer, 
		   "%s", 
		   molecule->name );
	  fscanf ( in, 
		   "%d %d %*d", 
		   &molecule->number_of_atoms, 
		   &molecule->number_of_bonds );
	  if ( molecule->number_of_atoms > MAX_NUMBER_OF_MOL2_ATOMS )
	    err_panic ( "read_mol2", "more than MAX_NUMBER_OF_MOL2_ATOMS" );
	  if ( molecule->number_of_bonds > MAX_NUMBER_OF_MOL2_BONDS )
	    err_panic ( "read_mol2", "more than MAX_NUMBER_OF_MOL2_BONDS" );
	}
      if ( strncmp ( linebuffer, "@<TRIPOS>ATOM", 13 ) == 0 )
	/* read all atom data */
	for ( i = 0; i < molecule->number_of_atoms; i++ )
	  {
	    /* we cannot use fscanf() here, since we only want to grep
	       six items out of a line which can have six or more entries */
	    fgets ( linebuffer, sizeof ( linebuffer ), in );
	    sscanf ( linebuffer, 
		     "%d %s %lf %lf %lf %s %*d %*s %*f",
		     &molecule->atoms[i].number,
		     molecule->atoms[i].name,
		     &molecule->atoms[i].pos[X],
		     &molecule->atoms[i].pos[Y],
		     &molecule->atoms[i].pos[Z],
		     molecule->atoms[i].type_str );
	    /* it is assumed that atoms are counted from 1 to n, but
	       the internal representation (e.g. for bonds) uses numbering
	       from 0 to n-1, so stop, if there are inconsitencies */
	    if ( molecule->atoms[i].number != i + 1 )
	      {
		fprintf ( stderr, 
			  "%s: i = %d, line = %d %s %s\n", 
			  filename,
			  i, 
			  molecule->atoms[i].number,
			  molecule->atoms[i].name,
			  molecule->atoms[i].type_str );
		err_panic ( "read_mol2", "atom numbers incorrect" );
	      }
	    assign_type_and_orbit ( molecule->atoms[i].type_str,
				    &molecule->atoms[i].type,
				    &molecule->atoms[i].orbit );
	  }
      if ( strncmp ( linebuffer, "@<TRIPOS>BOND", 13 ) == 0 )
	/* read all bond data */
	for ( i = 0; i < molecule->number_of_bonds; i++ )
	  {
	    fscanf ( in,
		     "%d %d %d %s",
		     &molecule->bonds[i].number,
		     &molecule->bonds[i].atom1,
		     &molecule->bonds[i].atom2,
		     molecule->bonds[i].type_str );
	    /* convert file numbering into internal numbering */
	    molecule->bonds[i].atom1--;
	    molecule->bonds[i].atom2--;
	    if ( molecule->bonds[i].type_str[0] == '1' )
	      molecule->bonds[i].type = SINGLE;
	    else
	      if ( molecule->bonds[i].type_str[0] == '2' )
		molecule->bonds[i].type = DOUBLE;
	      else
		if ( molecule->bonds[i].type_str[0] == '3' )
		  molecule->bonds[i].type = TRIPLE;
		else
		  if ( molecule->bonds[i].type_str[0] == '4' )
		    molecule->bonds[i].type = QUADRUPLE;
		  else
		    if ( molecule->bonds[i].type_str[0] == '7' )
		      molecule->bonds[i].type = DELOCALIZED;
		    else
		      if ( molecule->bonds[i].type_str[0] == 'a' 
			 || molecule->bonds[i].type_str[0] == 'A' )
		      molecule->bonds[i].type = AROMATIC;
		    else
		      if ( molecule->bonds[i].type_str[0] == 'u' 
			   || molecule->bonds[i].type_str[0] == 'U' )
			molecule->bonds[i].type = UNKNOWN;
		      else
			{
			  fprintf ( stderr, 
				    "bond type: %s\n",
				    molecule->bonds[i].type_str );
			  err_panic ( "read_mol2",
				      "unknown bond type" );
			}
	  }
    }
  while ( fgets ( linebuffer, sizeof ( linebuffer ), in ) != NULL );
  fclose ( in );
}
