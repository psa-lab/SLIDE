#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"
#include "read_rule.h"

int  read_flex_defn ( char              *filename,
		      flex_bond_defn_pt bond_rule )
{
  char  linebuffer[MAX_MOL2_LINELENGTH];
  FILE  *in;
  int   number_of_rules;
  
  in = fopen ( filename, "r" );
  if ( in == NULL ) 
    {
      sprintf ( linebuffer, "unable to open file %s", filename );
      err_panic ( "read_flex_defn", linebuffer );
    }
  number_of_rules = 0;
  while ( fgets ( linebuffer, sizeof ( linebuffer ), in ) != NULL )
    {
      /* superscribe `\n` with `\0` */
      linebuffer[strlen(linebuffer)-1] = '\0';
      if ( strncmp ( linebuffer, "name", 4 ) == 0 )
	sscanf ( linebuffer, "%*s %s", bond_rule[number_of_rules].name );      
      if ( strncmp ( linebuffer, "definition", 10 ) == 0 )
	/* this is the first atom to be bonded to this flexible bond */
	{
	  if ( number_of_rules >= MAX_NUMBER_OF_FLEX_BOND_RULES )
	    err_panic ( "read_flex_defn",
			"too many flexible-bond rules" );
	  parse_definition_line ( linebuffer, 
				  &bond_rule[number_of_rules].atom[0] );
	  
	  /*
	  print_rule ( &bond_rule[number_of_rules].atom[0] );
	  */
	  /* get the definition for the other atom bonded to this bond */
	  fgets ( linebuffer, sizeof ( linebuffer ), in );
	  /* superscribe `\n` with `\0` */
	  linebuffer[strlen(linebuffer)-1] = '\0';
	  if ( strncmp ( linebuffer, "definition", 10 ) != 0 )
	    err_panic ( "read_flex_defn", 
			"only one atom for flexible-bond rule defined" );
	  parse_definition_line ( linebuffer, 
				  &bond_rule[number_of_rules].atom[1] );
	  /*
	    print_rule ( &bond_rule[number_of_rules].atom[1] );
	    */
	  number_of_rules++;
	}
    }
  fclose ( in );
  return number_of_rules;
}
	  
