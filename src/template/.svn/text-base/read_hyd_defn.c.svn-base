#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "types.h"
#include "basics.h"
#include "read_rule.h"

void  read_hyd_defn ( char         *filename,
		      hyd_defn_pt  hyd_rules )
{
  char    linebuffer[MAX_MOL2_LINELENGTH];
  FILE    *in;
  rule_pt acceptor_rules,
          donor_rules;
  int     number_of_donor_rules,
          number_of_acceptor_rules;
  
  in = fopen ( filename, "r" );
  if ( in == NULL ) 
    {
      perror ( "ERROR  " );
      sprintf ( linebuffer, "unable to open file %s", filename );
      err_panic ( "read_hyd_defn", linebuffer );
    }
  acceptor_rules = hyd_rules->acceptor_rules;
  donor_rules = hyd_rules->donor_rules;
  number_of_donor_rules = 0;
  number_of_acceptor_rules = 0;
  while ( fgets ( linebuffer, sizeof ( linebuffer ), in ) != NULL )
    if ( linebuffer[0] != '#' )
      {
	/* superscribe `\n` with `\0` */
	linebuffer[strlen(linebuffer)-1] = '\0';
	if ( strncmp ( linebuffer, "donor", 5 ) == 0 )
	  /* this is the definition of a donor rule */
	  {
	    if ( number_of_donor_rules >= MAX_NUMBER_OF_HYD_RULES )
	      err_panic ( "read_hyd_defn",
			  "too many donor rules" );
	    parse_definition_line ( linebuffer, 
				    &donor_rules[number_of_donor_rules] );
	    /*
	      print_rule ( &donor_rules[number_of_donor_rules] );
	      */
	    number_of_donor_rules++;
	  }
	if ( strncmp ( linebuffer, "acceptor", 8 ) == 0 )
	  /* this is the definition of a acceptor rule */
	  {
	    if ( number_of_acceptor_rules >= MAX_NUMBER_OF_HYD_RULES )
	      err_panic ( "read_hyd_defn",
			  "too many acceptor rules" );
	    parse_definition_line ( linebuffer, 
				    &acceptor_rules[number_of_acceptor_rules]);
	    /*
	    print_rule ( &acceptor_rules[number_of_acceptor_rules] );
	    */
	    number_of_acceptor_rules++;
	  }
    }
  hyd_rules->number_of_donor_rules = number_of_donor_rules;
  hyd_rules->number_of_acceptor_rules = number_of_acceptor_rules;
  fclose ( in );
}
	  
