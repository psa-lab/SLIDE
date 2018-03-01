/*
 *
 *    find_hyd_atoms.c     Volker Schnecke     Fri Feb 13 14:42:40 EST 1998
 *
 *    functions:  check_hyd_atom()
 *                find_hyd_atoms()
 *
 *    checks all atoms in the structure to identify Hydrogen-bond
 *    donors and acceptors
 *    
 */

#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "check_rule.h"

/*  This routine checks the atom with index 'index' in the molecule
 *  structure for donor or acceptor properties as described in 'rules'.
 *  The possible return values are DONOR, ACCEPTOR, or NOTHING 
 */
int  check_hyd_atom ( molecule_pt  molecule,
		      hyd_defn_pt  rules,
		      int          index )
{
  atom_pt  atoms;
  rule_pt  rule;
  int      result,
           temp_hyd;
  int      i, j;

  atoms = molecule->atoms;
  temp_hyd = NOTHING;
  /* check acceptor rules first */
  for ( i = 0; 
	i < rules->number_of_acceptor_rules && temp_hyd == NOTHING; 
	i++ )
    {
      rule = &rules->acceptor_rules[i];
      if ( check_atom ( rule->type,
			rule->orbit,
			atoms[index].type,
			atoms[index].orbit ) )
	/* atom matches current acceptor rule atom, so check required and
	   prohibited rules */
	{
	  result = SUCCESS;
	  for ( j = 0; 
		j < MAX_RULES_PER_BONDED_ATOM && result == SUCCESS; 
		j++ )
	    if ( rule->prohibited[j] != NULL )
	      /* the third parameter is the index of a neighbored atom
		 that should not be included in the checks; this is only
		 relevant when checking flexible bond rules, so pass
		 index = -1 here */
	      result = check_prohibited_rule ( molecule,
					       index,
					       -1,
					       rule->prohibited[j] );
	  for ( j = 0; 
		j < MAX_RULES_PER_BONDED_ATOM && result == SUCCESS; 
		j++ )
	    if ( rule->required[j] != NULL )
	      result = check_required_rule ( molecule,
					     index,
					     -1,
					     rule->required[j] );
	  if ( result == SUCCESS )
	    /* all rules fulfilled, i.e. this atom is an ACCEPTOR, but
	       don't return, check first, if this atom can donate, too */
	    temp_hyd = ACCEPTOR;
	}
    }
  /* if current atom is no acceptor, it might still be a donor... */
  for ( i = 0; i < rules->number_of_donor_rules; i++ )
    {
      rule = &rules->donor_rules[i];
      if ( check_atom ( rule->type,
			rule->orbit,
			atoms[index].type,
			atoms[index].orbit ) )
	{
	    result = SUCCESS;
	    for ( j = 0; 
		j < MAX_RULES_PER_BONDED_ATOM && result == SUCCESS; 
		j++ )
	    if ( rule->prohibited[j] != NULL )
	      result = check_prohibited_rule ( molecule,
					       index,
					       -1,
					       rule->prohibited[j] );
	  for ( j = 0; 
		j < MAX_RULES_PER_BONDED_ATOM && result == SUCCESS; 
		j++ )
	    if ( rule->required[j] != NULL )
	      result = check_required_rule ( molecule,
					     index,
					     -1,
					     rule->required[j] );
	  if ( result == SUCCESS )
	    /* hit, so this atom is a donor, but has it already been 
	       identified as an acceptor? */
	    {
	      if ( temp_hyd == ACCEPTOR )
		return DONEPTOR;
	      else
		return DONOR;
	    }
	}
    }
  /* definitely no donor, in that case the routine would have been left
     before, so it might be nothing or a acceptor, temp_hyd has the answer */
  return temp_hyd;
}

/*  This routine checks all atoms in 'molecule' and assigns the 'hyd'
 *  properties to them
 */
void  find_hyd_atoms ( molecule_pt molecule,
		       hyd_defn_pt rules )
{
  int     i;

  for ( i = 0; i < molecule->number_of_atoms; i++ )
    /* only have a closer look at Nitrogens and Oxygens */
    if ( molecule->atoms[i].type == N 
	 || molecule->atoms[i].type == O
	 || molecule->atoms[i].type == S
	 || molecule->atoms[i].type == F 
	 || molecule->atoms[i].type == CL )
      molecule->atoms[i].hyd = check_hyd_atom ( molecule, 
						rules, 
						i );
}
