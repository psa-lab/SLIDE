#ifndef _READ_FLEX_DEFN_H
#define _READ_FLEX_DEFN_H

extern void  parse_definition_line ( char    *line,
				     rule_pt rule );

extern int  read_flex_defn ( char              *filename,
			     flex_bond_defn_pt bond_rule );

#endif
