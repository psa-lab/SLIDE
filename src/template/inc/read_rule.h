#ifndef _READ_RULE_H
#define _READ_RULE_H

extern char  *get_rule_item ( char         *line, 
			      rule_list_pt rule );
  
extern void  parse_definition_line ( char    *line,
				     rule_pt rule );

extern void  print_rule_item ( rule_list_pt rule );

extern void  print_rule ( rule_pt  rule );

#endif
