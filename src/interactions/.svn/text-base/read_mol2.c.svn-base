#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"
#include "assign_type.h"

void  read_mol2 ( char        *filename,
		  molecule_pt molecule )
{
  char  linebuffer[MAX_MOL2_LINELENGTH];
  char  linebuffer2[MAX_MOL2_LINELENGTH];
  FILE  *in;
  int   i;
  int   rv;
  char msg[1024];
  char *rv2;

  molecule->number_of_atoms =0; 
  molecule->number_of_bonds =0;
  molecule->name[0] = 0;  
  in = fopen ( filename, "r" );
  if ( in == NULL ) 
    {
      sprintf ( linebuffer, "unable to open file %s", filename );
      err_panic ( "read_mol2", linebuffer );
    }
  do
    fgets ( linebuffer, sizeof ( linebuffer ), in );
  while ( linebuffer[0] != '@' );
  do
    {
      /*printf ("x%s",linebuffer);*/
      if ( strncmp ( linebuffer, "@<TRIPOS>MOLECULE", 17 ) == 0 )
	{
	  /* read molecule name and number of atoms and bonds */
	  rv2 = fgets ( linebuffer, sizeof ( linebuffer ), in );
	  /*	  printf ("return value %d", rv);*/
	  if (rv2 == NULL){
	    sprintf(msg, "Could not read the molecule name (%s)", filename);
            err_panic("read_mol2", msg);
	  }
	  
	  /* printf ("%s",linebuffer);*/
	  rv = sscanf(linebuffer, "%s", molecule->name );
	  if (rv < 1){
	    sprintf(msg, "Could not read the molecule name (%s)", filename);
            err_panic("read_mol2", msg);
	  }
	  rv = fscanf(in, "%d %d %*d", 
		   &molecule->number_of_atoms, 
		   &molecule->number_of_bonds );
          if(rv < 2){
            sprintf(msg, "Could not read the number of atoms and number of "
                    "bonds (%s) Molecule: %s", filename, molecule->name);
            err_panic("read_mol2", msg);
          }
	  /*printf ("%d %d\n",molecule->number_of_atoms,
	    molecule->number_of_bonds);*/
	  if ( molecule->number_of_atoms > MAX_NUMBER_OF_MOL2_ATOMS )
	    {
	      sprintf(msg, "more than MAX_NUMBER_OF_MOL2_ATOMS in "
                      "%s\n\tMolecule: %s", filename, molecule->name);
	      err_panic ( "read_mol2", msg);
	    }
	  if ( molecule->number_of_bonds > MAX_NUMBER_OF_MOL2_BONDS ){
            sprintf(msg, "more than MAX_NUMBER_OF_MOL2_BONDS in "
                    "%s\n\tMolecule: %s", filename, molecule->name);
	    err_panic ( "read_mol2", msg);
	  }
	}
      if ( strncmp ( linebuffer, "@<TRIPOS>ATOM", 13 ) == 0 )
	{ 
	    /* read all atom data */
	for ( i = 0; i < molecule->number_of_atoms; i++ )
	  {
	    /* we cannot use fscanf() here, since we only want to grep
	       six items out of a line which can have six or more entries */
	    fgets ( linebuffer, sizeof ( linebuffer ), in );
	    /*printf ("%s",linebuffer);*/
	    rv = sscanf(linebuffer, 
		     "%d %s %f %f %f %s %d %s %f",
		     &molecule->atoms[i].number,
		     molecule->atoms[i].name,
		     &molecule->atoms[i].pos[X],
		     &molecule->atoms[i].pos[Y],
		     &molecule->atoms[i].pos[Z],
		     molecule->atoms[i].type_str,
 		     &molecule->atoms[i].substrid,
		     molecule->atoms[i].substrname,
 		     &molecule->atoms[i].charge );
            if(rv < 6){
              sprintf(msg, "Error reading atom line\n\t(%s)\n\tFile: %s\n"
                      "\tMolecule: %s\n", linebuffer, filename, molecule->name);
              err_panic("read_mol2", msg);
            }
	    molecule->atoms[i].hyd = NOTHING;
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
	}
      if ( strncmp ( linebuffer, "@<TRIPOS>BOND", 13 ) == 0 )
	/* read all bond data */
	for ( i = 0; i < molecule->number_of_bonds; i++ )
	  {
	    fgets ( linebuffer, sizeof ( linebuffer ), in );
	    rv = sscanf ( linebuffer,
		     "%d %d %d %s",
		     &molecule->bonds[i].number,
		     &molecule->bonds[i].atom1,
		     &molecule->bonds[i].atom2,
		     molecule->bonds[i].type_str );
            if(rv < 4){
              sprintf(msg, "Error reading bond line\n\t(%s)\n\tFile: %s\n"
                      "\tMolecule: %s\n", linebuffer, filename, molecule->name);
              err_panic("read_mol2", msg);
            }
             
	    /*fscanf ( in,
                     "%d %d %d %s %*s",
                     &molecule->bonds[i].number,
                     &molecule->bonds[i].atom1,
                     &molecule->bonds[i].atom2,
                     molecule->bonds[i].type_str );*/
	    
	    /*printf ("%s",linebuffer);*/
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
			{
			  if ( molecule->bonds[i].type_str[1] == 'r' 
			       || molecule->bonds[i].type_str[1] == 'R' )
			    molecule->bonds[i].type = AROMATIC;
			  else if ( molecule->bonds[i].type_str[1] == 'm' 
			       || molecule->bonds[i].type_str[1] == 'M' )
			    molecule->bonds[i].type = AMIDE;
			}
		      else
		      if ( molecule->bonds[i].type_str[0] == 'u' 
			   || molecule->bonds[i].type_str[0] == 'U' )
			/* it is  an unknown bond type, set the bond type
			   to DOUBLE, so that this bond will not be rotated */
			molecule->bonds[i].type = DOUBLE;
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
  molecule->number_of_added_hydrogens = 0;
}
