#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "defs.h"
#include "types.h"
#include <basics.h>
#include <err_handle.h>

#define column(n)       (linebuffer+(n)-1)

int  read_pdb_file ( char  *filename, pdb_atom_pt pdb_atoms )
{
  FILE         *in;
  pdb_atom_pt  atom;
  char         linebuffer[MAX_LINE_LENGTH];
  int          number_of_atoms;
  char c13, c14;

  in = slide_fopen ( filename, "r" );
  if ( in == NULL ) 
    {
      sprintf ( linebuffer, "unable to open file %s", filename );
      err_panic ( "read_pdb_file", linebuffer );
    }
  atom = pdb_atoms;
  number_of_atoms = 0;
  while ( fgets ( linebuffer, MAX_LINE_LENGTH, in ) != NULL )
    if ( strncmp ( linebuffer, "ATOM", 4 ) == 0 )
      {
	atom->type = UNKNOWN;
	strncpy ( atom->atom_number, &linebuffer[6], 5 );
	/*	sscanf ( &linebuffer[13], "%4s", atom->atom_name );*/
	strncpy ( atom->atom_name, &linebuffer[12], 4 );
	strncpy ( atom->residue_number, &linebuffer[22], 4 );
	strncpy ( atom->residue_name, &linebuffer[17], 3 );
	atom->pos[X] = atof ( &linebuffer[30] );
	atom->pos[Y] = atof ( &linebuffer[38] );
	atom->pos[Z] = atof ( &linebuffer[46] );
	c13 = atom->atom_name[0];
	c14 = toupper (atom->atom_name[1]);
	if (( c13 == ' ' || isdigit(c13)) && (c14 == 'H' || c14 == 'D'))
	  {
	    fprintf (stderr, "Error: Hydrogens in PDB file.\n"
		     "%s\n"
		     "Please remove hydrogens using pdbdehydrogen, or use the script temp_gen to generate the template\n", linebuffer);
	    exit(1);
	  }
	if ( ( strncmp ( atom->atom_name, " N", 2 ) == 0 
	       && strncmp ( atom->residue_name, "PRO", 3 ) != 0 )
	     || strncmp ( atom->atom_name, " O", 2 ) == 0 ) 
	  { 
	    if ( strncmp ( atom->atom_name, " N", 2 ) == 0 ) 
	      atom->act = DONOR;
	    if ( strncmp ( atom->atom_name, " O", 2 ) == 0 ) 
	      atom->act = ACCEPTOR;
	    if ( strncmp ( atom->atom_name, " OG", 3 ) == 0 
		 || strncmp ( atom->atom_name, " OH", 3 ) == 0 ) 
	      atom->act = DONEPTOR;
	  }
	else
	  atom->act = NOTHING;
	atom++;
	number_of_atoms++;
	if ( number_of_atoms == MAX_NUMBER_OF_ATOMS )
	  err_panic ( "read_pdb_file", 
		      "more than MAX_NUMBER_OF_ATOMS atoms" );
     }
  else
    if ( strncmp ( linebuffer, "HETATM", 6 ) == 0
	 && strncmp ( "HOH", (char *) linebuffer + 17, 3 ) != 0  )
      {
	atom->type = HETERO;
	strncpy ( atom->atom_number, &linebuffer[6], 5 );
	/*	sscanf ( &linebuffer[12], "%4s", atom->atom_name );*/
	strncpy ( atom->atom_name, &linebuffer[12], 4 );
	strncpy ( atom->residue_number, &linebuffer[22], 4 );
	strncpy ( atom->residue_name, &linebuffer[17], 3 );
	atom->pos[X] = atof ( &linebuffer[30] );
	atom->pos[Y] = atof ( &linebuffer[38] );
	atom->pos[Z] = atof ( &linebuffer[46] );
	atom->act = NOTHING;
	atom->hphil = 317;

	      /*****************************************************************
		handle the metal HETATM to count the metal hbond in scoring

                   metal_class  metal_name                    metal_hbond_distance
                    1            Ca, Na, K                       2.9 angstroms
                    2            Co, Cu, Fe, Mg, Mn, Ni, Zn      2.6 angstroms

                		       -- added by Litian He  12/02/2003  
	      ******************************************************************/

	      /* "CA" is a special case, both atom name and residue name should be "CA"  */
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'A' && 
		   *column ( 15 ) == ' ' && *column ( 18 ) == ' ' && 
		   *column ( 19 ) == 'C' && *column ( 20 ) == 'A' )
		 {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = CA_M;
		      atom->act = METAL_1;
		      atom->hphil = 635;
		      atom->rad = RAD_CA;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      /* we only need to check the atom name for all other metals */
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'O' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = CO;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_CO;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'U' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = CU;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_CU;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'F' && *column ( 14 ) == 'E' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = FE;
		      atom->act = METAL_2;
		      atom->hphil = 317;  /*  Sameer - is this deliberate or oversight ???*/
		      atom->hphil = 635;
		      atom->rad = RAD_FE;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'K' && *column ( 14 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = K;
		      atom->act = METAL_1;
		      atom->hphil = 635;
		      atom->rad = RAD_K;
		      
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'M' && *column ( 14 ) == 'G' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = MG;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_MG;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'M' && *column ( 14 ) == 'N' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = MN;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_MN;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'N' && *column ( 14 ) == 'A' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = NA;
		      atom->act = METAL_1;
		      atom->hphil = 635;
		      atom->rad = RAD_NA;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'N' && *column ( 14 ) == 'I' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = NI;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_NI;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }
	      if ( *column ( 13 ) == 'Z' && *column ( 14 ) == 'N' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atom->atom_name );
		      atom->type = ZN;
		      atom->act = METAL_2;
		      atom->hphil = 635;
		      atom->rad = RAD_ZN;
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atom->atom_name, linebuffer);
#endif
		  }

	      /******************* end of metal check ********************/


	if ( strncmp ( atom->atom_name, " N", 2 ) == 0 
	     || strncmp ( atom->atom_name, " O", 2 ) == 0 ) 
	  {
	    if ( strncmp ( atom->atom_name + 1, " DD", 3 ) == 0 )
	      atom->act = DONOR;
	    if ( strncmp ( atom->atom_name + 1, " 1D", 3 ) == 0 )
	      atom->act = DONOR;
	    if ( strncmp ( atom->atom_name + 1, " 2D", 3 ) == 0 )
	      atom->act = DONOR;
	    if ( strncmp ( atom->atom_name + 1, " AA", 3 ) == 0 )
	      atom->act = ACCEPTOR;
	    if ( strncmp ( atom->atom_name + 1, " NN", 3 ) == 0 )
	      atom->act = DONEPTOR;
	  }
	if ( strncmp ( atom->atom_name, " N", 2 ) == 0 )
	  atom->hphil = 350;
	if ( strncmp ( atom->atom_name, " O", 2 ) == 0 ) 
	  atom->hphil = 530;
	if ( strncmp ( atom->atom_name, " C", 2 ) == 0 ) 
	  atom->hphil = 80;
	atom++;
	number_of_atoms++;
	if ( number_of_atoms == MAX_NUMBER_OF_ATOMS )
	  err_panic ( "read_pdb_file", 
		      "more than MAX_NUMBER_OF_ATOMS atoms" );
     }

  fclose ( in );
  return number_of_atoms;
}
