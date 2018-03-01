#include <string.h>
#include "defs.h"
#include "types.h"
#include "hydro.h"
#include "basics.h"

void  assign_atom_and_residue_types ( pdb_atom_pt atom )
{
  switch ( atom->residue_name[0] )
    {
    case 'A':
      switch ( atom->residue_name[1] )
	{
	case 'C':
	  atom->residue_type = UNKNOWN;
	  break;
	case 'L':
	  atom->residue_type = ALA;
	  break;
	case 'R':
	  atom->residue_type = ARG;
	  break;
	case 'S':
	  switch ( atom->residue_name[2] )
	    {
	    case 'N':
	      atom->residue_type = ASN;
	      break;
	    case 'P':
	      atom->residue_type = ASP;
	      break;
	    default:
	      err_unknown_residue ( atom->residue_name, 
				    atom->residue_number );
	      break;
	    }
	  break;
	default:
	  err_unknown_residue ( atom->residue_name, 
				atom->residue_number );
	  break;
	}
      break;
    case 'C':
      if ( strncmp ( atom->residue_name, "CYS", 3 ) == 0 )
	atom->residue_type = CYS;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'G':
      if ( atom->residue_name[1] == 'L' )
	switch ( atom->residue_name[2] )
	  {
	  case 'N':
	      atom->residue_type = GLN;
	      break;
	  case 'U':
	    atom->residue_type = GLU;
	      break;
	  case 'Y':
	    atom->residue_type = GLY;
	      break;
	  default:
	    err_unknown_residue ( atom->residue_name, 
				  atom->residue_number );
	    break;
	  }
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'H':
      if ( strncmp ( atom->residue_name, "HIS", 3 ) == 0 )
	atom->residue_type = HIS;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'I':
      if ( strncmp ( atom->residue_name, "ILE", 3 ) == 0 )
	atom->residue_type = ILE;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'L':
      switch ( atom->residue_name[1] )
	{
	case 'E':
	  atom->residue_type = LEU;
	  break;
	case 'Y':
	  atom->residue_type = LYS;
	  break;
	default:
	  err_unknown_residue ( atom->residue_name, 
				atom->residue_number );
	  break;
	}
      break;
    case 'M':
      if ( strncmp ( atom->residue_name, "MET", 3 ) == 0 )
	atom->residue_type = MET;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'P':
      switch ( atom->residue_name[1] )
	{
	case 'H':
	  atom->residue_type = PHE;
	  break;
	case 'R':
	case 'C':   /* treat PCA (pyrrolidone carboxylic acid) like Proline */
	  atom->residue_type = PRO;
	  break;
	default:
	  err_unknown_residue ( atom->residue_name, 
				atom->residue_number );
	  break;
	}
      break;
    case 'S':
      if ( strncmp ( atom->residue_name, "SER", 3 ) == 0 )
	atom->residue_type = SER;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    case 'T':
      switch ( atom->residue_name[1] )
	{
	case 'H':
	  atom->residue_type = THR;
	  break;
	case 'R':
	  atom->residue_type = TRP;
	  break;
	case 'Y':
	  atom->residue_type = TYR;
	  break;
	default:
	  err_unknown_residue ( atom->residue_name, 
				atom->residue_number );
	  break;
	}
      break;
    case 'V':
      if ( strncmp ( atom->residue_name, "VAL", 3 ) == 0 )
	atom->residue_type = VAL;
      else
	err_unknown_residue ( atom->residue_name, 
			      atom->residue_number );
      break;
    default:
      err_unknown_residue ( atom->residue_name, 
			    atom->residue_number );
      break;
    }
  switch ( atom->atom_name[1] )
    {
    case 'A':
      switch ( atom->atom_name[2] )
	{
	case 'D':
	  if ( strncmp ( atom->residue_name, "ASP", 3 ) )
	    switch ( atom->atom_name[3] )
	      {
	      case '1':
		atom->type = OD1;
		break;
	      case '2':
		atom->type = ND2;
		break;
	      default:
		atom->type = UNKNOWN;
		err_unknown_atom ( atom->atom_name, 
				   atom->residue_name, 
				   atom->residue_number );
		break;
	      }
	  else
	    if ( strcmp ( atom->residue_name, "HIS" ) )
	      switch ( atom->atom_name[3] )
		{
		case '1':
 		  atom->type = ND1;
		  break;
		case '2':
		  atom->type = CD2;
		  break;
		default:
		  atom->type = UNKNOWN;
		  err_unknown_atom ( atom->atom_name, 
				     atom->residue_name, 
				     atom->residue_number );
		  break;
	      }
	    else
	      {
		atom->type = UNKNOWN;
		err_unknown_atom ( atom->atom_name, 
				   atom->residue_name, 
				   atom->residue_number );
	      }
	  break;
	case 'E':
	  if ( strcmp ( atom->residue_name, "GLN" ) )
	    switch ( atom->atom_name[3] )
	      {
	      case '1':
		atom->type = OE1;
		break;
	      case '2':
		atom->type = NE2;
		break;
	      default:
		atom->type = UNKNOWN;	
		err_unknown_atom ( atom->atom_name, 
				   atom->residue_name, 
				   atom->residue_number );
		break;
	      }
	  else
	    if ( strcmp ( atom->residue_name, "HIS" ) )
	      switch ( atom->atom_name[3] )
		{
		case '1':
		  atom->type = CE1;
		  break;
		case '2':
		  atom->type = NE2;
		  break;
		default:
		  atom->type = UNKNOWN;
		  err_unknown_atom ( atom->atom_name, 
				     atom->residue_name, 
				     atom->residue_number );
		  break;
	      }
	    else
	      {
		atom->type = UNKNOWN;
		err_unknown_atom ( atom->atom_name, 
				   atom->residue_name, 
				   atom->residue_number );
	      }
	  break;
	default:
	  atom->type = UNKNOWN;
	  err_unknown_atom ( atom->atom_name, 
			     atom->residue_name, 
			     atom->residue_number );
	  break;
	}
      break;
    case 'C':
      switch ( atom->atom_name[2] )
	{
	case ' ':
	  atom->type = C;
	  break;
	case 'A':
	  atom->type = CA;
	  break;
	case 'B':
	  atom->type = CB;
	  break;
	case 'D':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = CD;
	      break;
	    case '1':
	      atom->type = CD1;
	      break;
	    case '2':
	      atom->type = CD2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'E':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = CE;
	      break;
	    case '1':
	      atom->type = CE1;
	      break;
	    case '2':
	      atom->type = CE2;
	      break;
	    case '3':
	      atom->type = CE3;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'G':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = CG;
	      break;
	    case '1':
	      atom->type = CG1;
	      break;
	    case '2':
	      atom->type = CG2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'H':
	  switch ( atom->atom_name[3] )
	    {
	    case '2':
	    case '3': 
	      atom->type = CH2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'Z':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = CZ;
	      break;
	    case '2':
	      atom->type = CZ2;
	      break;
	    case '3':
	      atom->type = CZ3;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	default:
	  atom->type = UNKNOWN;
	  break;
	}
      break;
    case 'O':
      switch ( atom->atom_name[2] )
	{
	case ' ':
	  atom->type = O;
	  break;
	case 'D':
	  switch ( atom->atom_name[3] )
	    {
	    case '1':
	      atom->type = OD1;
	      break;
	    case '2':
	      atom->type = OD2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'E':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':        /* atom OE exists in residue PCA */
	    case '1':
	      atom->type = OE1;
	      break;
	    case '2':
	      atom->type = OE2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'G':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = OG;
	      break;
	    case '1':
	      atom->type = OG1;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'H':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = OH;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'X':
	  if ( atom->atom_name[3] == 'T' )
	    atom->type = OXT;
	  else
	    {
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	    }
	  break;
	case 'T':
	  atom->type = OXT;
	  break;
	default:
	  atom->type = UNKNOWN;
	  err_unknown_atom ( atom->atom_name, 
			     atom->residue_name, 
			     atom->residue_number );
	  break;
	}
      break;
    case 'N':
      switch ( atom->atom_name[2] )
	{
	case ' ':
	  atom->type = N;
	  break;
	case 'D':
	  switch ( atom->atom_name[3] )
	    {
	    case '1':
	      atom->type = ND1;
	      break;
	    case '2':
	      atom->type = ND2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'E':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = NE;
	      break;
	    case '1':
	      atom->type = NE1;
	      break;
	    case '2':
	      atom->type = NE2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'H':
	  switch ( atom->atom_name[3] )
	    {
	    case '1':
	      atom->type = NH1;
	      break;
	    case '2':
	      atom->type = NH2;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	case 'Z':
	  switch ( atom->atom_name[3] )
	    {
	    case ' ':
	      atom->type = NZ;
	      break;
	    default:
	      atom->type = UNKNOWN;
	      err_unknown_atom ( atom->atom_name, 
				 atom->residue_name, 
				 atom->residue_number );
	      break;
	    }
	  break;
	default:
	  atom->type = UNKNOWN;
	  err_unknown_atom ( atom->atom_name, 
			     atom->residue_name, 
			     atom->residue_number );
	  break;
	}
      break;
    case 'S':
      switch ( atom->atom_name[2] )
	{
	case 'D':
	  atom->type = SD;
	  break;
	case 'G':
	  atom->type = SG;
	  break;
	default:
	  atom->type = UNKNOWN;
	  err_unknown_atom ( atom->atom_name, 
			     atom->residue_name, 
			     atom->residue_number );
	  break;
	}
      break;
    default:
      atom->type = UNKNOWN;
      err_unknown_atom ( atom->atom_name, 
			 atom->residue_name, 
			 atom->residue_number );
      break;
    }
}

void  assign_hphil ( pdb_atom_pt  atoms,
		     int          number_of_atoms )
{
  int  i;
  
  for ( i = 0; i < number_of_atoms; i++ )
    {
      if ( (atoms[i].type != HETERO ) && ( atoms[i].act != METAL_1 ) && ( atoms[i].act != METAL_2 ))
	{
	  assign_atom_and_residue_types ( &atoms[i] );
	  /*	  atoms[i].hphil = hydro[atoms[i].residue_type][atoms[i].type];*/
	}
    }
}
