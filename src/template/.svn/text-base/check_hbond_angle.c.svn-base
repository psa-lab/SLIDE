#include <math.h>
#include <strings.h>
#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "distance.h"
#include "basics.h"

#define STRAIGHT_IN_PLANE    0
#define TWO_POSSIBILITIES    1
#define POINTING_TO_ACCEPTOR 2
#define TETRAHEDRAL          3


void  transform_point ( double             pos[3],
			transform_matrix_t matrix )
{
  double x[3];
  int    i, j;

  for ( i = 0; i < 3; i++ )
    x[i] = pos[i] - matrix[i][3];
  for ( i = 0; i < 3; i++ )
    {
      pos[i] = 0.0;
      for ( j = 0; j < 3; j++ )
	pos[i] += matrix[i][j] * x[j];
    }
  for ( i = 0; i < 3; i++ )
    if ( pos[i] > -0.0005 && pos[i] < 0.0005 )
      pos[i] = 0.0;
}

void  transform_point_back ( double             pos[3],
			     transform_matrix_t matrix )
{
  double x[3];
  int    i, j;

  for ( i = 0; i < 3; i++ )
    x[i] = pos[i];
  for ( i = 0; i < 3; i++ )
    {
      pos[i] = 0.0;
      for ( j = 0; j < 3; j++ )
	pos[i] += (float) matrix[j][i] * x[j];
    }
  for ( i = 0; i < 3; i++ )
    {
      pos[i] += (float) matrix[i][3];
      if ( pos[i] > -0.0005 && pos[i] < 0.0005 )
	pos[i] = 0.0;
    }

}
void  compute_yzx_matrix ( double             angle[3],
			   double             trans[3],
			   transform_matrix_t matrix )
{
  double  c[3];
  double  s[3];
  int     i;

  for ( i = 0; i < 3; i++ )
    {
      s[i] = sin ( angle[i] );
      c[i] = cos ( angle[i] );
      matrix[i][3] = trans[i];
    }
  matrix[0][0] = c[Y] * c[Z];
  matrix[0][1] = s[X] * s[Y] - c[Y] * s[Z] * c[X];
  matrix[0][2] = s[X] * s[Z] * c[Y] + c[X] * s[Y];

  matrix[1][0] = s[Z];
  matrix[1][1] = c[X] * c[Z];
  matrix[1][2] = (-1) * s[X] * c[Z];
  
  matrix[2][0] = (-1) * s[Y] * c[Z];
  matrix[2][1] = c[X] * s[Y] * s[Z] + s[X] * c[Y];
  matrix[2][2] = c[X] * c[Y] - s[X] * s[Y] * s[Z];
}

void  compute_transformation_matrix ( double             n[3],
				      double             ca[3],
				      double             cb[3],
				      transform_matrix_t matrix )
{
  double axis[3],
         angle[3],
         trans[3],
         help[3];
  int    i;
  
  for ( i = 0; i < 3; i++ )
   axis[i] = (double) cb[i] - (double)  ca[i];
  for ( i = 0; i < 3; i++ )
    trans[i] = (double) cb[i];
  if ( axis[Y] == 0 )
    angle[X] = 0;
  else
    {
      angle[X] = (-1) * atan ( axis[Z] / axis[Y] );
      if ( axis[Y] < 0 ) 
	angle[X] -= M_PI;
    }
  if ( axis[Y] * cos ( angle[X] ) - axis[Z] * sin ( angle[X] ) == 0 )
    angle[Z] = 0;
  else
    {
      angle[Z] = atan ( axis[X] / ( axis[Y] * cos ( angle[X] ) 
				    - axis[Z] * sin ( angle[X] ) ) );
      if ( axis[Y] * cos ( angle[X] ) - axis[Z] * sin ( angle[X] ) < 0 ) 
	angle[Z] += M_PI;
    }
  for ( i = 0; i < 3; i++ )
    help[i] = n[i] - trans[i];
  if ( help[X] * cos ( angle[Z] )
       - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
       + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) == 0 )
    angle[Y] = 0;
  else
    {
      angle[Y] = atan ( ( help[Y] * sin ( angle[X] )
			  + help[Z] * cos ( angle[X] ) )
			/ ( help[X] * cos ( angle[Z] )
			  - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
			  + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) ) );
      if ( help[X] * cos ( angle[Z] )
	   - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
	   + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) < 0 )
	angle[Y] += M_PI;
    }
  compute_yzx_matrix ( angle, trans, matrix );
}

/*
 *  This function computes the angle between two lines which
 *  are defined by the tuples (tail,joint) and (joint,head).
 */
double  compute_angle ( double  tail[3],
			double  joint[3],
			double  head[3] )
{
  double vector1[3],
         vector2[3],
         angle;
  int    i;

  for ( i = 0; i < 3; i++ )
    {
      vector1[i] = joint[i] - tail[i];
      vector2[i] = head[i] - joint[i];
    }
  /* formula to compute angle between two vectors
     is |vec1*vec2|/(|vec1|*|vec2|) and |vec1*vec2| is
     vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2] in
     this case, because of the orthogonal base */
  angle = acos ( ( vector1[0] * vector2[0] 
		   + vector1[1] * vector2[1]
		   + vector1[2] * vector2[2] )
		 / sqrt ( vector1[0] * vector1[0]
			  + vector1[1] * vector1[1]
			  + vector1[2] * vector1[2] )
		 / sqrt ( vector2[0] * vector2[0]
			  + vector2[1] * vector2[1]
			  + vector2[2] * vector2[2] ) );
  /* transform to deg */
  angle = 180 - ( angle * 180 / M_PI );
  return angle;
}

void  get_positions ( pdb_atom_pt atoms,
		      int         number_of_atoms,
		      int         index,
		      int         *types,
		      int         number_of_atoms_to_find,
		      double      positions[3][3] )
{
  int  close_atom_indices[10];
  int  number_of_close_atoms,
       found;
  int  i, j, k;

  number_of_close_atoms = 0;
  for ( i = 0; i < number_of_atoms; i++ )
    if ( distance ( atoms[i].pos, atoms[index].pos ) < 2.5 )
      close_atom_indices[number_of_close_atoms++] = i;  
  for ( i = 0; i < number_of_atoms_to_find; i++ )
    {
      positions[i][0] = 9999.9;
      for ( j = 0, found = FALSE; j < number_of_close_atoms && !found; j++ )
	if ( atoms[close_atom_indices[j]].type == types[i] )
	  {
	    for ( k = 0; k < 3; k++ )
	      positions[i][k] = atoms[close_atom_indices[j]].pos[k];
	    found = TRUE;
	  }
      if ( positions[i][0] == 9999.9 )
	{
	  printf ( "close atoms: " );
	  for ( k = 0; k < number_of_close_atoms; k++ )
	    printf ( "%d ", atoms[close_atom_indices[k]].type );
	  printf ( "looking for " );
	  for ( k = 0; k < number_of_atoms_to_find; k++ )
	    printf ( "%d ", types[k] );
	  printf ( "\n" );
	  err_panic ( "get_positions", "neighbored atom not found" );
	}
    }		    
}

double  compute_hydrogen_angle ( pdb_atom_pt  atoms,
				 int          number_of_atoms,
				 int          index,
				 point_pt     acceptor )
{
  transform_matrix_t transformation_matrix;
  double             atom_positions[3][3];
  double             hydrogen_position[3],
                     pos[3];
  double             length,
                     angle;
  int                atom_types[3];
  int                which_case;
  int                i; 

  if ( atoms[index].type == N )
    /* main-chain nitrogen */
    {
      atom_types[0] = CA;
      atom_types[1] = N;
      get_positions ( atoms,
		      number_of_atoms,
		      index,
		      atom_types,
		      2,
		      atom_positions );
      /* move CA to last position in vector points */
      for ( i = 0; i < 3; i++ )
	atom_positions[2][i] = atom_positions[0][i];
      /* search C in previous residue */
      atom_types[0] = C;
      get_positions ( atoms,
		      number_of_atoms,
		      index,
		      atom_types,
		      1,
		      atom_positions );
      length = 1.02;
      which_case = STRAIGHT_IN_PLANE;
    }
  else
    /* side-chain donor */
    switch ( atoms[index].residue_type )
      {
      case ARG:
	switch ( atoms[index].type )
	  {
	  case NE:
	    /* axis is avg(CD,CZ),NE, point is CZ (just one of the
	       three, since everything is planar, hydrogen position
	       if independent of acceptor position) */
	    atom_types[0] = CD;
	    atom_types[1] = NE;
	    atom_types[2] = CZ;
	    get_positions ( atoms,
			    number_of_atoms,
			    index,
			    atom_types,
			    3,
			    atom_positions );
	    length = 1.02;
	    which_case = STRAIGHT_IN_PLANE;
	    break;
	  case NH1:
	    /* axis is CZ, NH1, point is NH2 */
	    atom_types[0] = CZ;
	    atom_types[1] = NH1;
	    atom_types[2] = NH2;
	    get_positions ( atoms,
			    number_of_atoms,
			    index,
			    atom_types,
			    3,
			    atom_positions );
	    length = 1.02;
	    angle = 117;
	    which_case = TWO_POSSIBILITIES;
	    break;
	  case NH2:
	    /* axis is CZ, NH1, point is NH2 */
	    atom_types[0] = CZ;
	    atom_types[1] = NH2;
	    atom_types[2] = NH1;
	    get_positions ( atoms,
			    number_of_atoms,
			    index,
			    atom_types,
			    3,
			    atom_positions );
	    angle = 117;
	    length = 1.02;
	    which_case = TWO_POSSIBILITIES;
	    break;
	  default:
	    break;
	  }
	break;
      case ASN:
	/* axis is CG, ND2, point is OD1 */
	atom_types[0] = CG;
	atom_types[1] = ND2;
	atom_types[2] = OD1;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			3,
			atom_positions );
	angle = 115;
	length = 1.02;
	which_case = TWO_POSSIBILITIES;
	break;
      case GLN:
	/* axis is CD, NE2, point is OE1 */
	atom_types[0] = CD;
	atom_types[1] = NE2;
	atom_types[2] = OE1;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			3,
			atom_positions );
	angle = 115;
	length = 1.02;
	which_case = TWO_POSSIBILITIES;
	break;
      case HIS:
	if ( atoms[index].type == ND1 )
	  /* axis is avg(CG,CE1), ND1, point is CE1 (could be any, since
	     ring is planar) */
	  {
	    atom_types[0] = CG;
	    atom_types[1] = ND1;
	    atom_types[2] = CE1;
	    get_positions ( atoms,
			    number_of_atoms,
			    index,
			    atom_types,
			    3,
			    atom_positions );
	  }
	else
	  /* axis is avg(CD1, CE1), NE2, point is CE1 (could be any, since
	     ring is planar) */
	  {
	    atom_types[0] = CE1;
	    atom_types[1] = NE2;
	    atom_types[2] = CD2;
	    get_positions ( atoms,
			    number_of_atoms,
			    index,
			    atom_types,
			    3,
			    atom_positions );
	  }
	length = 1.02;
	which_case = STRAIGHT_IN_PLANE;
	break;
      case LYS:
	/* axis is CE, NZ, point is acceptor position */
	atom_types[0] = CE;
	atom_types[1] = NZ;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			2,
			atom_positions );
	for ( i = 0; i < 3; i++ )
	  atom_positions[2][i] = acceptor->pos[i];
	angle = 114;
	length = 1.02;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case SER:
	/* axis is CB, OG, point is acceptor position */
	atom_types[0] = CB;
	atom_types[1] = OG;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			2,
			atom_positions );
	for ( i = 0; i < 3; i++ )
	  atom_positions[2][i] = acceptor->pos[i];
	angle = 107;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case THR:
	/* axis is CB, OG1, point is acceptor position */
	atom_types[0] = CB;
	atom_types[1] = OG1;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			2,
			atom_positions );
	for ( i = 0; i < 3; i++ )
	  atom_positions[2][i] = acceptor->pos[i];
	angle = 106;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case TRP:
	/* axis is avg(CD1, CE2), NE1, point is CE2 (could be any, since
	   ring is planar) */
	atom_types[0] = CD1;
	atom_types[1] = NE1;
	atom_types[2] = CE2;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			3,
			atom_positions );
	length = 1.02;
	which_case = STRAIGHT_IN_PLANE;
	break;
      case TYR:
	/* axis is CZ, OH, point is acceptor position */
	atom_types[0] = CZ;
	atom_types[1] = OH;
	get_positions ( atoms,
			number_of_atoms,
			index,
			atom_types,
			2,
			atom_positions );
	for ( i = 0; i < 3; i++ )
	  atom_positions[2][i] = acceptor->pos[i];
	angle = 108;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      default:
	err_panic ( "compute_hydrogen_position",
		    "unknown residue" );
	break;
      }
  compute_transformation_matrix ( atom_positions[2],          /* X-axis */
				  atom_positions[0],          /* Y-axis */
				  atom_positions[1],          /* origin */
				  transformation_matrix );
  switch ( which_case )
    {
    case POINTING_TO_ACCEPTOR:
      /* acceptor position is in the positive X-part of the XY-plane */
      hydrogen_position[X] = length * sin ( ( 180 - angle ) / 180 * M_PI );
      hydrogen_position[Y] = length * cos ( ( 180 - angle ) / 180 * M_PI );
      hydrogen_position[Z] = 0;
      break;
    case TWO_POSSIBILITIES:
      /* two positions for the hydrogen are possible, both on the XY-plane
	 of the transformed system, one with positive X ordinate, the other
	 is it's mirror image with regard to the Y-axis */
      /* better don't mess around with the original acceptor entry... */
      for ( i = 0; i < 3; i++ )
	pos[i] = acceptor->pos[i];
      transform_point ( pos,
			transformation_matrix );
      /* default is 'right' hydrogen position */
      hydrogen_position[X] 
	= length * sin ( ( 180 - angle ) * M_PI / 180 );
      hydrogen_position[Y] 
	= length * cos ( ( 180 - angle ) * M_PI / 180 );
      if ( pos[X] < 0 )
	/* projection of accceptor is in XY-plane left of Y-axis, so
	   choose 'left' hydrogen position */
	hydrogen_position[X] *= (-1);
      hydrogen_position[Z] = 0;
      break;
    case STRAIGHT_IN_PLANE:
      angle = compute_angle ( atom_positions[0],
			      atom_positions[1],
			      atom_positions[2] );
      angle = ( 360 - angle ) / 2;
      /* Hydrogen is located on the XY-plane just between atom_positions[0]
	 and atom_positions[2], atom_positions[1] is the origin of the transformed 
	 system */	 
      hydrogen_position[X] = length * sin ( ( 180 - angle ) / 180 * M_PI );
      hydrogen_position[Y] = length * cos ( ( 180 - angle ) / 180 * M_PI );
      hydrogen_position[Z] = 0;
      for ( i = 0; i < 3; i++ )
	pos[i] = atom_positions[2][i];
      transform_point ( pos,
			transformation_matrix );
      if ( pos[X] > 0 )
	/* projection of the third point is in XY-plane with a positive
	   X-ordinate, thus the hydrogen needs a position with negative
	   X-ordinate */
	hydrogen_position[X] *= (-1);
      break;
    default:
      break;
    }
  /* transform hydrogen position back to world coordinate system */
  transform_point_back ( hydrogen_position, 
			 transformation_matrix );      
  angle = compute_angle ( atoms[index].pos,
			  hydrogen_position,
			  acceptor->pos );
  return angle;
}
				 

int  check_hbond_angle ( pdb_atom_pt  atoms,
			 int          number_of_atoms,
			 int          atom_index,
			 point_pt     point )
{
  double  angle;

  angle = 
    compute_hydrogen_angle ( atoms,
			     number_of_atoms,
			     atom_index,
			     point );
  if ( angle > MIN_HYDROGEN_ANGLE && angle < MAX_HYDROGEN_ANGLE )
    return TRUE;
  else
    return FALSE;
}


