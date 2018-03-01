#ifndef _READ_MOL2_H
#define _READ_MOL2_H

extern void  read_mol2 ( char        *filename,
			 molecule_pt molecule );

extern void  assign_type_and_orbit ( char  *str,
				     int   *type,
				     int   *orbit );
#endif
