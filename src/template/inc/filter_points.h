#ifndef  _FILTER_POINTS_H
#define  _FILTER_POINTS_H

extern void  filter_points ( point_pt     points,
			     int          number_of_points,
			     pdb_atom_pt  atoms,
			     int          number_of_atoms );

extern void  filter_ligand_points ( point_pt    points,
				    int         number_of_points,
				    molecule_pt ligand,
				    double      threshold );

#endif
