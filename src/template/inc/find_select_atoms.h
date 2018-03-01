#ifndef _FIND_SELECT_ATOMS_H
#define _FIND_SELECT_ATOMS_H

extern void find_select_atoms(pdb_atom_pt atoms, int number_of_atoms, pdb_atom_pt b_atoms, int number_of_b_atoms, pdb_atom_pt s_atoms, int *number_of_s_atoms);

extern void create_pdb_line(pdb_atom_pt atom, char line[MAX_LINE_LENGTH], int f);

#endif

