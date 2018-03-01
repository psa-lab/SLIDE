/** find_select_atoms looks for possible protein atoms that could hydrogen bond
to the binding_site atoms. It also checks the neighbors of the selected atoms
if they fall outside the binding site **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

#define THRESHOLD1 1.8  /** max single bond legnth **/
#define THRESHOLD2 2.6   /** max distance for atoms two bond lengths apart **/
#define THRESHOLD3 0.1
#define THRESHOLD4 3.0   /** max distance for hetatms 2 bonds apart **/

/** creates a valid pdb file string from pdb_atom_pt data structure **/
void create_pdb_line(pdb_atom_pt atom, char line[MAX_LINE_LENGTH], int f)
{
	if(f==0)
		sprintf(line, "ATOM  %5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f", atom->atom_number, atom->atom_name, atom->residue_name,
		    atom->chainID, atom->residue_number, atom->icode, atom->pos[0], atom->pos[1], atom->pos[2]);
	else if(f==1)
		sprintf(line, "HETATM%5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f", atom->atom_number, atom->atom_name, atom->residue_name,
		    atom->chainID, atom->residue_number, atom->icode, atom->pos[0], atom->pos[1], atom->pos[2]);
}

/** add an atom to a list of atoms **/
void insert_list(pdb_atom_pt b, pdb_atom_pt s)
{
	strcpy(s->atom_name, b->atom_name);
	strcpy(s->atom_number, b->atom_number);
	strcpy(s->residue_name, b->residue_name);
	strcpy(s->residue_number, b->residue_number);
	s->icode=b->icode;
	s->chainID=b->chainID;
	s->pos[0]=b->pos[0];
	s->pos[1]=b->pos[1];
	s->pos[2]=b->pos[2];
	s->type=b->type;
	s->residue_type=b->residue_type;
	s->hphil=b->hphil;
	s->act=b->act;
	s->class=b->class;
	s->rad=b->rad;
}

int get_neighbors (char* atom_name, char* residue_name, char chain_id, int residue_number, char* match_A, char* match_C, pdb_atom_pt atoms, int number_of_atoms, pdb_atom_pt s_atoms, int *number_of_s_atoms )
{
	int j, control;
	int Anum, Cnum;
	int residue_number_A;
	char pdbline_A[MAX_LINE_LENGTH];
	char atom_A_line[MAX_LINE_LENGTH];
	char atom_C_line[MAX_LINE_LENGTH];

	char atom_name_A[5], residue_name_A[5], residue_insertion_A;
	char chain_id_A;

	float x_A, y_A, z_A;

	control = 0;
	strncpy(atom_A_line,"NOT FOUND ", 10);
	strncpy(atom_C_line,"NOT FOUND ", 10);

	for(j=0; j<number_of_atoms; j++)
	{
		if(atoms[j].type==UNKNOWN) /** not hetero **/
		{
			strcpy(atom_name_A, atoms[j].atom_name);
			strcpy(residue_name_A, atoms[j].residue_name);
			chain_id_A = atoms[j].chainID;
			residue_number_A = atoi(atoms[j].residue_number);
			residue_insertion_A=atoms[j].icode;
			x_A=atoms[j].pos[0];
			y_A=atoms[j].pos[1];
			z_A=atoms[j].pos[2];

			create_pdb_line(&atoms[j], pdbline_A, 0);

			if ((chain_id_A == chain_id)&&
			    (residue_number_A == residue_number)&&
			    (!strncmp(residue_name_A, residue_name, 4))&&
			    (!strncmp(atom_name_A, match_A, 4)))
			{
				strncpy(atom_A_line, pdbline_A, 54 );
				atom_A_line[54]='\0';
				Anum=j;
			}

			if ((chain_id_A == chain_id)&&
			    (residue_number_A == (residue_number))&&
			    (!strncmp(residue_name_A, residue_name, 3))&&
			    (!strncmp(atom_name_A, match_C, 4)))
			{
				strncpy(atom_C_line, pdbline_A, 54 );
				atom_C_line[54]='\0';
				Cnum=j;
			}
		}
	}

	if((strncmp(atom_A_line, "NOT FOUND ", 10))&&
	    (strncmp(atom_C_line, "NOT FOUND ", 10))){
		insert_list(&atoms[Anum], &s_atoms[(*number_of_s_atoms)]);
		(*number_of_s_atoms)++;
		insert_list(&atoms[Cnum], &s_atoms[(*number_of_s_atoms)]);
		(*number_of_s_atoms)++;
	}
	else control = 1;

	return control;
}

int get_hetero_neighbors (char* atom_name, float x, float y, float z, pdb_atom_pt atoms, int number_of_atoms, pdb_atom_pt s_atoms, int *number_of_s_atoms) 
{
	int j, k, n, control;
	float x_A, y_A, z_A;
	float x_atom_A, y_atom_A, z_atom_A;
	float x_C, y_C, z_C;
	double distance, distance_A_C, distance_B_C;
	char pdbline_A[MAX_LINE_LENGTH];
	char pdbline_C[MAX_LINE_LENGTH];
	char atom_A_line[MAX_LINE_LENGTH];
	char atom_C_line[MAX_LINE_LENGTH];
	char atom_name_A[5], residue_name_A[5];
	char atom_name_C[5], residue_name_C[5];
	int Anum=-1, Cnum=-1;

	k=0;
	control=0;
	strncpy(atom_A_line,"NOT FOUND ", 10);
	strncpy(atom_C_line,"NOT FOUND ", 10);

	for(j=0; j<number_of_atoms; j++)
	{
		if(atoms[j].type==HETERO || atoms[j].act == METAL_1 || atoms[j].act == METAL_2)
		{
			strcpy(atom_name_A, atoms[j].atom_name);
			strcpy(residue_name_A, atoms[j].residue_name); /** not water **/
			x_A=atoms[j].pos[0];
			y_A=atoms[j].pos[1];
			z_A=atoms[j].pos[2];

			create_pdb_line(&atoms[j], pdbline_A, 1);
			if(k<2)
			{

				distance = sqrt(((x_A - x)*(x_A - x)) + 
				    ((y_A - y)*(y_A - y)) + 
				    ((z_A - z)*(z_A - z)));

				if ((distance < THRESHOLD1)&&(distance > THRESHOLD3) && (k==0))
				{
					strncpy(atom_A_line, pdbline_A, 54 );
					atom_A_line[54]='\0';
					x_atom_A = x_A;
					y_atom_A = y_A;
					z_atom_A = z_A;
					Anum=j;
					k++;
				}
				else if ((distance < THRESHOLD1)&&(distance > THRESHOLD3) && (k==1) && (strncmp(atom_A_line,
				    pdbline_A, 54)))
				{
					strncpy(atom_C_line, pdbline_A, 54 );
					Cnum=j;
					k++;
				}
			}
			else if (k >= 2)
				break;
		}
	}

	if (k<2){                    /* only one neighbor was found for atom B */
		/* now find a neighbor of A */
		n = 0;
		for(j=0; j<number_of_atoms; j++)
		{
			if(atoms[j].type==HETERO || atoms[j].act == METAL_1 || atoms[j].act == METAL_2 )
			{
				strcpy(atom_name_C, atoms[j].atom_name);
				strcpy(residue_name_C, atoms[j].residue_name); /** not water **/
				x_C=atoms[j].pos[0];
				y_C=atoms[j].pos[1];
				z_C=atoms[j].pos[2];

				create_pdb_line(&atoms[j], pdbline_C, 1);

				distance_A_C = sqrt(((x_atom_A - x_C)*(x_atom_A - x_C)) + 
				    ((y_atom_A - y_C)*(y_atom_A - y_C)) + 
				    ((z_atom_A - z_C)*(z_atom_A - z_C)));

				distance_B_C = sqrt(((x - x_C)*(x - x_C)) + 
				    ((y - y_C)*(y - y_C)) + 
				    ((z - z_C)*(z - z_C)));

				if ((distance_A_C < THRESHOLD1)&&(distance_A_C > THRESHOLD3)&&
				    (distance_B_C < THRESHOLD4)&&(distance_B_C > THRESHOLD3) && (n<1) )
				{
					/**insert_list(&atoms[i], &s_atoms[(*number_of_s_atoms)]);	
					          (*number_of_s_atoms)++;**/
					Cnum=j;
					strncpy(atom_C_line, pdbline_C, 54 );
					n++;
				}
			}
			/*if (n == 1)
			         break;*/
		}
	}
	if((strncmp(atom_A_line, "NOT FOUND ", 10))&&
	    (strncmp(atom_C_line, "NOT FOUND ", 10))){

		insert_list(&atoms[Anum], &s_atoms[(*number_of_s_atoms)]);
		(*number_of_s_atoms)++;
		insert_list(&atoms[Cnum], &s_atoms[(*number_of_s_atoms)]);
		(*number_of_s_atoms)++;
	}
	else
		control = 1;
	return control;

}

void find_select_atoms(pdb_atom_pt atoms, int number_of_atoms, pdb_atom_pt b_atoms, int number_of_b_atoms, pdb_atom_pt s_atoms, int *number_of_s_atoms) 
{
	int i, j, find;
	int Anum, Cnum;
	int residue_number_A, residue_number_B;
	float x_A, y_A, z_A;
	float x_B, y_B, z_B;
	float x, y, z;
	double distance;
	char pdbline_A[MAX_LINE_LENGTH];
	char pdbline_B[MAX_LINE_LENGTH];
	char atom_A_line[MAX_LINE_LENGTH];
	char atom_C_line[MAX_LINE_LENGTH];
	char atom_name_A[5], residue_name_A[5], residue_insertion_A;
	char match_A[5], chain_id_A;
	char atom_name_B[5], residue_name_B[5], residue_insertion_B;
	char match_C[5], chain_id_B;

	find=0;
	(*number_of_s_atoms)=0;
	/** start looking for oxygens and nitrogens in the binding site */
	for(i=0; i<number_of_b_atoms; i++)
	{
		if(b_atoms[i].type==UNKNOWN) /** not hetero **/
		{
			strcpy(atom_name_B, b_atoms[i].atom_name);
			strcpy(residue_name_B, b_atoms[i].residue_name);
			chain_id_B=b_atoms[i].chainID;
			residue_number_B = atoi(b_atoms[i].residue_number);
			residue_insertion_B=b_atoms[i].icode;
			x_B=b_atoms[i].pos[0];
			y_B=b_atoms[i].pos[1];
			z_B=b_atoms[i].pos[2];

			create_pdb_line(&b_atoms[i], pdbline_B, 0);

			if (!strncmp(atom_name_B, " O  ", 4))
			{  /* found a carbonyl oxygen.  */

				/** insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);	
				      (*number_of_s_atoms)++; **/

				strncpy(atom_A_line,"NOT FOUND ", 10);
				strncpy(atom_C_line,"NOT FOUND ", 10);


				/* let's find the neighbors */

				for(j=0; j<number_of_atoms; j++)
				{
					if(atoms[j].type==UNKNOWN) /** not hetero **/
					{
						strcpy(atom_name_A, atoms[j].atom_name);
						strcpy(residue_name_A, atoms[j].residue_name);
						chain_id_A = atoms[j].chainID;
						residue_number_A = atoi(atoms[j].residue_number);
						residue_insertion_A=atoms[j].icode;
						x_A=atoms[j].pos[0];
						y_A=atoms[j].pos[1];
						z_A=atoms[j].pos[2];

						create_pdb_line(&atoms[j], pdbline_A, 0);

						/*  calculate distance between atoms  A and B */

						distance = sqrt(((x_A - x_B)*(x_A - x_B)) + 
						    ((y_A - y_B)*(y_A - y_B)) + 
						    ((z_A - z_B)*(z_A - z_B)));

						if ((chain_id_A == chain_id_B)&&
						    (residue_number_A == residue_number_B)&&
						    (!strncmp(residue_name_A, residue_name_B, 3))&&
						    (!strncmp(atom_name_A, " C  ", 4))&&
						    (distance < THRESHOLD1)){

							strncpy(atom_A_line, pdbline_A, 54 );
							atom_A_line[54]='\0';
							x=x_A; y=y_A; z=z_A;
							Anum=j;
						}

						if ((chain_id_A == chain_id_B)&&
						    (distance < THRESHOLD2)&&
						    (!strncmp(atom_name_A, " N  ", 4))){

								distance = sqrt(((x_A - x)*(x_A - x)) + 
						   	   	((y_A - y)*(y_A - y)) + 
						       	((z_A - z)*(z_A - z)));
								if(distance<THRESHOLD1) {
							/* found backbone nitrogen of following residue */
								strncpy(atom_C_line, pdbline_A, 54);
								atom_C_line[54]='\0';
								Cnum=j;
							}
						}
					}
				}

				if((strncmp(atom_A_line, "NOT FOUND ", 10))&&
				    (strncmp(atom_C_line, "NOT FOUND ", 10))){
					insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
					insert_list(&atoms[Anum], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
					insert_list(&atoms[Cnum], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
				}
				else 
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n",pdbline_B);

			}

                       
			if (!strncmp(atom_name_B, " N  ", 4) 
                            && strcmp(residue_name_B, "PRO") ) {
				/**insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);	
				      (*number_of_s_atoms)++; **/

				/* Neighbors initially not known. Need this to check    */
				/* later if they were found or not.                     */

				strncpy(atom_A_line,"NOT FOUND ", 10);
				strncpy(atom_C_line,"NOT FOUND ", 10);

				/* found a main chain nitrogen   */
				/* now let's find the neighbors  */

				for(j=0; j<number_of_atoms; j++)
				{
					if(atoms[j].type==UNKNOWN) /** not hetero **/
					{
						strcpy(atom_name_A, atoms[j].atom_name);
						strcpy(residue_name_A, atoms[j].residue_name);
						chain_id_A = atoms[j].chainID;
						residue_number_A = atoi(atoms[j].residue_number);
						residue_insertion_A=atoms[j].icode;
						x_A=atoms[j].pos[0];
						y_A=atoms[j].pos[1];
						z_A=atoms[j].pos[2];

						create_pdb_line(&atoms[j], pdbline_A, 0);

						distance = sqrt(((x_A - x_B)*(x_A - x_B)) + 
						    ((y_A - y_B)*(y_A - y_B)) + 
						    ((z_A - z_B)*(z_A - z_B)));


						if ((chain_id_A == chain_id_B)&& 
						    (distance < THRESHOLD1)&&
						    (!strncmp(atom_name_A, " C  ", 4))){
							/* Found alpha carbon of previous residue */
							strncpy(atom_A_line, pdbline_A, 54);
							atom_A_line[54]='\0';
							x=x_A; y=y_A; z=z_A;
							Anum=j;
						}
						if ((chain_id_A == chain_id_B)&&
						    (distance <= THRESHOLD2)&&
						    (!strncmp(atom_name_A, " O  ", 4))){
							distance = sqrt(((x_A - x)*(x_A - x)) + 
						   	   	((y_A - y)*(y_A - y)) + 
						       	((z_A - z)*(z_A - z)));
								if(distance<THRESHOLD2) {
							/* Found carbonyl oxygen of previous residue */
								strncpy(atom_C_line, pdbline_A, 54);
								atom_C_line[54]='\0';
								Cnum=j;
							}
						}
					}
				}
				if((strncmp(atom_A_line, "NOT FOUND ", 10))&&
				    (strncmp(atom_C_line, "NOT FOUND ", 10))){
					insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
					insert_list(&atoms[Anum], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
					insert_list(&atoms[Cnum], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
				}
				else fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);

			}

			if ((!strncmp(atom_name_B, " OD1", 4))&&
			    (!strncmp(residue_name_B, "ASP", 3))){

				strcpy(match_A," CG ");
				strcpy(match_C," OD2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OD2", 4))&&
			    (!strncmp(residue_name_B, "ASP", 3))){

				strcpy(match_A," CG ");
				strcpy(match_C," OD1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OD1", 4))&&
			    (!strncmp(residue_name_B, "ASN", 3))){

				strcpy(match_A," CG ");
				strcpy(match_C," ND2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms); 
				if (find == 1) 
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " ND2", 4))&&
			    (!strncmp(residue_name_B, "ASN", 3))){

				strcpy(match_A," CG ");
				strcpy(match_C," OD1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OE1", 4))&&
			    (!strncmp(residue_name_B, "GLU", 3))){

				strcpy(match_A," CD ");
				strcpy(match_C," OE2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OE2", 4))&&
			    (!strncmp(residue_name_B, "GLU", 3))){

				strcpy(match_A," CD ");
				strcpy(match_C," OE1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OE1", 4))&&
			    (!strncmp(residue_name_B, "GLN", 3))){

				strcpy(match_A," CD ");
				strcpy(match_C," NE2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NE2", 4))&&
			    (!strncmp(residue_name_B, "GLN", 3))){

				strcpy(match_A," CD ");
				strcpy(match_C," OE1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NH1", 4))&&
			    (!strncmp(residue_name_B, "ARG", 3))){

				strcpy(match_A," CZ ");
				strcpy(match_C," NH2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NH2", 4))&&
			    (!strncmp(residue_name_B, "ARG", 3))){

				strcpy(match_A," CZ ");
				strcpy(match_C," NH1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NE", 3))&&
			    (!strncmp(residue_name_B, "ARG", 3))){

				strcpy(match_A," CZ ");
				strcpy(match_C," CD ");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " ND1", 4))&&
			    (!strncmp(residue_name_B, "HIS", 3))){

				strcpy(match_A," CE1");
				strcpy(match_C," NE2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NE2", 4))&&
			    (!strncmp(residue_name_B, "HIS", 3))){

				strcpy(match_A," CE1");
				strcpy(match_C," ND1");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NE1", 4))&&
			    (!strncmp(residue_name_B, "TRP", 3))){

				strcpy(match_A," CD1");
				strcpy(match_C," CE2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OH", 3))&&
			    (!strncmp(residue_name_B, "TYR", 3))){

				strcpy(match_A," CZ ");
				strcpy(match_C," CE2");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " NZ", 3))&&
			    (!strncmp(residue_name_B, "LYS", 3))){

				strcpy(match_A," CE ");
				strcpy(match_C," CD ");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OG", 3))&&
			    (!strncmp(residue_name_B, "SER", 3))){

				strcpy(match_A," CB ");
				strcpy(match_C," CA ");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}

			if ((!strncmp(atom_name_B, " OG", 3))&&
			    (!strncmp(residue_name_B, "THR", 3))){

				strcpy(match_A," CB ");
				strcpy(match_C," CA ");
				insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
				(*number_of_s_atoms)++;
				find = get_neighbors(atom_name_B, residue_name_B, chain_id_B, residue_number_B, match_A, match_C, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
				if (find == 1)
				{
					fprintf(stderr, "Warning: Neighbors of %.30s were not found\n", pdbline_B);
					(*number_of_s_atoms)--;
				}
			}
		}

		/* Finding ligand atoms that can form H-bonds *
		       * and will be part of the binding site.      *
		       * Waters are not included yet.               */

		if(b_atoms[i].type==HETERO || b_atoms[i].act == METAL_1 || b_atoms[i].act == METAL_2 ) /** HETATM **/
		{
			strcpy(atom_name_B, b_atoms[i].atom_name);
			strcpy(residue_name_B, b_atoms[i].residue_name);
			x_B=b_atoms[i].pos[0];
			y_B=b_atoms[i].pos[1];
			z_B=b_atoms[i].pos[2];

			create_pdb_line(&b_atoms[i], pdbline_B, 1);

			if (strncmp(residue_name_B, " HOH",4)) {    /*non-HOH ligand found*/
				if ((!strncmp(atom_name_B, " ONN",4))||
				    (!strncmp(atom_name_B, " OAA",4))||
				    (!strncmp(atom_name_B, " NAA",4))||
				    (!strncmp(atom_name_B, " NDD",4))) {

					insert_list(&b_atoms[i], &s_atoms[(*number_of_s_atoms)]);
					(*number_of_s_atoms)++;
					find = get_hetero_neighbors(atom_name_B, x_B, y_B, z_B, atoms, number_of_atoms, s_atoms, number_of_s_atoms);
					if (find == 1)
					{
						fprintf(stderr, "Warning: Neighbors of %.30s were not found\n",pdbline_B);
						(*number_of_s_atoms)--;
					}
				}
			}
		}
	}
}
