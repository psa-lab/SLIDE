/** This file checks for bumps between the template points and
    atoms in pdb file. The subset of atoms from the pdb file that 
    are close to the template points is defined by the radfile.
                      --RSK                                **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

double sqre(double x)
{
	return x*x;
}

double find_dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt( sqre(x1-x2) + sqre(y1-y2) + sqre(z1-z2));
}

/** copy atomB to atomA **/
void pdbeq(pdb_atom_pt atomA, pdb_atom_pt atomB)
{
	strcpy(atomA->atom_name, atomB->atom_name);
	strcpy(atomA->atom_number, atomB->atom_number);
	strcpy(atomA->residue_name, atomB->residue_name);
	strcpy(atomA->residue_number, atomB->residue_number);
	atomA->icode = atomB->icode;
	atomA->chainID = atomB->chainID;

	atomA->pos[0] = atomB->pos[0];
	atomA->pos[1] = atomB->pos[1];
	atomA->pos[2] = atomB->pos[2];

	atomA->type = atomB->type;
	atomA->residue_type = atomB->residue_type;
	atomA->hphil = atomB->hphil;
	atomA->act = atomB->act;
	atomA->class = atomB->class;
	atomA->rad = atomB->rad;
}

/** for a surface point x, y, z -  find binding site atoms and rad atoms **/
void get_atoms( double x, double y, double z, pdb_atom_pt atoms, int number_of_atoms, pdb_atom_pt b_atoms, int *number_of_b_atoms, pdb_atom_pt r_atoms, int *number_of_r_atoms, short bselected[MAX_NUMBER_OF_ATOMS], short rselected[MAX_NUMBER_OF_ATOMS])
{

	int j;
	double dist;

	for(j=0;j<number_of_atoms;j++)
	{
		dist=find_dist(x, y, z, atoms[j].pos[0], atoms[j].pos[1], atoms[j].pos[2]);

		if(dist<MAXRADDIST && rselected[j]==0)
		{
			pdbeq(&r_atoms[(*number_of_r_atoms)], &atoms[j]);
			(*number_of_r_atoms)++;
			rselected[j]=1;
		}
		if(dist<MAXBINDDIST && bselected[j]==0)
		{
			pdbeq(&b_atoms[(*number_of_b_atoms)], &atoms[j]);
			(*number_of_b_atoms)++;
			bselected[j]=1;
		}
	}
}

/** check if x,y,z bumps into any of the NN points **/
int check_bump(double x, double y, double z,  pdb_atom_pt r_atoms, int nrad,
double maxdist)
{
	int j;
	double dist;
	/*double mindist = 10000;*/

	for(j=0;j<nrad;j++)
	{
		dist=find_dist(x, y, z, r_atoms[j].pos[0], r_atoms[j].pos[1], r_atoms[j].pos[2]);
		/*if(dist<mindist) mindist=dist;*/

		if(dist<maxdist)
		{
			/**	printf( "(%.3f) %.3f %.3f %.3f - [%d] %.3f %.3f %.3f\n",  dist, x, y, z, j, NN[0][j], NN[1][j], NN[2][j]); **/
			break;
		}
	}

	if(j==nrad) { /*printf("(%.3f) ", mindist);*/
		return -1;
	}
	else return j;
}
