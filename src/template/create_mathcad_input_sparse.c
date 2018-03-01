/* This program reads in the coordinates of binding-site atoms  */
/* from a pdb file and creates a mathcad input file. The pdb    */
/* file should not contain empty lines!                         */
/* SPARSE version */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "defs.h"
#include "types.h"

void create_mathcad_input_sparse(pdb_atom_pt s_atoms, int number_of_s_atoms)
{
	FILE *outfile;
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	int i, j;
	float l;             /* H-bond length */
	char atom_name[6];
	char res_name[6];
	char strbuf[MAX_LINE_LENGTH];
	char output_line[MAX_LINE_LENGTH];
	char output2[MAX_LINE_LENGTH];
	char output3[MAX_LINE_LENGTH];

	if((outfile = fopen(MCDINFILE, "w"))==NULL) {
		fprintf(stderr, "File: Cannot open %s\n", MCDINFILE);
		exit(1);
	}

	j = 0;
	for(i=0; i<number_of_s_atoms; i++)
	{
		if(s_atoms[i].type==UNKNOWN) {
			/* AHA!  We have an ATOM record */

			/* we read in the atom name and x,y,z coordinates  by   */
			/* character position: atom name starts at position 14  */
			/* followed by residue name. After that 10 characters   */
			/* are ignored. The following three strings will be     */
			/* the coordinates.                                     */

			strcpy(atom_name, s_atoms[i].atom_name);
			strcpy(res_name, s_atoms[i].residue_name);
			x2 = s_atoms[i].pos[0];
			y2 = s_atoms[i].pos[1];
			z2 = s_atoms[i].pos[2];

			/* It is not a mistake: first line contains atom B,      */
			/* the second atom A  and the third atom C.              */

			i++;
			x1 = s_atoms[i].pos[0];
			y1 = s_atoms[i].pos[1];
			z1 = s_atoms[i].pos[2];

			i++;
			x3 = s_atoms[i].pos[0];
			y3 = s_atoms[i].pos[1];
			z3 = s_atoms[i].pos[2];

			/* we make a string of the nine coordinates necessary  */
			/* to calculate the position of the template point     */
			/* corresponding to  atom B and then print it out.     */

			sprintf(output_line,"%d %7.3lf %7.3lf %7.3lf %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f", j, x1, y1, z1,
			    x2, y2, z2, x3, y3, z3);

			if((strcmp(atom_name," O  ")!=0)&&(strcmp(atom_name," N  ")!=0)){

				if((strcmp(res_name, "ASP") == 0)||
				    (strcmp(res_name, "GLU") == 0)||
				    (strcmp(res_name, "LYS") == 0)||
				    (strcmp(res_name, "THR") == 0)||
				    (strcmp(res_name, "TYR") == 0)) {
					l = 2.9;
				}
				else
					l = 3.0;
			}
			else l = 3.0;

			/* first column: 1 if TEMPLATE POINT (!) can be an acceptor,    */
			/* 2 if donor, 3 if "doneptor".                                 */
			/* number in second column: angle alpha (A-B-template point)    */
			/* third: anlge beta between plane ABC and B-template point     */
			/* forth: distance l between B and template point (bond length) */

			if(strcmp(atom_name, " O  ") == 0) {    /* main chain O */
				fprintf(outfile, "2 140   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 120   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 160   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 200   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 180   0 %3.1f %s\n", l, output_line);

			}

			if(strcmp(atom_name, " N  ") == 0) {    /* main chain N */
				fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 260   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 230  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 250  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 230 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 250 -30 %3.1f %s\n", l, output_line);
			}

			if((strncmp(atom_name, " OD",3) == 0)||
			    (strncmp(atom_name, " OE",3) == 0)) {
				/* Glu, Gln, Asp and Asn O  */

				fprintf(outfile, "2 110   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 110  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 110 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 130   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 200   0 %3.1f %s\n", l, output_line);
			}

			if(strncmp(atom_name, " NH",3) == 0){        /* Arg */

				fprintf(outfile, "1 120   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240  30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240 -30 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 145   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
			}

			if(strncmp(atom_name, " ND",3) == 0){

				/* ND could be from His or Asn */

				if(strcmp(res_name, "HIS") == 0){
					fprintf(outfile, "3 220   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 240   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 260   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 230  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 250  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 230 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 250 -30 %3.1f %s\n", l, output_line);
				}
				if(strcmp(res_name, "ASN") == 0) {
					fprintf(outfile, "1 120   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 120  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 120 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 100   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 145   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
				}
			}

			if(strncmp(atom_name, " NE",3) == 0){

				/* NE could be from His,Trp,Arg or Gln */

				if((strcmp(res_name, "TRP") == 0)||
				    (strcmp(res_name, "ARG") == 0)){

					fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 230  15 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 250  15 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 230 -15 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 250 -15 %3.1f %s\n", l, output_line);
				}

				if(strcmp(res_name, "HIS") == 0) {
					fprintf(outfile, "3 220   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 240   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 260   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 230  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 250  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 230 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "3 250 -30 %3.1f %s\n", l, output_line);

				}
				if(strcmp(res_name, "GLN") == 0) {
					fprintf(outfile, "1 120   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 120  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240  30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 120 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 240 -30 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 100   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 145   0 %3.1f %s\n", l, output_line);
					fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
				}
			}

			if(strcmp(atom_name, " OH ") == 0) {       /*  Tyr  */

				fprintf(outfile, "3 120   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 110  15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 130  15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 230  15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 250  15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 110 -15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 130 -15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 230 -15 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 250 -15 %3.1f %s\n", l, output_line);
			}

			if(strncmp(atom_name, " OG",3) == 0) {   /* Ser, Thr */

				fprintf(outfile, "3 105 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105  70 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 290 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 310 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125  70 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 290 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 310 %3.1f %s\n", l, output_line);
			}

			if(strcmp(atom_name, " NZ ") == 0) {       /* Lys */

				fprintf(outfile, "1 100 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100  70 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100 290 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100 310 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120  70 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 290 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 310 %3.1f %s\n", l, output_line);
			}

		}


		/** this part will create the template points for heteroatoms **/

		else if(s_atoms[i].type==HETERO || s_atoms[i].act == METAL_1 || s_atoms[i].act == METAL_2) {
			/* AHA!  We also have heteroatoms as part of the binding site! */

			strcpy(atom_name, s_atoms[i].atom_name);
			strcpy(res_name, s_atoms[i].residue_name);
			x2 = s_atoms[i].pos[0];
			y2 = s_atoms[i].pos[1];
			z2 = s_atoms[i].pos[2];

			/* It is not a mistake: first line contains atom B,      */
			/* the second atom A  and the third atom C.              */

			i++;
			x1 = s_atoms[i].pos[0];
			y1 = s_atoms[i].pos[1];
			z1 = s_atoms[i].pos[2];

			i++;
			x3 = s_atoms[i].pos[0];
			y3 = s_atoms[i].pos[1];
			z3 = s_atoms[i].pos[2];


			/* we make a string of the nine coordinates necessary  */
			/* to calculate the position of the template point     */
			/* corresponding to  atom B and then print it out.     */

			sprintf(output_line,"%d %7.3lf %7.3lf %7.3lf",i, x1, y1, z1);
			sprintf(output2," %7.3lf %7.3lf %7.3lf", x2, y2, z2);
			sprintf(output3," %7.3lf %7.3lf %7.3lf", x3, y3, z3);

			strcat(output_line, output2);
			strcat(output_line, output3);
			l = 3.0;

			if(strcmp(atom_name, " OAA") == 0){
				/* acceptor oxygen, similar to main chain carbonyl O */
				fprintf(outfile, "2 140   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 140 -50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220 -50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 160   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 200   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 180   0 %3.1f %s\n", l, output_line);
			}
			if(strcmp(atom_name, " ONN") == 0){
				/* doneptor oxygen, as in Ser of Thr */
				fprintf(outfile, "3 105  60 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 300 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105  40 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105  80 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 280 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 105 320 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125  60 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 180 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 300 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125  40 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125  80 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 280 %3.1f %s\n", l, output_line);
				fprintf(outfile, "3 125 320 %3.1f %s\n", l, output_line);
			}
			if(strcmp(atom_name, " N1D") == 0){
				/* nitrogen that can donate one H, similar to Trp NE */
				fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 260   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 220  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 260  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 220 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 260 -25 %3.1f %s\n", l, output_line);
			}
			if(strcmp(atom_name, " N2D") == 0){
				/* nitrogen that can donate two H-s, similar to Arg NH */
				fprintf(outfile, "1 120   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240  50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 120 -50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 240 -50 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 100   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 145   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "1 220   0 %3.1f %s\n", l, output_line);
			}
			if(strcmp(atom_name, " NAA") == 0){
				/* acceptor nitrogen, similar geometry as for His NE/ND */
				fprintf(outfile, "2 220   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 240   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 260   0 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 240  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 260  25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 220 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 240 -25 %3.1f %s\n", l, output_line);
				fprintf(outfile, "2 260 -25 %3.1f %s\n", l, output_line);
			}
		}
		j++;
	}
	fclose(outfile);
}
