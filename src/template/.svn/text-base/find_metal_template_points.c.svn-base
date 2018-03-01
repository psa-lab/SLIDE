#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "types.h"
#include "basics.h"

void find_metal_template_points ( pdb_atom_pt bs_atoms,
				 int         number_of_bs_atoms,
				 point_pt metal_template_points,
				 int         *number_of_metal_template_points )
{
  FILE *selfile;
  int i, j, k, number_of_metals;
  float distance, xtemp[20][32], ytemp[20][32], ztemp[20][32];
  int unbump, number_of_metal_template_points_per_metal;
  int mp;

  unbump = 0;

  number_of_metals = 0;
  mp = 0;
  i = j = k = 0;
  distance = 0.0;
  
  /*  if((selfile=fopen( SEL_ATOMS, "w"))==NULL)
    { 
      fprintf(stderr,"ERROR: Cannot find %s\n", SEL_ATOMS); 
      exit(1); 
    }
  */
  
  for (i = 0;i<number_of_bs_atoms;i++)
    {
      if(/*bs_atoms[i].type==HETERO ||*/ bs_atoms[i].act == METAL_1 || bs_atoms[i].act == METAL_2 )
	{
	  j = 0;
	  /*	  printf("HETERO%5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f %f\n", bs_atoms[i].atom_number, bs_atoms[i].atom_name, bs_atoms[i].residue_name, bs_atoms[i].chainID, bs_atoms[i].residue_number, bs_atoms[i].icode, bs_atoms[i].pos[0], bs_atoms[i].pos[1], bs_atoms[i].pos[2], bs_atoms[i].rad);*/
	  /*	  fprintf(selfile, "HETERO%5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f %f\n", bs_atoms[i].atom_number, bs_atoms[i].atom_name, bs_atoms[i].residue_name, bs_atoms[i].chainID, bs_atoms[i].residue_number, bs_atoms[i].icode, bs_atoms[i].pos[0], bs_atoms[i].pos[1], bs_atoms[i].pos[2], bs_atoms[i].rad);*/

	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.299121 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.230511 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.925954 + bs_atoms[i].pos[Z];
	  /*	  fprintf(selfile, "A*  %8.3f%8.3f%8.3f\n", xtemp[number_of_metals][j], ytemp[number_of_metals][j], ztemp[number_of_metals][j]);*/
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.810603 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.052722 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.583217 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.399627 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.869602 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.289980 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.870867 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.254322 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.420607 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.025394 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.747346 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.663949 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.399627 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.869602 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.289980 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.111304 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.990768 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.077386 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.152206 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.782372 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.603925 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.456100 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.447749 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.769085 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.973609 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.227942 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.011271 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.639578 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.745007 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.189482 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.872956 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.449337 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.189848 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.152206 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.782372 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.603925 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.973609 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.227942 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.011271 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.522281 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.784941 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.333300 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.299121 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.230511 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.925954 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.618959 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.648112 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.443666 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.639578 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.745007 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.189482 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.522281 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.784941 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.333300 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.122903 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.250537 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.960273 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.406204 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.315076 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.857744 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.111304 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.990768 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.077386 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.685140 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.310920 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.658719 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.025394 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.747346 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.663949 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.872956 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.449337 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.189848 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.685140 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.310920 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.658719 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.618959 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.648112 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.443666 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.810603 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.052722 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*0.583217 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.122903 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.250537 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.960273 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*-0.456100 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*-0.447749 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.769085 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.406204 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.315076 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.857744 + bs_atoms[i].pos[Z];
	  ++j;
	  xtemp[number_of_metals][j] = bs_atoms[i].rad*0.870867 + bs_atoms[i].pos[X];
	  ytemp[number_of_metals][j] = bs_atoms[i].rad*0.254322 + bs_atoms[i].pos[Y];
	  ztemp[number_of_metals][j] = bs_atoms[i].rad*-0.420607 + bs_atoms[i].pos[Z];
	  ++j;
	  ++number_of_metals;
	}
      number_of_metal_template_points_per_metal = j;
    }
  for (k = 0; k < number_of_metals; k++)
    {
      for (j = 0; j < number_of_metal_template_points_per_metal; j++)
	{
	  unbump = 0;
	  for (i = 0; i < number_of_bs_atoms; i++)
	    {
	      if(bs_atoms[i].act != METAL_1 && bs_atoms[i].act != METAL_2 )
		{
		  distance = sqrt((xtemp[k][j]-bs_atoms[i].pos[X])*(xtemp[k][j]-bs_atoms[i].pos[X])+(ytemp[k][j]-bs_atoms[i].pos[Y])*(ytemp[k][j]-bs_atoms[i].pos[Y])+(ztemp[k][j]-bs_atoms[i].pos[Z])*(ztemp[k][j]-bs_atoms[i].pos[Z]));
		  if (distance < 2.5)
		    {
		      //printf("bump\n");
		      unbump = 1;
		      break;
		    }
		}
	    }
	  if (unbump == 0)
	    {
	      /*	      fprintf(selfile, "A*  %8.3f%8.3f%8.3f\n", xtemp[k][j], ytemp[k][j], ztemp[k][j]);	      */
	      metal_template_points[mp].pos[X] = xtemp[k][j];
	      metal_template_points[mp].pos[Y] = ytemp[k][j];
	      metal_template_points[mp].pos[Z] = ztemp[k][j];
	      mp++;
	    }
	}
    }
  printf("%d metal acceptor points\n", mp);
  (*number_of_metal_template_points) = mp;
  /*  fclose(selfile);*/
}
