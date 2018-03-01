/* This program creates template points around metals that are         */
/* in a binding site and are supposed to interact with the ligand.     */
/* Metal ions types included: Ca, Co, Cu, Fe, K, Mg, Mn, Na, Zn.       */
/* Input: 1) pdb file with coordinates of the metals                   */
/*        2) pdb file of the target binding site without the metal     */
/*        3) file with the coordinates of 32 points evenly             */
/*           distributed on the surface of a sphere: points_on_sphere  */
/*	     A sample file is provided in: "slide/src/template/"       */
/* Output:   coordinates of acceptor points around the metal ion       */
/*           at least 2.5 A from any protein atom (except the metal)   */   

/***  MIZ January 31, 2005  ***/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BUFSIZE 256

#define RAD_CA 2.4 
#define RAD_CO 1.9
#define RAD_CU 2.1 
#define RAD_FE 2.2
#define RAD_K  2.4 
#define RAD_MG 2.1
#define RAD_MN 2.2
#define RAD_NA 2.4
#define RAD_ZN 2.1  /* 2.2 */


int main(int argc, char *argv[]) {

  FILE *metal_file_pdb;
  FILE *bind_site_pdb;
  FILE *template_file;

  float xm[20], ym[20], zm[20]; 
  float x, y, z;
  float rad;
  float xp[10000], yp[10000], zp[10000]; /* coordinates of protein atom    */
  float xtemp[20][32], ytemp[20][32], ztemp[20][32];       /* coordinates of metal template points */
  float xtemp2, ytemp2, ztemp2; /* coordinates of old template file */
  float distance;                     /* square pf protein_atom - metal template point distance */
  float distance2;                     /* protein_atom - metal template point distance */

  char strbuf[BUFSIZE+1];
  char metal_file[80];
  char bind_site[80];
  char template[80];
  char points[80];
  char atom_name[10];
  char target_atom_name[10];
  int i, j, k;
  int no_metal_points;
  int no_temp_points;
  int no_target_atoms;
  int unbump;

  /* Check the command line arguments */
  if (argc != 4) {
    fprintf(stderr, 
	    "Usage: %s <target_metal.pdb> <binding_site.pdb> <original template>\n", argv[0]);
    return -1;
  }
  /* Get the command line arguments */
  strncpy(metal_file, argv[1], 79);
  strncpy(bind_site, argv[2], 79); 
  strncpy(template, argv[3], 79);


  metal_file_pdb = fopen(metal_file, "r");
  if (metal_file_pdb == 0) {
    fprintf(stderr, "Unable to open input file: %s\n", metal_file);
    return -1;
  }

  bind_site_pdb = fopen(bind_site, "r");
  if (bind_site_pdb == 0) {
    fprintf(stderr, "Unable to open input file: %s\n", bind_site);
    return -1;
  }

  template_file = fopen(template, "r");
  if (template_file == 0) {
    fprintf(stderr, "Unable to open input file: %s\n", template);
    return -1;
  }

  i = 0;
  no_metal_points = 0;
  while(fgets(strbuf, sizeof(strbuf), metal_file_pdb)!=NULL){
    sscanf(strbuf, "%*12c %s", atom_name);
    sscanf(strbuf, "%*31c %f %f %f", &xm[no_metal_points], &ym[no_metal_points], &zm[no_metal_points]);
           printf("%s %f %f %f\n", atom_name,  xm[no_metal_points], ym[no_metal_points], zm[no_metal_points]); 
    i = 0;
    if((!strcmp(atom_name, "CA"))||(!strcmp(atom_name, "Ca"))) {
      rad = RAD_CA;
    }
    if((!strcmp(atom_name, "CO"))||(!strcmp(atom_name, "Co"))) {
      rad = RAD_CO;
    }
    if((!strcmp(atom_name, "CU"))||(!strcmp(atom_name, "Cu"))) {
      rad = RAD_CU;
    }
    if((!strcmp(atom_name, "FE"))||(!strcmp(atom_name, "Fe"))) {
      rad = RAD_FE;
    }
    if((!strcmp(atom_name, "K"))) {
      rad = RAD_K;
    }
    if((!strcmp(atom_name, "MG"))||(!strcmp(atom_name, "Mg"))) {
      rad = RAD_MG;
    }
    if((!strcmp(atom_name, "MN"))||(!strcmp(atom_name, "Mn"))) {
      rad = RAD_MN;
    }
    if((!strcmp(atom_name, "NA"))||(!strcmp(atom_name, "Na"))) {
      rad = RAD_NA;
    }
    if((!strcmp(atom_name, "ZN"))||(!strcmp(atom_name, "Zn"))) {
      rad = RAD_ZN;
    }
    /*     printf("%f\n", rad);*/

    /*    while(fgets(strbuf, sizeof(strbuf), points_on_sphere)!=NULL){
      sscanf(strbuf, "%*5c %f %f %f", &x, &y, &z);
 
      xtemp[i] = rad*x + xm;
      ytemp[i] = rad*y + ym; 
      ztemp[i] = rad*z + zm;*/
      /*            printf("%d  %8.3f %8.3f %8.3f\n", i, xtemp[i], ytemp[i], ztemp[i]); */
      /*      printf("xtemp[i] = rad*%f + xm;\n", x); 
      printf("ytemp[i] = rad*%f + ym;\n", y); 
      printf("ztemp[i] = rad*%f + zm;\n", z);
      printf("++i;\n");
      ++i;*/
      xtemp[no_metal_points][i] = rad*0.299121 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.230511 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.925954 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.810603 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.052722 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.583217 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.399627 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.869602 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.289980 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.870867 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.254322 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.420607 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.025394 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.747346 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.663949 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.399627 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.869602 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.289980 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.111304 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.990768 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.077386 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.152206 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.782372 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.603925 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.456100 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.447749 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.769085 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.973609 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.227942 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.011271 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.639578 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.745007 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.189482 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.872956 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.449337 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.189848 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.152206 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.782372 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.603925 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.973609 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.227942 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.011271 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.522281 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.784941 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.333300 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.299121 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.230511 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.925954 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.618959 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.648112 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.443666 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.639578 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.745007 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.189482 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.522281 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.784941 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.333300 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.122903 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.250537 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.960273 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.406204 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.315076 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.857744 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.111304 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.990768 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.077386 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.685140 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.310920 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.658719 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.025394 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.747346 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.663949 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.872956 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.449337 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.189848 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.685140 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.310920 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.658719 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.618959 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.648112 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.443666 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.810603 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.052722 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*0.583217 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.122903 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.250537 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.960273 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*-0.456100 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*-0.447749 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.769085 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.406204 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.315076 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.857744 + zm[no_metal_points];
      ++i;
      xtemp[no_metal_points][i] = rad*0.870867 + xm[no_metal_points];
      ytemp[no_metal_points][i] = rad*0.254322 + ym[no_metal_points];
      ztemp[no_metal_points][i] = rad*-0.420607 + zm[no_metal_points];
      ++i;
      ++no_metal_points;
      
      /* } */
  }
  /*  printf ("# metal points: %d\n", no_metal_points);*/
  no_temp_points = i;
  /*  printf("# temp points: %d\n", no_temp_points);*/
  /*  printf("DEBUG INFO: CHK PT 1: # of template points %d\n", no_temp_points);*/
  j = 0;
  while (fgets(strbuf, sizeof(strbuf), bind_site_pdb)!=NULL){
    if (!strncmp(strbuf, "ATOM", 4)) {        
      
      sscanf(strbuf, "%*12c %s", target_atom_name);
      /*      printf ("j/name %d %s\n", j, target_atom_name);*/
      if(strncmp(target_atom_name, "H",  1)){ 
	sscanf(strbuf, "%*30c %f %f %f", &xp[j], &yp[j], &zp[j]);
	++j;
      }
    }
  }
  
  no_target_atoms = j;
  /*  printf ("# metal points %d\n", no_metal_points);
  printf ("# temp points %d\n", no_temp_points);
  printf ("# target_atoms %d\n", no_target_atoms);*/
  for (k = 0; k < no_metal_points; ++k){
    for (i = 0; i < no_temp_points; ++i){ 
      unbump = 0; 
      for (j = 0; j < no_target_atoms; ++j){ 	  
	/*      distance = ((xtemp[i]-xp[j])*(xtemp[i]-xp[j])+
		(ytemp[i]-yp[j])*(ytemp[i]-yp[j])+
		(ztemp[i]-zp[j])*(ztemp[i]-zp[j]));*/
	distance2 = sqrt((xtemp[k][i]-xp[j])*(xtemp[k][i]-xp[j])+
			 (ytemp[k][i]-yp[j])*(ytemp[k][i]-yp[j])+
			 (ztemp[k][i]-zp[j])*(ztemp[k][i]-zp[j]));
	/*      printf ("distance %f %f\n", distance, distance2);*/
	if (distance2 < 2.5 ){
	  /*		printf ("%d  %8.3f %8.3f %8.3f\n",i, xp[j], yp[j], zp[j]); 
			printf ("%d  %8.3f %8.3f %8.3f\n",i, xtemp[i], ytemp[i], ztemp[i]); 
			printf ("%d  %8.3f\n", i, distance2); 	*/
	  unbump = 1;   
	  break;
	}
      }
      if (unbump == 0){
	printf ("A*  %7.3f %7.3f %7.3f\n", xtemp[k][i], ytemp[k][i], ztemp[k][i]); 
      }     
    }
  }

  /* start of new template pruning */
  while (fgets(strbuf, sizeof(strbuf), template_file)!=NULL){  
    sscanf(strbuf, "%*4c %f %f %f", &xtemp2, &ytemp2, &ztemp2);
    unbump = 0;
    j = 0;
    distance = 10.0;
    for (j = 0; j < no_metal_points; ++j){ 
      distance = (xtemp2-xm[j])*(xtemp2-xm[j])+
	(ytemp2-ym[j])*(ytemp2-ym[j])+
	(ztemp2-zm[j])*(ztemp2-zm[j]);
      
      if (distance < 4.0 ){ 
	unbump = 1;
      }
    }
    if (unbump == 0){
      printf ("%s", strbuf); 
    }     
  }
  /* end of new template pruning */
  
  
  /* Clean up and exit */
  fclose(metal_file_pdb);
  fclose(bind_site_pdb);
  /* fclose(points_on_sphere);*/
  fclose(template_file);
  return (0);
}
