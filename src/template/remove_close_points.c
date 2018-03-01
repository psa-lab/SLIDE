/* This program removes template points closer than 2.0 A to           */
/* binding-site metals metals                                          */

/* Input: -  template file                                             */
/*        -  pdb file of the metal                                     */
/* Output: Template file with points further than 2.0 A from the metal */

/***  MIZ January 31, 2005  ***/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#define BUFSIZE 256


int main(int argc, char *argv[]) {

  FILE *template_file;
  FILE *metal_file;

  float xm[10], ym[10], zm[10];   /* coordinates of metal atoms    */
  float xtemp, ytemp, ztemp;      /* coordinates of template points */
  float distance;

  char strbuf[BUFSIZE+1];
  char metal[80];
  char template[80];
  char atom_name[10];
  int i, j;
  int no_temp_points;
  int no_metals;
  int unbump;

  /* Check the command line arguments */
  if (argc != 3) {
    fprintf(stderr, 
	    "Usage: %s <metal_file.pdb> <template_file>\n", argv[0]);
    return -1;
  }

  /* Get the command line arguments */
  strncpy(metal, argv[1], 79);
  strncpy(template, argv[2], 79); 


  metal_file = fopen(metal, "r");
  if (metal_file == 0) {
    fprintf(stderr, "Unable to open metal file: %s\n", metal);
    return -1;
  }

  template_file = fopen(template, "r");
  if (template_file == 0) {
    fprintf(stderr, "Unable to open template file: %s\n", template);
    return -1;
  }


  i = 0;

  while(fgets(strbuf, sizeof(strbuf), metal_file)!=NULL){

    sscanf(strbuf, "%*12c %s", atom_name);
    sscanf(strbuf, "%*31c %f %f %f", &xm[i], &ym[i], &zm[i]);

    /*   printf("%f %f %f\n", xm[i], ym[i], zm[i]); */
    ++i;
  }

  no_metals = i;
  /*  printf("%d\n", no_metals);*/

  while (fgets(strbuf, sizeof(strbuf), template_file)!=NULL){  
    sscanf(strbuf, "%*4c %f %f %f", &xtemp, &ytemp, &ztemp);
    unbump = 0;
    j = 0;
    distance = 10.0;
    for (j = 0; j < no_metals; ++j){ 
      distance = (xtemp-xm[j])*(xtemp-xm[j])+
	         (ytemp-ym[j])*(ytemp-ym[j])+
	         (ztemp-zm[j])*(ztemp-zm[j]);

      if (distance < 4.0 ){ 
	unbump = 1;
      }
    }
    if (unbump == 0){
    printf ("%s", strbuf); 
    }     
  }

  /* Clean up and exit */
  fclose(metal_file);
  fclose(template_file);
  return (0);
}
