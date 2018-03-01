/* this program will take a pdb file containing                       *
 * binding site atoms as an input. Will pick two neighboring atoms    *
 * for each atom capable of H-bonding according to a set of rules.    *
 * The output will be a file containing three lines from the pdb      *
 * file of the protein for each H-bonding atom: first line contains   *
 * atom B (the one forming H-bonds), second line contains atom A,     *
 * third line atom C.                                                 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>

#include "defs.h"
#include "types.h"

#include "transform.h"
#include "cube.h"
#include "basics.h"
#include "nbh.h"
#include "read_pdb2.h"
#include "find_select_atoms.h"
#include "create_mathcad_input.h"
#include "create_mathcad_input_sparse.h"
#include "create_mathcad_input_minimal.h"
#include "collapse.h"
#include "quicksort.h"
#include "check_hbond_angle.h"
#include "distance.h"
#include "find_metal_template_points.h"
#include <mymalloc.h>

void  find_hbond_template_points ( char  pdbfilename[256],
                   point_pt     surface_points,
                   int          number_of_surface_points,
                   point_pt     center_points,
                   int          number_of_center_points,
                   point_pt     acceptor_template_points,
                   int          *number_of_acceptor_points,
                   point_pt     donor_template_points,
                   int          *number_of_donor_points,
                   point_pt     doneptor_template_points,
                   int          *number_of_doneptor_points, 
                   int          hbond_pt_density)
{
	pdb_atom_pt atoms;                    /** atoms in the protein **/
	pdb_atom_pt b_atoms;                  /** binding site atoms **/
	pdb_atom_pt r_atoms;                  /** rad atoms **/
	pdb_atom_pt s_atoms;                  /** selected atoms **/
	point_pt metal_template_points;

	int number_of_metal_template_points;

	int i, j, k;
	int       number_of_b_atoms,     /** total number of binding site atoms **/
	          number_of_atoms,       /** total number of atoms **/
	          number_of_r_atoms,     /** total number of rad atoms **/
              number_of_s_atoms;     /** total number of selected atoms **/

	FILE *infile, *mcfile, *selfile;
	char buf[BUFSIZE];
	int n=0;
	double bond_length[MAXSIZE];       /** bond length         **/
	int alp, bet;                      /** temp alpha, beta    **/
	int type[MAXSIZE], label[MAXSIZE]; /** type (A,D,N), label **/
	double alpha[MAXSIZE], beta[MAXSIZE]; /** alpha beta       **/
	double A0[DIM][MAXSIZE], B0[DIM][MAXSIZE], C0[DIM][MAXSIZE];
	double RR[DIM][MAXSIZE];            /** resulting matrix   **/ 
	double distance;
	int    met_index;
	int    too_close;

	double borders[8][3];               /** border points of cube **/
	double s;

	/* int ttt=0;*/

	int ap = 0, dp = 0, np = 0; 
	char temp[16];

	int cpair_cnt=0;                    /** number of pairs **/
	int clt;
	int flag;
	char atm[5];
	char mol[5];

	point_t PP[MAX_HB_TEMP_PTS];  /** initial list of points **/
	point_t temp_pt;
	struct _spoint *SP, *curr, *currj; /** pointer to clustered list **/
	struct _cpair cpair[MAX_HB_TEMP_PTS*MAX_HB_TEMP_PTS]; /** stores array of
                                                   closest pairs of atoms **/

    short bselected[MAX_NUMBER_OF_ATOMS];
    short rselected[MAX_NUMBER_OF_ATOMS];
	double x, y, z;

	/**printf("%d, %d, %d, %d\n", number_of_atoms, number_of_surface_points, number_of_center_points, hbond_pt_density);
		for(j=0; j<number_of_center_points; j++)
		printf("%.3f %.3f %.3f\n", center_points[j].pos[0], center_points[j].pos[1], center_points[j].pos[2]); 
	**/

	atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
	b_atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
	s_atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
	r_atoms = (pdb_atom_pt) mymalloc ( MAX_NUMBER_OF_ATOMS * sizeof (pdb_atom_t) );
	metal_template_points = (point_pt) mymalloc ( MAX_NUMBER_OF_TEMPLATE_POINTS * sizeof (point_t) );

	number_of_metal_template_points = 0;

	number_of_atoms=0; number_of_b_atoms=0; number_of_r_atoms=0;
	number_of_atoms = read_pdb2( pdbfilename, atoms);

	/**  for (i = 0;i<number_of_atoms;i++) 
	   		 printf("%s %s %s %d\n",atoms[i].residue_name, atoms[i].residue_number, atoms[i].atom_name, atoms[i].act);
	**/

   
	for(j=0;j<number_of_atoms;j++)
	{ bselected[j]=0; rselected[j]=0; }

	for(i=0;i<number_of_surface_points;i++)
	{
		if(surface_points[i].type==OK)
		{
			x=surface_points[i].pos[0];
			y=surface_points[i].pos[1];
			z=surface_points[i].pos[2];
			/**printf("%.2f %.2f %.2f\n", x, y, z); **/

			get_atoms( x, y, z, atoms, number_of_atoms, b_atoms, &number_of_b_atoms, r_atoms, &number_of_r_atoms, bselected, rselected); 
		}
	}
	/*	printf("%d, %d, %d, %d\n", number_of_b_atoms, number_of_r_atoms, number_of_atoms, number_of_s_atoms );*/


	/**for (i = 0;i<number_of_b_atoms;i++) {
	   		 printf("%s, %s, %s, %c, %s, %c, ", b_atoms[i].atom_number, b_atoms[i].atom_name, b_atoms[i].residue_name, b_atoms[i].chainID, b_atoms[i].residue_number, b_atoms[i].icode);
	    	printf("%.3f %.3f %.3f, ",b_atoms[i].pos[0], b_atoms[i].pos[1], b_atoms[i].pos[2]);
	   		 printf("%d %d %d %d %d\n",b_atoms[i].type, b_atoms[i].residue_type, b_atoms[i].hphil, b_atoms[i].act, b_atoms[i].class);
	  	}
	**/

	/** select atoms that can hydrogen bond **/
	find_select_atoms(atoms, number_of_atoms, b_atoms, number_of_b_atoms, s_atoms, &number_of_s_atoms);
	/*printf("%d, %d, %d, %d\n", number_of_b_atoms, number_of_r_atoms, number_of_atoms, number_of_s_atoms );*/
	find_metal_template_points (b_atoms, number_of_b_atoms, metal_template_points, &number_of_metal_template_points);

#ifdef DEBUG
	/** create intermediate file for debugging **/
	if((selfile=fopen( SEL_ATOMS, "w"))==NULL)
	{ 
		fprintf(stderr,"ERROR: Cannot find %s\n", SEL_ATOMS); 
		exit(1); 
	}

	for (i = 0;i<number_of_s_atoms;i++)
	  {
	    if(s_atoms[i].type==UNKNOWN)
	      fprintf(selfile, "ATOM  %5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f\n", s_atoms[i].atom_number, s_atoms[i].atom_name, s_atoms[i].residue_name, s_atoms[i].chainID, s_atoms[i].residue_number, s_atoms[i].icode, s_atoms[i].pos[0], s_atoms[i].pos[1], s_atoms[i].pos[2]);
	    else if(s_atoms[i].type==HETERO || s_atoms[i].act == METAL_1 || s_atoms[i].act == METAL_2 )
	      fprintf(selfile, "HETERO%5s %4s %3s %1c%4s%1c   %8.3f%8.3f%8.3f\n", s_atoms[i].atom_number, s_atoms[i].atom_name, s_atoms[i].residue_name, s_atoms[i].chainID, s_atoms[i].residue_number, s_atoms[i].icode, s_atoms[i].pos[0], s_atoms[i].pos[1], s_atoms[i].pos[2]);
	  }
	fclose(selfile);
#endif

	/** creates another intermediate file MCDINFILE which is used in the next
        step **/
    if(hbond_pt_density==DENSE)
		create_mathcad_input(s_atoms, number_of_s_atoms);
    else if(hbond_pt_density==SPARSE)
		create_mathcad_input_sparse(s_atoms, number_of_s_atoms);
    else if(hbond_pt_density==MINIMAL)
		create_mathcad_input_minimal(s_atoms, number_of_s_atoms);

	/** reads the intermediate file created in the last step **/
	if((infile=fopen(MCDINFILE, "r"))==NULL)
	{ 
		fprintf(stderr,"ERROR: Cannot find %s\n", MCDINFILE); 
		exit(1); 
	}

	while(fgets(buf, BUFSIZE, infile)!=NULL)
	{
		if((++n)>MAXSIZE)
		{ 
			fprintf(stderr, "ERROR: Max size of possible selected points exceeded\n"); 
			exit(1); 
		}

		/** A0, B0, C0 initial coords of atoms A, B and C **/
		sscanf( buf, "%d%d%d%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf", &type[n-1], &alp, &bet, &bond_length[n-1], &label[n-1], &A0[0][n-1], &A0[1][n-1], &A0[2][n-1], &B0[0][n-1], &B0[1][n-1], &B0[2][n-1], &C0[0][n-1], &C0[1][n-1], &C0[2][n-1]);
		alpha[n-1]=(alp * PI)/180.0;
		beta[n-1]=(bet * PI)/180.0;

		/**printf("%d %d %d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", type[n-1], alpha[n-1], beta[n-1], bond
		    _length[n-1], label[n-1], A0[0][n-1], A0[1][n-1], A0[2][n-1], B0[0][n-1], B0[1][n-1], B0[2][n-1], C0[0][n
		    -1], C0[1][n-1], C0[2][n-1]);**/

		/** transform matrix **/
		mat_trans(RR, A0, B0, C0, alpha[n-1], beta[n-1], bond_length[n-1], n);
	}
	fclose(infile);
	/** delete intermediate file **/
	unlink(MCDINFILE);

#ifdef DEBUG
	/** create intermediate debugging output **/
	if((mcfile=fopen( MCDOUTFILE, "w"))==NULL)
	{ 
		fprintf(stderr,"ERROR: Cannot find %s\n", MCDOUTFILE); 
		exit(1); 
	}

	for(j=0;j<n;j++)find_select_atom
	{
		if(type[j]==1) fprintf(mcfile, "A* ");
		else if(type[j]==2) fprintf(mcfile, "D* ");
		else if(type[j]==3) fprintf(mcfile, "N* ");
		fprintf(mcfile, "%.3f %.3f %.3f\n", RR[0][j], RR[1][j], RR[2][j]);
	}
	fclose(mcfile);
#endif

	strcpy(atm,"O"); 
	strcpy(mol,"HOH");
	curr=SP=NULL;

	for(j=0; j<number_of_center_points; j++)
	{
		borders[j][0] = center_points[j].pos[0];	
		borders[j][1] = center_points[j].pos[1];	
		borders[j][2] = center_points[j].pos[2];
	}

	k=0;
	for(j=0;j<n;j++)
	{
		/** only output points if they are inside the box and dont bump into
		other protein molecules **/
		if((point_in_cube(borders, RR[0][j], RR[1][j], RR[2][j])==1) &&
		    (check_bump(RR[0][j], RR[1][j], RR[2][j], r_atoms, number_of_r_atoms, NBH_THRESHOLD)==-1))
		{
			if((k+1)>MAX_HB_TEMP_PTS)
			{ 
				fprintf(stderr, "ERROR: Max number of HB template points exceeded\n"); 
				exit(1); 
			}

			/**if(type[j]==1) printf("A* ");
				else if(type[j]==2) printf("D* ");
				else if(type[j]==3) printf("N* ");
				printf("%.3f %.3f %.3f\n", RR[0][j], RR[1][j], RR[2][j]);**/

			if(type[j]==1) PP[k].type = ACCEPTOR;
			else if(type[j]==2)  PP[k].type = DONOR;
			else if(type[j]==3) PP[k].type = DONEPTOR;

			PP[k].pos[0]=RR[0][j];
			PP[k].pos[1]=RR[1][j];
			PP[k].pos[2]=RR[2][j];


			/** fill up the cluster data structure
				initally each point is considered to be
				a single cluster **/

			if(!curr)
			{
				SP=(struct _spoint *)malloc(sizeof(struct _spoint));
				curr=SP;
			}
			else
			{
				curr->next=(struct _spoint *)malloc(sizeof(struct _spoint));
				curr=curr->next;
			}
			curr->next=NULL;
			curr->tpnt.pos[0] = PP[k].pos[0];
			curr->tpnt.pos[1] = PP[k].pos[1];
			curr->tpnt.pos[2] = PP[k].pos[2];
			curr->tpnt.type = PP[k].type;
			curr->oindx[0] = k;
			curr->nindx = 1;
			k++;
		}
	}

	/*printf("Total = %d\n",n=k );*/
	/** change n to new k **/
	n=k;

	/*for(curr=SP; curr!=NULL; curr=curr->next)
	{
		if(curr->tpnt.type==ACCEPTOR) printf("A* ");
		else if(curr->tpnt.type==DONOR) printf("D* ");
		else if(curr->tpnt.type==DONEPTOR) printf("N* ");

		printf("%8.3f %8.3f %8.3f [ ", curr->tpnt.pos[0], curr->tpnt.pos[1], curr->tpnt.pos[2]);
		for(i=0; i<	curr->nindx; i++)
			printf("%d ", curr->oindx[i]);
		printf("]\n");
	}*/

	flag = 4; /** start by accepted as default **/

	/** each iteration reduces the number of clusters by one **/
	while(flag!=5)
	{
		
		/** fill in n^2 closest pair cpair structure... **/
		cpair_cnt=0;
		for(curr=SP; curr!=NULL; curr=curr->next)
		{
			for(currj=curr->next; currj!=NULL; currj=currj->next)
			{
				if(curr!=currj)
				{
					cpair[cpair_cnt].pnts[0]=curr;
					cpair[cpair_cnt].pnts[1]=currj;
					cpair[cpair_cnt].dist = find_distance(curr->tpnt, currj->tpnt);
					cpair_cnt++;
				}
			}
		}

		/** we use sort program of our own and sort **/
		/*q_sort(cpair, cpair_cnt, sizeof(struct _cpair), &comp_cpair);*/
		quicksort_cpair(cpair, 0, cpair_cnt-1);


		/** look for the closest pair that can be clustered **/
		for(i=0; i<cpair_cnt; i++)
		{
			/** cannot clusters points (clusters) that are far away **/
			if(cpair[i].dist<HBOND_THRESHOLD)
			{

				/** printf("%d: %.2f %.2f, %.2f %.2f, %.2f %.2f : %f\n", i, cpair[i].pnts[0]->tpnt.pos[0], cpair[i].pnts[1]->tpnt.pos[0], cpair[i].pnts[0]->tpnt.pos[1], cpair[i].pnts[1]->tpnt.pos[1], cpair[i].pnts[0]->tpnt.pos[2], cpair[i].pnts[1]->tpnt.pos[2], cpair[i].dist);
				printf("( [ ");
				for(j=0;j<cpair[i].pnts[0]->nindx;j++) printf("%d ", cpair[i].pnts[0]->oindx[j]);
				printf("] [ ");
				for(j=0;j<cpair[i].pnts[1]->nindx;j++) printf("%d ", cpair[i].pnts[1]->oindx[j]);
				printf("] ) ");**/

				/** check for proper type **/
				if( cpair[i].pnts[0]->tpnt.type == DONEPTOR || 
				    cpair[i].pnts[1]->tpnt.type == DONEPTOR || 
				    ( cpair[i].pnts[0]->tpnt.type == ACCEPTOR && cpair[i].pnts[1]->tpnt.type == DONOR) ||
				    ( cpair[i].pnts[1]->tpnt.type == ACCEPTOR && cpair[i].pnts[0]->tpnt.type == DONOR) ||
				    ( cpair[i].pnts[1]->tpnt.type == ACCEPTOR && cpair[i].pnts[0]->tpnt.type == ACCEPTOR) ||
				    ( cpair[i].pnts[1]->tpnt.type == DONOR && cpair[i].pnts[0]->tpnt.type == DONOR))
				{
					/** check if the new point would not be too far from 
									others in the cluster **/
					if( valid_pt(cpair[i].pnts[0], cpair[i].pnts[1], PP, &temp_pt))
					{
						flag = 4;
						/**printf("ACCEPTED\n");**/
						break;
					} else { 
						flag = 3; 
						/**printf("Rejected: Large drift\n"); **/
					}
				} else { 
					flag = 2; 
					/**printf("Rejected: Type mismatch\n"); **/
				}
			} else { 
				flag = 1; 
			}
		}

		/** .. and finally cluster **/
		if(flag==4) /* ith pair accepted for clustering **/
		{

			/** find new mean point. change coordinates of Point 1 and Point 2 
			to the new mean point **/
			cpair[i].pnts[0]->tpnt.pos[0] = temp_pt.pos[0];
			cpair[i].pnts[0]->tpnt.pos[1] = temp_pt.pos[1];
			cpair[i].pnts[0]->tpnt.pos[2] = temp_pt.pos[2];

			/** set the proper type for the new cluster **/
			if(cpair[i].pnts[0]->tpnt.type != cpair[i].pnts[1]->tpnt.type)
				cpair[i].pnts[0]->tpnt.type = cpair[i].pnts[1]->tpnt.type = DONEPTOR;

			/** merge the list of points in each cluster **/
			cmerge( cpair[i].pnts[0]->oindx, cpair[i].pnts[1]->oindx, &(cpair[i].pnts[0]->nindx), &(cpair[i].pnts[1]->nindx));

			/** delete the second point(cluster) from the list **/
			n--;
			if(cpair[i].pnts[1] == SP)
			{
				SP = SP->next;
				free(cpair[i].pnts[1]);
			}
			else
			{
				for(curr=SP; curr->next!=cpair[i].pnts[1]; curr=curr->next);
				curr->next = curr->next->next;
				free(cpair[i].pnts[1]);
			}
			/** down to single cluster, so set flag for exit*/
			if(n==1) {
				flag=5; 
				SP->next=NULL;
			}
		}

		/** looked through all clusters and didnt find one, 
			so set flag for exit **/
		if(i==cpair_cnt) flag=5;

	} /** end while loop **/

	/**printf("Total = %d\n", n);**/
	clt=10;

	ap = dp = np = 0;
	for(curr=SP; curr!=NULL; curr=curr->next)
	{
		x = curr->tpnt.pos[0];
		y = curr->tpnt.pos[1];
		z = curr->tpnt.pos[2];
		too_close = FALSE;
		for (met_index = 0; met_index < number_of_metal_template_points; met_index++ )
		  {
		    distance = sqrt((x-metal_template_points[met_index].pos[X])*(x-metal_template_points[met_index].pos[X])+(y-metal_template_points[met_index].pos[Y])*(y-metal_template_points[met_index].pos[Y])+(z-metal_template_points[met_index].pos[Z])*(z-metal_template_points[met_index].pos[Z]));
		    /*		    printf ("distance = %f\n", distance);*/
		    if (distance < HBOND_THRESHOLD )
		      {
			/*printf ("too close distance = %f\n", distance);*/
			too_close = TRUE;
			break;
		      }
		  }
				
		if(curr->tpnt.type==ACCEPTOR) 
		{
			/**printf("A* "); toneroma 07JUN07 - want to keep acceptor points even if they are close to metals**/
			acceptor_template_points[ap].pos[0]=x;
			acceptor_template_points[ap].pos[1]=y;
			acceptor_template_points[ap].pos[2]=z;
			acceptor_template_points[ap].type=ACCEPTOR;
			ap++;
		}
		else if(curr->tpnt.type==DONOR && too_close == FALSE)
		{
			/**printf("D* "); toneroma 07JUN07 - want to remove donor points that are too close to metals**/
			donor_template_points[dp].pos[0]=x;
			donor_template_points[dp].pos[1]=y;
			donor_template_points[dp].pos[2]=z;
			donor_template_points[dp].type=DONOR;
			dp++;
		}
		else if(curr->tpnt.type==DONEPTOR && too_close == FALSE)
		{
			/**printf("N* "); toneroma 07JUN07 - want to keep acceptor points that are close to metals, so remove donor part of doneptor**/
			doneptor_template_points[np].pos[0]=x;
			doneptor_template_points[np].pos[1]=y;
			doneptor_template_points[np].pos[2]=z;
			doneptor_template_points[np].type=DONEPTOR;
			np++;
		}
		else if(curr->tpnt.type==DONEPTOR && too_close == TRUE)
		{
			/**printf("A* "); toneroma 07JUN07 - want to keep acceptor points that are close to metals, so remove donor part of doneptor**/
			acceptor_template_points[ap].pos[0]=x;
			acceptor_template_points[ap].pos[1]=y;
			acceptor_template_points[ap].pos[2]=z;
			acceptor_template_points[ap].type=ACCEPTOR;
			ap++;
		}

	}
	
	/* set metal template points as acceptor points - toneroma 07JUN07 */
		for (met_index = 0; met_index < number_of_metal_template_points; met_index++ )
	  {
	    acceptor_template_points[ap].pos[X]=metal_template_points[met_index].pos[X];
	    acceptor_template_points[ap].pos[Y]=metal_template_points[met_index].pos[Y];
	    acceptor_template_points[ap].pos[Z]=metal_template_points[met_index].pos[Z];
	    acceptor_template_points[ap].type=ACCEPTOR;
	    ap++;
	  }	 
	
	(*number_of_acceptor_points) = ap;
	(*number_of_donor_points) = dp;
	(*number_of_doneptor_points) = np;

}

int   check_hbond_point ( pdb_atom_pt  atoms,
			  int          number_of_atoms,
			  point_pt     point )
{
  double dist;
  int    i;
  int    retval = NOTHING;
  
  for ( i = 0; i < number_of_atoms; i++ )
    if ( atoms[i].act != NOTHING )
      {
	dist = distance ( point->pos,
			  atoms[i].pos );
	if ( dist > MIN_HBOND_LENGTH && dist < MAX_HBOND_LENGTH )
	  {
	    if ( ( point->type == ACCEPTOR
		   && ( atoms[i].act == ACCEPTOR
			|| atoms[i].act == DONEPTOR ) )
		 || ( point->type == DONOR
		      && ( atoms[i].act == DONOR
			   && ( atoms[i].type == HETERO || atoms[i].act == METAL_1 || atoms[i].act == METAL_2 
				|| check_hbond_angle ( atoms,
						       number_of_atoms,
						       i,
						       point ) ) ) ) )
	      retval = DONEPTOR;
	    else
	      switch ( atoms[i].act )
		{
		case DONOR:
		  if ( atoms[i].type == HETERO || atoms[i].act == METAL_1 || atoms[i].act == METAL_2
		       || check_hbond_angle ( atoms,
					      number_of_atoms,
					      i,
					      point ) )
		    retval = ACCEPTOR;
		  break;
		case ACCEPTOR:
		  retval = DONOR;
		  break;
		case DONEPTOR:
		  if ( atoms[i].type == HETERO || atoms[i].act == METAL_1 || atoms[i].act == METAL_2
		       || check_hbond_angle ( atoms,
					      number_of_atoms,
					      i,
					      point ) )
		    retval = DONEPTOR;
		  else
		    retval = DONOR;
		  break;
		default:
		  err_panic ("check_hbond_point",
			     "unknown atom act" );
		  break;
		}
	  }
      }
  point->type = retval;
  return point->type;
}

