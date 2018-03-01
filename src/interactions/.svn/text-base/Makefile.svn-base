#
# Makefile for compute_interaction_centers, check_connectivity,
#    split_mol2, and generate_rasmol_script     
#
#                        Matt Tonero      Tue Feb 06, 2007
#
#                        Volker Schnecke  Thu Sep 16 11:08:50 EDT 1999
#
# Before using this makefile, the environment variable $SLIDE_DIR has to be
# correctly set to the root directory of the global SLIDE installation.
# BINPATH is the path to directory that includes the binaries of the
# program and all Perl scripts that come with SLIDE
#BINPATH = $(SLIDE_DIR)/bin
BINPATH = ../../bin

#DEMO_CODE = check_key.o
 
CC = gcc
COPTS = -O2 #-m32
#COPTS = -g
INCLUDES = -Iinc 
LIBS = -lm
CFLAGS = $(COPTS) $(INCLUDES)

OBJS = 	read_mol2.o read_flex_defn.o err_handle.o mymalloc.o \
	adj_list.o find_flexible_bonds.o find_cycles.o read_hyd_defn.o \
	read_rule.o assign_type.o check_rule.o find_hyd_atoms.o \
	assign_hydrogens.o bitstrings.o find_carbon_centers.o \
	find_carbon_ring_centers.o check_connectivity.o split_molecule.o \
        complete_link_clustering.o distance.o basics.o refine.o $(DEMO_CODE)

default: $(BINPATH)/generate_rasmol_script $(BINPATH)/split_mol2 \
	$(BINPATH)/check_connectivity $(BINPATH)/compute_interaction_centers

$(BINPATH)/compute_interaction_centers:	 $(OBJS) $(SRCS) main_interactions.c
	$(CC) $(CFLAGS) -o $(BINPATH)/compute_interaction_centers $(OBJS) \
		main_interactions.c $(LIBS)

$(BINPATH)/check_connectivity:	 $(OBJS) $(SRCS) main_connectivity.c
	$(CC) $(CFLAGS) -o $(BINPATH)/check_connectivity $(OBJS) \
		main_connectivity.c $(LIBS)

$(BINPATH)/split_mol2:	 $(OBJS) $(SRCS) main_split.c
	$(CC) $(CFLAGS) -o $(BINPATH)/split_mol2 $(OBJS) \
		main_split.c $(LIBS)

$(BINPATH)/generate_rasmol_script:	 $(OBJS) $(SRCS) main_rasmol_script.c
	$(CC) $(CFLAGS) -o $(BINPATH)/generate_rasmol_script $(OBJS) \
		main_rasmol_script.c $(LIBS)

.c.o: 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -f *.o $(BINPATH)/compute_interaction_centers \
		$(BINPATH)/check_connectivity $(BINPATH)/generate_rasmol_script
