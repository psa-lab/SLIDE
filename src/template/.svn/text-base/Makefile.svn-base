#
# Makefile for unbiased_template and ligand_based_template
#                     Matt Tonero         Tue Feb 06, 2007
#
#                     Volker Schnecke     Wed Sep 15 14:46:44 EDT 1999
#
# Before using this makefile, the environment variable $SLIDE_DIR has to be
# correctly set to the root directory of the global SLIDE installation.
 
# BINPATH is the path to directory that includes the binaries of the
# program and all Perl scripts that come with SLIDE

#BINPATH = $(SLIDE_DIR)/bin
BINPATH = ../../bin
#DEMO_CODE = check_key.o

# Set the SLIDE version string here -- note the \s are needed with the quotes
# to denote VERSION as a string
VERSION = -DVERSION=\"3.4\"
 
CC = gcc
#COPTS = -O2 -m32 -gstabs+ #-g
COPTS = -O2 #-Wall#-m32 #-Wall #-gstabs+ #-g
#COPTS = -Wall -g
INCLUDES = -Iinc -I../utils/inc
CFLAGS = $(COPTS) $(INCLUDES) $(VERSION)
LIBS = ../utils/libslideutils.a -lm

OBJS =  read_pdb_file.o read_borders_file.o \
	distance.o find_hbond_template_points.o assign_hphil.o \
	find_hphob_template_points.o generate_points.o filter_points.o \
	complete_link_clustering.o check_hbond_angle.o \
	find_hyd_atoms.o read_hyd_defn.o read_flex_defn.o read_mol2.o \
	bitstrings.o assign_type.o check_rule.o find_carbon_ring_centers.o \
	adj_list.o assign_hydrogens.o read_rule.o \
	math_util.o transform_points.o get_matrix.o assign_atom_class.o \
	create_mathcad_input.o create_mathcad_input_minimal.o \
	create_mathcad_input_sparse.o  cube.o find_select_atoms.o \
	nbh.o read_pdb2.o transform.o collapse.o quicksort.o \
	find_carbon_centers.o find_metal_template_points.c $(DEMO_CODE)

default: $(BINPATH)/ligand_based_template $(BINPATH)/template $(BINPATH)/gen_lig_box $(BINPATH)/metal_template $(BINPATH)/remove_close_points

$(BINPATH)/template: $(OBJS) main_template.c
	$(CC) $(CFLAGS) -o $(BINPATH)/template \
		main_template.c $(OBJS) $(LIBS)

$(BINPATH)/ligand_based_template: $(OBJS) main_ligand_based.c
	$(CC) $(CFLAGS) -o $(BINPATH)/ligand_based_template \
		main_ligand_based.c $(OBJS) $(LIBS)

$(BINPATH)/gen_lig_box: $(OBJS) gen_lig_box.c
	$(CC) $(CFLAGS) -o $(BINPATH)/gen_lig_box read_mol2.o \
		assign_type.o bitstrings.o gen_lig_box.c $(LIBS)

#
# 2005:03:03 
#
# TO FACILITATE INSTALLATION OF metal_template AND remove_close_points
#

$(BINPATH)/metal_template: metal_template.o
	$(CC) $(CFLAGS) -o $(BINPATH)/metal_template metal_template.c $(LIBS)

$(BINPATH)/remove_close_points: remove_close_points.o
	$(CC) $(CFLAGS) -o $(BINPATH)/remove_close_points \
			remove_close_points.c $(LIBS)

.c.o: 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -f *.o core

veryclean:
	rm -f *.o $(BINPATH)/unbiased_directed_template \
         $(BINPATH)/ligand_based_template \
	 $(BINPATH)/random_template \
	 $(BINPATH)/template core

touch:  
	touch *.c


