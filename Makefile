#
# Top-Level Makefile for SLIDE   
# John B. Johnston         
# Thurs June 02, 2011
#
# Before making SLIDE, the environment variable $SLIDE_DIR has to be
# correctly set to the root directory of the global SLIDE installation.

all:
	perl ./install_slide.pl
clean:
	-rm bin/*
	cd src/interactions && make clean
	cd src/slide && make clean
	cd src/template && make clean
	cd src/utils && make clean


