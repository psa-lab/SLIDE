#
# Makefile for utils   Jeff Van Voorst     Thu Apr 30

CC = gcc
COPTS = -O2 # Use for distribution
#COPTS = -O2 -g -pg # Use for profiling with gprof
#COPTS = -O2 -g -pg -fprofile-arcs -ftest-coverage # Use for profiling with gprof and gcov
#COPTS = -g -Wall -Wundef -Wshadow # use for debugging

INCLUDES = -Iinc
# Set the SLIDE version string here -- note the \s are needed with the quotes
# to denote VERSION as a string
VERSION = -DVERSION=\"3.4\"
CFLAGS = $(COPTS) $(INCLUDES) $(VERSION)
LIBS =  -lm

OBJS =	basics.o slide_itimer.o err_handle.o

libslideutils.a: $(OBJS)
	rm -f libslideutils.a && \
        ar -r libslideutils.a *.o && \
	ranlib libslideutils.a

.c.o: 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -f *.a *.o *.gcda *.gcno 

touch:  
	touch *.c


