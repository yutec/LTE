######################################################################
#
#  Copyright (c) 2006-2014 by Ziena Optimization LLC
#
#  Makefile for Unix platforms.
#
#  Edit the macro for OPTPROB to solve one of the example problems.
#    To build on Linux:    gmake
#    To build on MacOSX:   gnumake
#
#  This makefile builds both statically linked and dynamically linked
#  versions of each example.  The dynamic versions (and also the static 
#  versions for MacOSX) cannot execute unless
#  ../../lib is added to $LD_LIBRARY_PATH (on MacOSX use the environment
#  variable $DYLD_LIBRARY_PATH).
#    bash shells:  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib
#    tcsh shells:  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:../../lib
#
######################################################################

# Set the location of KNITRO.
KNDIR = /usr/local/knitro/
KNRELEASE = 1000

# Set up platform-specific make parameters.
UNAME = $(shell uname -s)

# These are parameters for Linux platforms.
ifeq ($(findstring Linux,${UNAME}), Linux)
  CC = gcc 
  LD = g++ 
  LIBS = -ldl -lm -lgsl -lgslcblas
  KNLIB_STATIC  = $(KNDIR)/lib/libknitro$(KNRELEASE).a
  KNLIB_DYNAMIC = $(KNDIR)/lib/libknitro.so

  CFLAGS = -c -O3 -march=native -std=c++11

  # Try to detect if this is a 64-bit platform.
  MNAME = $(shell uname -m)
  ifeq ($(findstring x86_64,${MNAME}), x86_64)
    CC = gcc -fopenmp 
    LD = g++ -fopenmp
    CFLAGS = -c -O3 -march=native -m64 -std=c++11
  endif
endif

# These are parameters for MacOSX platforms.
ifeq ($(findstring Darwin,${UNAME}), Darwin)
# gcc compiler options
#  CC = gcc -arch x86_64 -fopenmp 
#  LD = g++ -arch x86_64 
# clang compiler options
  CC = gcc -arch x86_64 -stdlib=libstdc++
  LD = g++ -arch x86_64 -stdlib=libstdc++
  LIBS = $(KNDIR)/lib/libiomp5.dylib -ldl 
  KNLIB_STATIC  = $(KNDIR)/lib/libknitro$(KNRELEASE).a
  KNLIB_DYNAMIC = $(KNDIR)/lib/libknitro.dylib
  CFLAGS = -c -O
endif

######################################################################

a.out: main.o nfpIteration.o functions.o tools.o extern.o $(KNLIB_STATIC)
	$(LD) -o $@ $^ $(LIBS)

main.o: main.cpp
	$(CC) $(CFLAGS) -I$(KNDIR)/include $<

nfpIteration.o: nfpIteration.cpp
	$(CC) $(CFLAGS) -I$(KNDIR)/include $<

functions.o: functions.cpp
	$(CC) $(CFLAGS) -I$(KNDIR)/include $<

tools.o: tools.cpp
	$(CC) $(CFLAGS) -I$(KNDIR)/include $<

extern.o: extern.cpp
	$(CC) $(CFLAGS) -I$(KNDIR)/include $<

clean:
	-rm -f *.o *.out *~ *.1 

