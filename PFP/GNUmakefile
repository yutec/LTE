## -*- mode: make; tab-width: 8; -*-
##
## Simple Makefile
##
## TODO: 
##  proper configure for non-Debian file locations,   [ Done ]
##  allow RHOME to be set for non-default R etc

## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)


## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
#RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS := 		$(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) -O3 -march=native -m64
LDLIBS := 		$(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)
CC = g++
Cl = -lm -lgsl -lgslcblas -L/usr/local/lib

a.out: main.o mcmc.o fxpMap.o functions.o tools.o extern.o
	$(CC) -o $@ main.o mcmc.o fxpMap.o functions.o tools.o extern.o $(Cl) $(CPPFLAGS) $(LDLIBS)

main.o: main.cpp
	$(CC) -c main.cpp $(CPPFLAGS) $(CXXFLAGS) 

mcmc.o: mcmc.cpp
	$(CC) -c mcmc.cpp $(CPPFLAGS) $(CXXFLAGS)

fxpMap.o: fxpMap.cpp
	$(CC) -c fxpMap.cpp $(CPPFLAGS) $(CXXFLAGS)

functions.o: functions.cpp
	$(CC) -c functions.cpp $(CPPFLAGS) $(CXXFLAGS) 

tools.o: tools.cpp
	$(CC) -c tools.cpp $(CPPFLAGS) $(CXXFLAGS) 

extern.o: extern.cpp
	$(CC) -c extern.cpp $(CPPFLAGS) $(CXXFLAGS) 

clean:
	-rm -f *.o *.out *.log *~ *.1 outMCMC.txt estimResult.txt
