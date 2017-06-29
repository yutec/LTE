# The pseudo fixed point (PFP) algorithm

This program runs the second & third Monte Carlo experiments of Sun & Ishihara (2017) (https://ssrn.com/abstract=2140617) that estimate the dynamic BLP
demand using the Nested Fixed Point, Pseudo Fixed Point, and MPEC algorithms. 

Compiling and running the program requires GSL (GNU Scientific Libraray) package (http://www.gnu.org/software/gsl). Specify the header path and link to the GSL library in the Makefile. You will also need GCC or an equivalent compiler, GNU Make, and R (https://cran.r-project.org) installed with Coda & RInside libraries. For NFP & MPEC, Artely's KNITRO is required (https://www.artelys.com). In addition, MPEC makes use of ADOL-C (https://projects.coin-or.org/ADOL-C) and ColPack (http://www.cs.odu.edu/~dnguyen/dox/colpack/html) for automatic differentiation and sparsity detection, respectively. 

The Makefile compiles the codes into a command-line console execution file "a.out"
in the Linux & Mac OS. The Makefile uses the GCC compiler while standard
C/C++ compilers should work with no problem.

A sample of data simulation codes written for Matlab is available in the subfolder "data." Users can run it to generate all the necessary data in the txt format for the Monte Carlo exercise. The C++ codes (extern.cpp or main.cpp) need to be updated with the correct path for the input data files as well as the global constants defined therein. For MPEC, it is strongly recommended to have more than 20GBs of physical RAM for successful experiment in the large data setting. 

Information on the authors can be found at http://people.stern.nyu.edu/mishihar/ and https://sites.google.com/site/yutecsun.
