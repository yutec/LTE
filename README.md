# Pseudo fixed point algorithm for dynamic BLP demand estimation

This program runs the second & third Monte Carlo experiments of Sun & Ishihara (2017) that estimate the dynamic BLP
demand using the Nested Fixed Point, Pseudo Fixed Point, and MPEC algorithms. 

Compiling and running the program requires GSL (GNU Scientific Libraray) package,
which can be obtained at "http://www.gnu.org/software/gsl/". Specify the header path and link to the GSL library in the Makefile. You will also need GCC or an equivalent compiler, GNU Make, and R (https://cran.r-project.org) installed with Coda & RInside libraries. For NFP & MPEC, Artely's KNITRO is required. In addition, MPEC makes use of ADOL-C (https://projects.coin-or.org/ADOL-C) and ColPack (http://www.cs.odu.edu/~dnguyen/dox/colpack/html/). 

The Makefile compiles the codes into a command-line console execution file "a.out"
in the Linux & Mac OS. The Makefile uses the GCC compiler while standard
C/C++ compilers should work with no problem.

A sample of data simulation codes writtin for Matlab is available in the subfolder "data." Users can run it to generate all the necessary data in the txt format for the Monte Carlo exercise. The C++ codes need to be updated with the correct path for the input data files. For MPEC, it is strongly recommended to have more than 20GBs of physical RAM for successful experiment with the large data setting. 
