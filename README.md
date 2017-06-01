# LTE

This program runs the second & third Monte Carlo experiments of dynamic BLP
demand estimation using the Pseudo Fixed Point algorithm (Sun & Ishihara, 2017).

Compiling and running the program requires GSL (GNU Scientific Libraray) package,
which can be obtained at "http://www.gnu.org/software/gsl/". Specify the path to
the headers and the link to the GSL library in the Makefile. You will also need GCC or an equivalent compiler, GNU Make, and R installed with Coda & RInside libraries. For NFP & MPEC, Artely's KNITRO is required. In addition, MPEC makes use of ADOL-C and ColPack libraries. 

The Makefile compiles the codes into a command-line console execution file "a.out"
in the Linux & Mac OS. The Makefile uses the GCC compiler while standard
C/C++ compilers should work with no problem.
