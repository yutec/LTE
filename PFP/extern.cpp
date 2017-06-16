#include "extern.h"

const int
n_param     =      6,     // # of all input parameters
n_product   =     50,     // # of products
n_period    =    100,     // # of periods for observations
n_inst      =     66,     // # of all instruments
n_theta1    =      5,     // # of linear parameters in utility
n_theta2    =      3,     // # of nonlinear parameter in utility
n_person    =     50,     // # of simulated consumers
n_grid      =     50,     // # of grid points in state space
n_mh        = 100000,     // # of MCMC iteration
n_history   =   2000,     // # of pseudo fixed point solutions stored in history
n_mc        =     20,     // # of Monte Carlo experiments
mc_last     =     20,     // index of the last MC experiment
mc_init     =      1,     // index of the first MC experiment
n_seed      =      5,
n_thin      =     10,
n_theta     = n_theta1 + n_theta2,
n_obs       = n_period * n_product; // # of total sample observations

const int
diagPrint = 0;  // diagnostics printout flag

const float
power = 0.65;

const double
beta = 0.99,
heidel = 1.5,
tol_blp_pfp = 1e-4,
tol_blm_pfp = 1e-6;

int **idx;

const int rcID[]={0,2,3,5}; // Index for variables of random coefficient

const char
directoryInput[1024]="./data/";
