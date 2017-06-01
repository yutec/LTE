#include "extern.h"

const int
n_product   =     50,     // # of products
n_period    =    100,     // # of periods for observations
n_inst      =     66,   // # of all instruments
n_theta1    =      5,     // # of linear parameters in utility
n_theta2    =      3,     // # of nonlinear parameter in utility
n_person    =     50,     // # of simulated consumers
n_grid      =     50,     // # of grid points in state space
n_mc        =     20,     // # of Monte Carlo experiments
mc_last     =      5,
mc_init     =      1,
n_seed      =      5,
n_theta     = n_theta1 + n_theta2,
n_obs       = n_period * n_product; // # of total sample observations

const int diagPrint = 0;  // diagnostics printout flag

const double
beta = 0.99,
tol_blp = 1e-12,
tol_blm = 1e-14;

int **idx;

const int rcID[]={0,2,3,5}; // Index for variables of random coefficient