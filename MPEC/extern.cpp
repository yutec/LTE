#include "extern.h"

const int
n_product   =     50,     // # of products
n_period    =    100,     // # of periods for observations
n_inst      =     66,     // # of all instruments
n_theta1    =      5,     // # of linear parameters in utility
n_theta2    =      3,     // # of nonlinear parameter in utility
n_person    =     50,     // # of simulated consumers
n_mc        =     20,     // # of total Monte Carlo experiments
mc_last     =      5,     // index of MC experiment after the last
mc_init     =      0,     // # of MC experiments to be skipped
n_seed      =      5,
n_theta     = n_theta1 + n_theta2,
n_obs       = n_period * n_product, // # of total sample observations
dimValfun   = n_person * n_grid,
dimOptim    = n_theta + n_obs + dimValfun + n_inst,
dimConst    = n_obs + dimValfun + n_inst;

const int diagPrint = 0;  // diagnostics printout flag

const double
beta = 0.99,
tol_blp = 1e-12,
tol_blm = 1e-14;

int **idx;

const int rcID[]={1,2,4}; // Index for variables of random coefficient