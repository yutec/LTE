extern const int
n_param,    // # of all parameters in param_input.txt
n_product,  // # of products
n_period,   // # of periods for observations
n_inst,     // # of total instruments
n_theta1,   // # of linear parameters in utility
n_theta2,   // # of nonlinear parameter in utility
n_person,   // # of simulated consumers
n_grid,     // # of grid points in state space
n_mc ,      // # of Monte Carlo experiments
mc_last,
mc_init,    // index of the first MC experiment to start
n_seed,
n_theta,
n_obs;      // # of total sample observations

extern const int diagPrint;

extern const double
beta,
tol_blp,  // convergence threshold for BLP fixed point operator
tol_blm;  // convergence threshold for Bellman fixed point operator

extern int
**idx;  // index mapping time t & product j into vector of length n_obs

extern const int
rcID[4];
