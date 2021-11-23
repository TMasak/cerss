library.R contains functions implementing all the algorithms from the main paper or the supplement.
demo.R shows how the functions can be used in form of a sample simulation run.
mortality_analysis.R can be used to recreate all the results of the real data analysis
testing.R performs the bootstrap separability test of Aston et al. (2017) and its version, testing
    the proposed separable-plus-banded model for the real data

All the results in the simulation sections of the paper were calculated using a cluster. The scripts are
provided in the folder 'scripts'.
One can manually check that a selected datum is really an output of a given script by running the script
with a the parameters set correspondingly.
However, running all the simulations on a single machine would take thousands of hours (a large portion
of that is due to the empirical covariance estimator calculation).
- script_norms.R is used to calculate the 'bias' curves for all the plots separately

