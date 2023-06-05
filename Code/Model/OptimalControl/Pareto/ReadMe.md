# Running code for optimal control

Include TechForParalysis\Code\Model in your Matlab path.

-If you are using a Matlab version prior to 2020, you should also add **TechForParalysis\Code\Model\OptimalControl** to the Matlab path to use the IPOPT version included.

If you are using Matlab after 2020, you will need to install the version of IPOPT given as a toolbox here: https://github.com/ebertolazzi/mexIPOPT. This can also be installed directly from Matlab using the Add-on Explorer.

## Elbow kinematics

To optimise elbow trajectory and muscle activations.

1. First run find_das3elbow_feas_solutions.m to find a range of feasible solutions. This 'runs' an optimisation with the cost function set to zero just to generate a set of feasible solutions for an initial population.
2. Run create_initial_population_elbow.m to select the successful outputs from the previous step and save them in a matrix (mat file init_pop_elbow_pareto). This also plots histograms and a scatterplot of the two cost functions.
3. Run das3elbow_optimize_pareto.m by entering ```result = das3elbow_optimize_pareto('output_filename');``` at the command line
