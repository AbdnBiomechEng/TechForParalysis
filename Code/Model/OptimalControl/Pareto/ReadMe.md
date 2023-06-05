# Running code for optimal control

Include TechForParalysis\Code\Model and TechForParalysis\Code\Model\OptimalControl in your Matlab path.

## Elbow kinematics

To optimise elbow trajectory and muscle activations.

1. First run find_das3elbow_feas_solutions.m to find a range of feasible solutions. This 'runs' an optimisation with the cost function set to zero just to generate a set of feasible solutions for an initial population.
2. Run create_initial_population_elbow.m selects the successful outputs from the previous step and saves them in a matrix (mat file init_pop_elbow_pareto).
3. Run das3elbow_optimize_pareto.m by entering ```result = das3elbow_optimize_pareto("output_filename");``` at the command line

