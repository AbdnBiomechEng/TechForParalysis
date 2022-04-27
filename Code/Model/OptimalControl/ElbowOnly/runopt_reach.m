function runopt_reach    
	% Run a sequence of optimizations, with mesh refinement
    
    % Add seed for reproducibility
    rng(123);

    factor = 2;    
    maxnodes = 40;		  % end close to this number of nodes
    nodes = 5;            % start with this number of nodes

    % Which model to use
    OptSetup.model_file = 'model_struct_A.mat'; 
    
    % Add external force at the centre of mass of the hand 
    % (defined in the global frame: +X is laterally, +Y is upwards and +Z is posteriorly)
    OptSetup.hand_force = [0;0;0];
    
    % Should angular velocities be zero at the start and end of movement?
    % (true/false)
    OptSetup.start_at_rest = true;
    OptSetup.end_at_rest = true;
    
    % Set optimisation options
    OptSetup.N = nodes;
	OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
    OptSetup.OptimTol = 1e-3;
    OptSetup.FeasTol = 1e-3;
    OptSetup.initialguess = 'random';  % initial guess (see options in das3elbow_optimise_reach.m) 
    OptSetup.Wdata = 1;        % weight for the kinematic term in the cost function
    OptSetup.Weffort = 0.2;    % weight for the energy consumption term in the cost function
    OptSetup.equality_constraints = 1;

    % Create folder for results, if it does not already exist
    folder_name = 'reach';
    if ~exist(folder_name,'dir')
        mkdir(folder_name);
    end
    filename = [folder_name '/A'];
    tic
    das3elbow_optimize_reach([filename '_' num2str(nodes)],OptSetup);
    disp(['Time elapsted: ' num2str(toc) ' seconds'])
	% now do a series of optimizations, each time increasing number of nodes by a certain factor
	while (nodes < maxnodes)
        prev_nodes = nodes;
		% the following multiplies the number of nodes by the factor
		nodes = prev_nodes*factor;
		% don't exceed max nodes, and do the last optimization with tighter tolerances
		if (nodes >= maxnodes) 
			nodes = maxnodes;
			OptSetup.OptimTol = 1e-3;
            OptSetup.FeasTol = 1e-3;
		end
		% redo the optimization with previous result as initial guess
        OptSetup.N = nodes;
        OptSetup.initialguess = [filename '_' num2str(prev_nodes)];
        tic
        das3elbow_optimize_reach([filename '_' num2str(nodes)],OptSetup);
        disp(['Time elapsted: ' num2str(toc) ' seconds'])

	end

end