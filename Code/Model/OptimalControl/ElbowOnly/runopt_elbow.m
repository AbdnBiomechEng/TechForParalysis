function runopt_elbow    
	% run a sequence of optimizations, with mesh refinement
    
    % add seed for reproducibility
    rng(123);

    factor = 2;
    x0 = load('equilibrium.mat');
    
    % Dynamic movement 

    maxnodes = 10;		 % end close to this number of nodes
    nodes = 10;          % start with this number of nodes

    % Elbow flexion-extension from 5 to 90 degrees 
    data = [x0.x(1:14) x0.x(1:14)]'; % two time points
    data(1,13) = 5*pi/180;  % flexion in first time point is set at 5 degrees
    data(2,13) = 90*pi/180; % flexion in second time point is set at 90 degrees
    t = [0;1];  % the time vector

    % Static posture - just one node
%     nodes = 1;
%     maxnodes = 1;
%     data = x0.x(1:14)';
%     data(13) = 90*pi/180;
%     t = 0;

    % Add external force at the centre of mass of the hand 
    % (defined in the global frame: +X is laterally, +Y is upwards and +Z is posteriorly)
    OptSetup.hand_force = [0;0;0];
    
    % Which degrees of freedom to lock (if none: empty vector [])
    OptSetup.lockeddofs = 1:12;
    
    % Should angular velocities be zero at the start and end of movement?
    % (true/false)
    OptSetup.start_at_rest = true;
    OptSetup.end_at_rest = true;
    
    % Set optimisation options
    OptSetup.N = nodes;
	OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
    OptSetup.OptimTol = 1e-3;
    OptSetup.FeasTol = 1e-1;
    OptSetup.initialguess = 'random';  % initial guess (see options in das3_optimise_elbow.m) 
    OptSetup.Wdata = 1;        % weight for the kinematic term in the cost function
    OptSetup.Weffort = 0.1;    % weight for the energy consumption term in the cost function
    OptSetup.equality_constraints = 1;
    OptSetup.solver = 'IPOPT'; % IPOPT or fmincon (fmincon not added yet)

    % Create folder for results, if it does not already exist
    folder_name = 'flexion5_90';
    if ~exist(folder_name,'dir')
        mkdir(folder_name);
    end
    filename = [folder_name '/output'];
    
    das3_optimize_elbow(data,t,[filename '_' num2str(nodes)],OptSetup);
	% now do a series of optimizations, each time increasing number of nodes by a certain factor
	while (nodes < maxnodes)
        prev_nodes = nodes;
		% the following multiplies the number of time intervals (nodes-1) by the factor
		nodes = max(nodes+1, round((nodes-1)*factor+1));		% at least increase number of nodes by 1!
		% don't exceed max nodes, and do the last optimization with tighter tolerances
		if (nodes >= maxnodes) 
			nodes = maxnodes;
			OptSetup.OptimTol = 1e-3;
            OptSetup.FeasTol = 1e-1;
		end
		% redo the optimization with previous result as initial guess
        OptSetup.N = nodes;
        OptSetup.initialguess = [filename '_' num2str(prev_nodes)];
        das3_optimize_elbow(data,t,[filename '_' num2str(nodes)],OptSetup);
	end

end