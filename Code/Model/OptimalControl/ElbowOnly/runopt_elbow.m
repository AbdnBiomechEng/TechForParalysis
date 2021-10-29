function runopt_elbow

	% run a sequence of optimizations, with mesh refinement

	maxnodes = 15;		% end close to this number of nodes
    nodes = 5;          % start with this number of nodes
    OptSetup.N = nodes;
	OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
    OptSetup.OptimTol = 1e-3;
    OptSetup.FeasTol = 1e-1;
    OptSetup.initialguess = 'random';  % initial guess (see options in das3_optimise_elbow.m) 
    OptSetup.Wdata = 100;    % weight for the kinematic term in the cost function
    OptSetup.Weffort = 1;    % weight for the energy consumption term in the cost function
    OptSetup.equality_constraints = 1;
    OptSetup.solver = 'IPOPT';

    factor = 2;

    x0 = load('equilibrium.mat');
    
    % Dynamic: elbow flexion-extension from 5 to 90 degrees
    data = [x0.x(1:14) x0.x(1:14)]'; % two time points
    data(1,13) = 5*pi/180;  % flexion in first time point is set at 5 degrees
    data(2,13) = 90*pi/180; % flexion in second time point is set at 90 degrees
    t = [0;1];  % the time vector

    % Static posture - just one node
    % nodes = 1;
    % maxnodes = 1;
    % data = x0.x(1:14)';
    % data(13) = 0.64;
    % t = 0;

    % Create folder for results, if it does not already exist
    folder_name = 'flexion90';
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