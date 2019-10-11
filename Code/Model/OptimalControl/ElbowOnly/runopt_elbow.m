function runopt_elbow

	% run a sequence of optimizations, with mesh refinement

	maxnodes = 17;		% end close to this number of nodes
    nodes = 17;
    OptSetup.N = nodes;
	OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
    OptSetup.OptimTol = 1e-3;
    OptSetup.FeasTol = 1e-1;
    OptSetup.initialguess = 'eqLce';
    OptSetup.Wdata = 1;
    OptSetup.Weffort = 0;
    OptSetup.equality_constraints = 1;
    OptSetup.solver = 'IPOPT';

    factor = 2;

    % data: static position, with elbow going from zero to 90 degrees
    x0 = load('equilibrium.mat');
    data = [x0.x(1:14) x0.x(1:14) x0.x(1:14)]'; % three time points
    data(1,13) = 0*pi/180; % change flexion in first time point to 0 degrees
    data(2,13) = 90*pi/180; % change flexion in second time point to 90 degrees
    data(3,13) = 0*pi/180; % change flexion in third time point to 0 degrees
    t = [0;0.5;1];

    % Create folder for results, if it does not already exist
    folder_name = 'flexion-no-effort';
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