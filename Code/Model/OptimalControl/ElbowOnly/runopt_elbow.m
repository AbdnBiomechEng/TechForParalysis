function runopt_elbow

	% run a sequence of optimizations, with mesh refinement

	maxnodes = 40;		% end close to this number of nodes
    nodes = 8;
    OptSetup.N = nodes;
	OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
    OptSetup.OptimTol = 1e-1;
    OptSetup.FeasTol = 1e-1;
    OptSetup.initialguess = 'init';
    OptSetup.Wdata = 1000;
    OptSetup.Weffort = 1;
    OptSetup.equality_constraints = 1;
    OptSetup.solver = 'IPOPT';

    factor = 2;

    % data: static position, with elbow going from zero to 90 degrees
    x0 = load('equilibrium.mat');
    data = [x0.x(1:14) x0.x(1:14)]'; % two time points
    data(2,13) = 90*pi/180; % change flexion in second time point to 90 degrees
    t = [0;1];
            
    filename = 'opt_results/elbow_flexion';
    
    das3_optimize_elbow(data,t,[filename '_' num2str(nodes)],OptSetup);
	% now do a series of optimizations, each time increasing number of nodes by a certain factor
	while (nodes < maxnodes)
        prev_nodes = nodes;
		% the following multiplies the number of time intervals (nodes-1) by the factor
		nodes = max(nodes+1, round((nodes-1)*factor+1));		% at least increase number of nodes by 1!
		% don't exceed max nodes, and do the last optimization with tighter tolerances
		if (nodes >= maxnodes) 
			nodes = maxnodes;
			OptSetup.OptimTol = 1e-5;
            OptSetup.FeasTol = 1e-5;
		end
		% redo the optimization with previous result as initial guess
        OptSetup.N = nodes;
        OptSetup.initialguess = [filename '_' num2str(prev_nodes)];
        das3_optimize_elbow(data,t,[filename '_' num2str(nodes)],OptSetup);
	end

end