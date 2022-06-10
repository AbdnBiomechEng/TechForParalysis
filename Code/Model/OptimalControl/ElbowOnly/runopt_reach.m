function runopt_reach    
	% Run a sequence of optimizations, with mesh refinement
    
    factor = 1.6;    
    maxnodes = 40;		   % end close to this number of nodes
    nodes = 10;            % start with this number of nodes

    
    % Which model to use
    OptSetup.model_file = 'model_struct_A_FES.mat'; 
    
    
    % Set up FES stimulation
    OptSetup.stim_mus_groups = [1:3, 4*ones(1,9), 5:24, 4*ones(1,5), 25:29];
    
   
    % Input data
    OptSetup.t = [0;1];  % the time vector
    OptSetup.data_init = zeros(2,8);

    % Thoracohumeral angles from 5 to 60 degrees of flexion
    OptSetup.data_init(1,4:6) = [90,5,0]*pi/180;
    OptSetup.data_init(2,4:6) = [90,60,0]*pi/180;

    % Elbow flexion-extension from 60 to 0 degrees 
    OptSetup.data_init(1,7) = 60*pi/180;  
    OptSetup.data_init(2,7) = 0*pi/180; 

    % Pronation-supination from 120 to 90 degrees
    OptSetup.data_init(1,8) = 120*pi/180;  
    OptSetup.data_init(2,8) = 90*pi/180;  

    
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
%    OptSetup.initialguess = 'reach10/A_paralysed_50';  % initial guess (see options in das3elbow_optimise_reach.m) 
    OptSetup.Wdata = 10;        % weight for the kinematic term in the cost function
    OptSetup.Weffort = 0.2;     % weight for the energy consumption term in the cost function
    OptSetup.equality_constraints = 1;

    
    % Create folder for results, if it does not already exist
    folder_name = 'A_reach_controller';
    if ~exist(folder_name,'dir')
        mkdir(folder_name);
    end
    filename = [folder_name '/flexion_test'];
    
    for itest=1:20
        tic
        das3elbow_optimize_reach([filename '_' num2str(itest)],OptSetup);
        disp(['Time elapsted: ' num2str(toc) ' seconds'])
    end

    return;
    
    % Do a series of optimizations, each time increasing number of nodes by a certain factor
    tic
    das3elbow_optimize_reach([filename '_' num2str(nodes)],OptSetup);
    disp(['Time elapsted: ' num2str(toc) ' seconds'])
    
	while (nodes < maxnodes)
        prev_nodes = nodes;
		% the following multiplies the number of nodes by the factor
		nodes = round(prev_nodes*factor);
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