function runopt
% Run a sequence of optimizations, with mesh refinement

% Add seed for reproducibility
rng(123);

%% Variable names
% DOFs in the model
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};

% Instead of clavicular and scapular angles, we may have
% thoraco-humeral angles angles instead (YZY sequence)
thor_hum = {'Trunk_Hum_y','Trunk_Hum_z','Trunk_Hum_yy'};

% There may also be external forces at the hand
% defined in the global frame: +X is laterally, +Y is upwards and +Z is posteriorly
hand_forces = {'hand_force_x','hand_force_y','hand_force_z'};

%% Read in data and set optimisation parameters
% It should contain "time" and some of the variables above
data = readtable('input_data.csv');

% Choose which columns to use (time, locked and tracking angles, and hand
% forces if applicable)
% or leave empty to use all columns in the input data
input_variables = {};

% Choose which degrees of freedom to lock
lockeddofs = {'TH_x','TH_y','TH_z','SC_y'};

% Choose which model parameter file to use
OptSetup.model = 'model_struct.mat';

% Should angular velocities be zero at the start and end of movement?
% (true/false)
OptSetup.start_at_rest = true;
OptSetup.end_at_rest = true;

% Set number of nodes
maxnodes = 12;		  % end close to this number of nodes
nodes = 3;            % start with this number of nodes
factor = 1.5;         % increase number of nodes by this factor

% Set optimisation options
OptSetup.N = nodes;
OptSetup.MaxIter = 10000;	% max number of iterations for each optimization
OptSetup.OptimTol = 1e-3;
OptSetup.FeasTol = 1e-3;
OptSetup.initialguess = 'random';  % initial guess (see options in das3_optimise.m)

% Cost function
OptSetup.Wdata = 1;        % weight for the kinematic term in the cost function
OptSetup.Weffort = 0.05;    % weight for the energy consumption term in the cost function
OptSetup.Wscap = 0.1;      % weight for scapulo-thoracic gliding plane (if missing, assumed to be constraint)
OptSetup.Whum = 0.1;       % weight for glenohumeral stability (if missing, assumed to be constraint)

% Create folder for results, if it does not already exist
folder_name = 'full_arm';
% Filename
filename = [folder_name '/output_scaphum'];

%% Check input
if isempty(input_variables)
    input_variables = data.Properties.VariableNames;
end

% Check column names
if ~ismember('time',input_variables)
    error('One of the columns in the input data should be "time"');
end
all_column_names = ['time',dofnames,thor_hum,hand_forces];
if ~all(ismember(input_variables,all_column_names))
    error(['Error reading in data: ' input_variables{~ismember(input_variables,all_column_names)} ' is not a valid column name']);
end
OptSetup.input_variables = input_variables;
data = data(:,input_variables);

% Check locked dofs
if ~all(ismember(lockeddofs,dofnames))
    error(['Error settting locked DOFs: ' lockeddofs{~ismember(lockeddofs,[dofnames,thor_hum])} ' is not a valid DOF name']);
end

OptSetup.lockeddof_indata = get_indeces(lockeddofs,input_variables);
OptSetup.lockeddof_inx = get_indeces(lockeddofs,dofnames);

% Find tracking dofs
tracking_dofs = input_variables(ismember(input_variables,dofnames));
tracking_dofs = tracking_dofs(~ismember(tracking_dofs,lockeddofs));
OptSetup.tracking_indata = get_indeces(tracking_dofs,input_variables);
OptSetup.tracking_inx = get_indeces(tracking_dofs,dofnames);

% Do we have hand forces?
handF = input_variables(ismember(input_variables,hand_forces));
OptSetup.handF_indata = get_indeces(handF,input_variables);
OptSetup.handF_inF = get_indeces(handF,hand_forces);

% Do we have thoraco-humeral angles?
thorhum = input_variables(ismember(input_variables,thor_hum));
OptSetup.thorhum_indata = get_indeces(thorhum,input_variables);
OptSetup.thorhum_inx = get_indeces(thorhum,thor_hum);

% Check output folder and create if needed
if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

%% Run optimisation
tic
das3_optimize(data,[filename '_' num2str(nodes)],OptSetup);
disp(['Time elapsted: ' num2str(toc) ' seconds'])

% now do a series of optimizations, each time increasing number of nodes by a certain factor
prev_nodes = nodes;
nodes = ceil(prev_nodes*factor);
while (nodes <= maxnodes)
    % redo the optimization with previous result as initial guess
    OptSetup.N = nodes;
    OptSetup.initialguess = [filename '_' num2str(prev_nodes)];
    tic
    das3_optimize(data,[filename '_' num2str(nodes)],OptSetup);
    disp(['Time elapsted: ' num2str(toc) ' seconds'])
    prev_nodes = nodes;
    nodes = ceil(prev_nodes*factor);
end


    function i_invector = get_indeces(variables,invector)
        % function that calculates the indeces of a cell array in a larger
        % array
        nvar = length(variables);
        i_invector = zeros(nvar,1);
        for ivar=1:nvar
            i_invector(ivar) = find(ismember(invector,variables{ivar}));
        end
    end
end