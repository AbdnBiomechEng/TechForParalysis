function CMC = setup_ComputedMuscleControl(varargin)

% function CMC = setup_ComputedMuscleControl(data, varagin)
%
% This function will create a Setup_CMC.XML file for use with OpenSIM
% based on the data from the C3D file and the names of the model files
% which will be used to scale the model.
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file)      
%         'MOTFile' - desired kinematics file MOT file for ID
%         'GRFFile' - filename string of the XML file containing GRF
%               information 
%
%         OPTIONAL RECOMMENDED parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'CMCTaskFile' - filename string of the Tasks XML file
%         'CMCForceFile' - filename string of the Actuator XML file
%         'CMCConstraintsFile' - File containing the constraints on the
%               controls.
%         'RRAControlsFile' - File containing the controls output by RRA. 
%               These can be used to place constraints on the residuals during CMC.
%
%         OPTIONAL parameters
%         'OutputPrecision' - number between 1-50 (default 8)
%         'LowPassFilterFreq' - number between 1 and 60 for low pass filter 
%                       of kinematics (default -1 = none)
%         'DirectoryName' - string which is the directory name to be made
%                           for the results (default 'CMCResults'
%         'ReplaceForceSet' - 'true' or 'false' to deteremine whether 
%                          the model actuator set is replaced with those
%                          from actuator file (default 'true' -replaced)
%         'InitialTime' - initial time for ID to run from (defaults to
%                       start of the MOT file time)
%         'FinalTime' - final time for ID to run until (defaults to
%                       time at end of the MOT file)
%         'AuxillaryStates' - whether or not to compute equilibrium values for
%       		    states other than the coordinates or speeds (default
%       		    - 'false')
%         'MaxIntegratorSteps' - maximum number of integrator steps
%                   (default=20000)
%         'MaxStepSize' - maximum integrator step size (default=1)
%         'MinStepSize' - maximum integrator step size (default=1e-6)
%         'ErrorTol' - inegrator error tolerance (default=1e-5)
%         'OptimizerAlgorithm' - 'ipopt' or 'cfsqp' (default - 'ipopt')
%         'OptimizerDx' - optimizer derivative dx value (default = 1e-4)
%         'OptimConvergCrit' - optmiser convergence criterion value
%               (default = 1e-6)
%         'OptimMaxIterations' - Maximum number of iterations for the optimize
%               (default - 30000)
%         'CMCTimeWindow' - Time window over which the desired actuator
%               forces are achieved (default - 0.01);
%         'UseCurvatureFilter' - Flag (true or false) indicating whether or 
%               not to use the curvature filter. Setting this flag to true
%               can reduce oscillations in the computed muscle excitations
%               (default - 'false')	    
%         'FastTarget' - Flag (true or false) indicating whether to    
%		    use the fast CMC optimization target. The fast target requires 
%		    the desired accelerations to be met (default - 'false')
%
% Output - CMC (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  CMC = setup_CMC('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,'GRFFile',GRFFile,...
%        'CMCTaskFile',CMCTaskFile,'CMCActuatorFile',CMCActuatorFile,'CMCControlFile',CMCControlFile,...
%        'AdjCOMRes','true','OptimMaxIter',20000,'LowPassFilterFreqLoad',-1,...
%        'LowPassFilterFreq',-1);
% 
% Written by Glen Lichtwark (The University of Queensland)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy
ModelFile = [];

ReplaceForceSet = 'false';
CMCForceFile = [];

DirectoryName = 'CMCResults';
OutputPrecision = 12;
InitialTime = [];
FinalTime = [];

AuxillaryStates = 'true';
MaxIntegratorSteps = 20000;
MaxStepSize = 1;
MinStepSize = 1e-6;
ErrorTol = 1e-5;

GRFFile = [];
MOTFile = [];

CMCTaskFile = [];
CMCConstraintsFile = [];
RRAControlsFile = ''; 

LowPassFilterFreq = -1;

CMCTimeWindow = 0.01;
UseCurvatureFilter = 'true';
FastTarget = 'true';

OptimizerAlgorithm = 'ipopt';
OptimizerDx = 1e-4;
OptimConvergCrit = 1e-6;
OptimMaxIterations = 30000;

data = [];

if ~isempty(varargin)
    if rem(length(varargin),2)
        error('Incorrect input arguments - must specify property and input')
    end
    for i = 1:2:length(varargin)
       n = varargin{i+1};
       eval([varargin{i} '= n;']); 
    end    
end

% define the initial and final times for inverse dynamics from the data structure 
% if this is passed to function and these aren't prescribed otherwise
if ~isempty(data)
    if isfield(data,'time')
        if isempty(InitialTime)
            InitialTime = data.time(1);
        end
        if isempty(FinalTime)
            FinalTime = data.time(end);
        end
    end
end

if isempty(MOTFile)
    error('MOTFile must be included for RRA analysis')
end

%setup root and RRATool
root = 'OpenSimDocument';

V.ATTRIBUTE.Version = '20302';
% define some names for outputting... use the data in the data structure to
% limit the filename size to important parts if data structure is passed
if ~isempty(data)
    Model = data.Name;
    [~,trial, ~] = fileparts(data.C3D_Filename);
else [~, Model, ~] = fileparts(ModelFile);
    [~, trial, ~] = fileparts(MOTFile);
end

V.CMCTool.ATTRIBUTE.name = [Model '_' trial '_CMC'];

% define the model file
if ~isempty(ModelFile)
    V.CMCTool.model_file = ModelFile;
else error('Please specify a model file')
end

% define the force set and determine whether this is to replace or
% append to current actuator set
V.CMCTool.replace_force_set = ReplaceForceSet;
V.CMCTool.force_set_files = CMCForceFile;

% define results directory and precision
V.CMCTool.results_directory = ['./' DirectoryName];
V.CMCTool.output_precision = OutputPrecision;

% Define times to perform analysis over 
V.CMCTool.initial_time = num2str(InitialTime,6);
V.CMCTool.final_time = num2str(FinalTime,6);

V.CMCTool.solve_for_equilibrium_for_auxiliary_states = AuxillaryStates;

% define integrator settings
V.CMCTool.maximum_number_of_integrator_steps = MaxIntegratorSteps;
V.CMCTool.maximum_integrator_step_size = MaxStepSize;
V.CMCTool.minimum_integrator_step_size = MinStepSize;
V.CMCTool.integrator_error_tolerance = ErrorTol;

% define the input files
V.CMCTool.external_loads_file = GRFFile;
V.CMCTool.desired_kinematics_file = MOTFile;

% define input and output files for CMC
V.CMCTool.task_set_file = CMCTaskFile;
V.CMCTool.constraints_file = CMCConstraintsFile;
V.CMCTool.rra_controls_file = RRAControlsFile;

% define the filtering low pass frequency
V.CMCTool.lowpass_cutoff_frequency = LowPassFilterFreq;

% CMC properties
V.CMCTool.cmc_time_window = CMCTimeWindow;
V.CMCTool.use_curvature_filter = UseCurvatureFilter;
V.CMCTool.use_fast_optimization_target = FastTarget;

% optimiser settings
V.CMCTool.optimizer_algorithm = OptimizerAlgorithm;
V.CMCTool.optimizer_derivative_dx = OptimizerDx;
V.CMCTool.optimizer_convergence_criterion = OptimConvergCrit;
V.CMCTool.optimizer_max_iterations = OptimMaxIterations;

V.CMCTool.optimizer_print_level = 0;
V.CMCTool.use_verbose_printing = 'false';

fileout = [Model '_Setup_ComputedMuscleControl.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

CMC = V;
        