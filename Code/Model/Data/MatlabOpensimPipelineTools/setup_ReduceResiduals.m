function RRA = setup_ReduceResiduals(varargin)

% function RRA = setup_RRA(data)
%
% This function will create a Setup_RRA.XML file for use with OpenSIM
% based on the data from the C3D file and the names of the model files
% which will be used to scale the model.
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file)      
%         'MOTFile' - desired kinematics file MOT file for ID
%         'GRFFile' - filename string of the XML file containing GRF
%         information 
%
%         OPTIONAL RECOMMENDED parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'RRATaskFile' - filename string of the Tasks XML file
%         'RRAForceFile' - filename string of the Actuator XML file
%         'RRAConstraintsFile' - File containing the constraints on the
%               controls.
%         'RRAControlsFile' - File containing the controls output by RRA. 
%               These can be used to place constraints on the residuals during CMC.
%
%         OPTIONAL parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'OutputPrecision' - number between 1-50 (default 8)
%         'LowPassFilterFreq' - number between 1 and 60 for low pass filter 
%                       of kinematics (default -1 = none)
%         'DirectoryName' - string which is the directory name to be made
%                           for the results (default 'RRAResults'
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
%         'MinStepSize' - maximum integrator step size (default=1e-8)
%         'ErrorTol' - inegrator error tolerance (default=1e-5)
%         'FineErrorTol' - fine inegrator error tolerance (default=0.00001)
%         'OptimizerAlgorithm' - 'ipopt' or 'cfsqp' (default - 'ipopt')
%         'OptimizerDx' - optimizer derivative dx value (default = 0.0001)
%         'OptimConvergCrit' - optmiser convergence criterion value
%                       (default = 1e-6)
%         'AdjCOMRes' - adjust the COM to reduce residuals (default - 'true')
%         'AdjustedCOMBody' - the body to reference the adjusted COM to if
%                       for RRA1 (default - 'torso' -->
%                       heaviest body)
%         'InitialTimeCOMAdjustment' - time to start the analysis for
%               performing COM adjustment (default = -1, uses Initial 
%                times for simulation)'
%         'FinalTimeCOMAdjustment' - time to finish the analysis for
%               performing COM adjustment (default = -1, uses Final 
%                times for simulation)'
%         'OutputModelFile' - Name of the model to output (defaults to
%               subject name + '_SCALED.osim' if data structure sent with C3D info) 
%
%
% Output - RRA (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  RRA = setup_RRA('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,'GRFFile',GRFFile,...
%        'RRATaskFile',RRATaskFile,'RRAActuatorFile',RRAActuatorFile,'RRAControlFile',RRAControlFile,...
%        'AdjCOMRes','true','OptimMaxIter',20000,'LowPassFilterFreq',-1);
% 
% Written by Glen Lichtwark (Griffith University)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy
ModelFile = [];

ReplaceForceSet = 'true';
RRAForceFile = [];

DirectoryName = 'RRAResults';
OutputPrecision = 12;
InitialTime = [];
FinalTime = [];

AuxillaryStates = 'false';
MaxIntegratorSteps = 20000;
MaxStepSize = 1;
MinStepSize = 1e-8;
ErrorTol = 1e-5;
FineErrorTol = 0.00001;

GRFFile = [];
MOTFile = [];

RRATaskFile = [];
RRAConstraintsFile = [];
RRAControlsFile = ''; 

LowPassFilterFreq = -1;

OptimizerAlgorithm = 'ipopt';
OptimizerDx = 0.0001;
OptimConvergCrit = 1e-6;

AdjCOMRes = 'true';
AdjustedCOMBody = 'torso';
InitialTimeCOMAdjustment = -1; 
FinalTimeCOMAdjustment = -1;

OutputModelFile = [];

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
    [~,trial, ~] = fileparts(data.TRC_Filename);
else [~, Model, ~] = fileparts(ModelFile);
    [~, trial, ~] = fileparts(MOTFile);
end

V.RRATool.ATTRIBUTE.name = [Model '_' trial '_RRA'];

% define the model file
if ~isempty(ModelFile)
    V.RRATool.model_file = ModelFile;
else error('Please specify a model file')
end

% define the force set and determine whether this is to replace or
% append to current actuator set
V.RRATool.replace_force_set = ReplaceForceSet;
V.RRATool.force_set_files = RRAForceFile;

% define results directory and precision
V.RRATool.results_directory = DirectoryName;
V.RRATool.output_precision = OutputPrecision;

% Define times to perform analysis over 
V.RRATool.initial_time = num2str(InitialTime,6);
V.RRATool.final_time = num2str(FinalTime,6);

V.RRATool.solve_for_equilibrium_for_auxiliary_states = AuxillaryStates;

% define integrator settings
V.RRATool.maximum_number_of_integrator_steps = MaxIntegratorSteps;
V.RRATool.maximum_integrator_step_size = MaxStepSize;
V.RRATool.minimum_integrator_step_size = MinStepSize;
V.RRATool.integrator_error_tolerance = ErrorTol;
V.RRATool.integrator_fine_tolerance = FineErrorTol;

% define the input files
V.RRATool.external_loads_file = GRFFile;
V.RRATool.desired_kinematics_file = MOTFile;

% define input and output files for RRA
V.RRATool.task_set_file = RRATaskFile;
V.RRATool.constraints_file = RRAConstraintsFile;
V.RRATool.rra_controls_file = RRAControlsFile;

% define the filtering low pass frequency
V.RRATool.lowpass_cutoff_frequency = LowPassFilterFreq;

% optimiser settings
V.RRATool.optimizer_algorithm = OptimizerAlgorithm;
V.RRATool.optimizer_derivative_dx = OptimizerDx;
V.RRATool.optimizer_convergence_criterion = OptimConvergCrit;

V.RRATool.adjust_com_to_reduce_residuals = AdjCOMRes;
V.RRATool.adjusted_com_body = AdjustedCOMBody;

V.RRATool.initial_time_for_com_adjustment  = InitialTimeCOMAdjustment;
V.RRATool.final_time_for_com_adjustment  = FinalTimeCOMAdjustment;

if isempty(OutputModelFile)
    V.RRATool.output_model_file = [Model '_' trial '_RRA_adjusted.osim'];
else V.RRATool.output_model_file = OutputModelFile;
end

V.RRATool.use_verbose_printing = 'false';

fileout = [Model '_Setup_ReduceResiduals.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

RRA = V;
        