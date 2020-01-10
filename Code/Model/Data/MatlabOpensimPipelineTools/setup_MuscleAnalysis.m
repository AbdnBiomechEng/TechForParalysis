function muscle_analysis = setup_MuscleAnalysis(varargin)

% function muscle_analysis = setup_MuscleAnalysis(varargin)
%
% This function will create a Setup_Muscle_Analysis.XML file for use with OpenSIM 
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file)    
%
%         OPTIONAL INPUT FILES to drive analysis
%         'ExternalLoadsFile' - filename string of the XML file containing force data -
%           pointer to STO file --> all external forces can be added in
%           this file (not just GRF forces).
%         'StatesFile' - model states file (STO) from CMC analysis (e.g.
%           ResultsCMC/subject1_walk1_states.sto') NOTE: IF THE STATES FILE
%           IS USED DO NOT USE THE COORDINATES FILE!        
%         'MOTFile' - coordinates file (STO) from CMC analysis (e.g.
%           ResultsCMC/subject1_walk1_Kinematics q.sto')
%         'ForceSetFile' - File containing forces acting on model not
%               contained in the states file
%         'SpeedsFile' - Storage file (.sto) containing the time history of the generalized
%		    speeds for the model. If coordinates_file is used in place of
%		    states_file, these can be optionally set as well to give the speeds.
%		    If not specified, speeds will be computed from coordinates by
%		    differentiation.
%
%         OPTIONAL MUSLCE ANALYSIS PARAMETERS
%         'AnalysisOn' - whether to use the muscle analysis (default -
%                   true')
%         'StartTime' - defaults to Analyze initial time
%         'EndTime - defaults to Analyze final time
%         'StepInterval' - how often (in integrator steps) to store results during simulation
%         'InDegrees' - whether results are printed in degrees or radians
%               (default - true)
%         'MuscleList' - String array containing names of muscles for which to perform the analysis
%               (default - 'all', other e.g. 'gaslat_r gaslat_l soleus_r soleus_l')
%         'MomentArmCoordinates' - Cell array of generalized coordinates
%               for which to calculate moment arms (default - empty, possible 'all')
%         'ComputeMoments' - whether to compute coordinate moments (default
%               = 'true')
%
%         OPTIONAL CONTROLLER parameters
%         'ActuatorList' -  list of actuators that this controller will control.The keyword ALL
%				indicates the controller will controll all the acuators in the model
%         'EnableController' - whether or not controller is on (defaults to
%               'true' if controls file is provided, 'false' if not)
%         'ControlsFile' - actuator control file (XML) 
%                --> drives the muscles and other actuators
%
%         OPTIONAL ANALYZE parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'ReplaceForceSet' - 'true' or 'false' to deteremine whether 
%              the model actuator set is replaced with those
%              from actuator file (default 'false')
%         'SolveForEquilibrium' -  'true' or 'false' (default - 'true')
%         'InitialTime' - initial time for ID to run from (defaults to
%              start of the MOT file time)
%         'FinalTime' - final time for ID to run until (defaults to
%              time at end of the MOT file)
%         'DirectoryName' - string which is the directory name to be made
%              for the results (default 'PI_Results')
%         'OutputPrecision' - number between 1-50 (default 8)
%         'MaxNumSteps' - maximum number of integrator steps (default -
%               20000)
%         'MaxStepSize' - maximum integrator step size (default=1)
%         'MinStepSize' - maximum integrator step size (default=1e-5)
%         'ErrorTol' - inegrator error tolerance (default=1e-5)
%         'LowPassFilterForKinematics' - -Low-pass cut-off frequency for filtering
%               the coordinates_file data
%
% Output - MuscleAnalysis (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  MuscleAnalysis = setup_muscle_analysis('ModelFile','c:\data\model123.osim', ...
%                 'ForceSetFile',CMC_actuator.xml, 'MOTFile', walk1_ik.mot, ...
%                 'ExternalLoadsFile', walk1_grf.mot, 'StatesFile', ...
%                 '/CMC/states.sto', 'ControlsFile','/CMC/controls.xml', ...
%                 'MuscleList','all','ComputeMoments','true');
%                    
%
% Written by Glen Lichtwark (The University of Queensland)

% setup default input files to emtpy

%Model specification
ModelFile = [];
ReplaceForceSet = 'false';
ForceSetFile = [];

% output
DirectoryName = 'MuscleAnalysis';
OutputPrecision = 8;

%Time interval
InitialTime = [];
FinalTime = [];

%solve for equilibrium for auxiliary states
SolveForEquilibrium = 'true';

%Integrator Settings
MaxNumSteps = 20000;
MaxStepSize = 1;
MinStepSize = 1e-8;
ErrorTol = 1e-5;

%MuscleAnalysis settings
AnalysisOn = 'true';
StartTime = [];
EndTime = [];
StepInterval = 1;
InDegrees = 'true';
MuscleList = 'all';
MomentArmCoordinates = ' ';
ComputeMoments = 'true';

%ControllerSettings
ActuatorList = 'all';
EnableController = [];
ControlsFile = [];

% storage files containing motion, states and forces
ExternalLoadsFile = [];
StatesFile = [];
MOTFile = []; %MOTFile
SpeedsFile = []; 

% filter settings for coordinates file
LowPassFilterForKinematics = -1;

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

% define the initial and final times for muscle analysis from the data structure 
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

%setup root and AnalyzeTool
root = 'OpenSimDocument';

V.ATTRIBUTE.Version = '20302';
[~, name, ~] = fileparts(ModelFile);
V.AnalyzeTool.ATTRIBUTE.name = name;

if ~isempty(ModelFile)
    V.AnalyzeTool.model_file = ModelFile;
else error('Please specify a model file')
end
V.AnalyzeTool.replace_force_set = ReplaceForceSet;
V.AnalyzeTool.force_set_files = ForceSetFile;

%output
V.AnalyzeTool.results_directory = DirectoryName;
V.AnalyzeTool.output_precision = OutputPrecision;

%Time interval
V.AnalyzeTool.initial_time = num2str(InitialTime,6);
V.AnalyzeTool.final_time = num2str(FinalTime,6);

% solve for equilibrium for auxiliary states
V.AnalyzeTool.solve_for_equilibrium_for_auxiliary_states = SolveForEquilibrium;

%Integrator Settings
V.AnalyzeTool.maximum_number_of_integrator_steps = MaxNumSteps;
V.AnalyzeTool.maximum_integrator_step_size = MaxStepSize;
V.AnalyzeTool.minimum_integrator_step_size = MinStepSize;
V.AnalyzeTool.integrator_error_tolerance = ErrorTol;

%Analysis set
V.AnalyzeTool.AnalysisSet.ATTRIBUTE.name = 'Analyses';
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.ATTRIBUTE.name = 'MuscleAnalysis';
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.on = AnalysisOn;
if isempty(StartTime)
    V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.start_time = InitialTime;
else V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.start_time = StartTime;
end
if isempty(EndTime)
    V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.end_time = FinalTime;
else V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.end_time = EndTime;
end
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.step_interval = StepInterval;
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.in_degrees = InDegrees;
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.muscle_list = MuscleList;
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.moment_arm_coordinate_list = MomentArmCoordinates;
V.AnalyzeTool.AnalysisSet.objects.MuscleAnalysis.compute_moments = ComputeMoments;
   
%controller set
V.AnalyzeTool.ControllerSet.ATTRIBUTE.name = 'Controllers';
% V.AnalyzeTool.ControllerSet.objects.ControlSetController.actuator_list = ActuatorList;
V.AnalyzeTool.ControllerSet.objects = [];
V.AnalyzeTool.ControllerSet.groups = [];
if ~isempty(EnableController)
    V.AnalyzeTool.ControllerSet.objects.ControlSetController.enable_controller = EnableController;
    if ~isempty(ControlsFile)
        V.AnalyzeTool.ControllerSet.objects.ControlSetController.controls_file = ControlsFile;
    end
end


% storage files containing motion, states and forces
V.AnalyzeTool.external_loads_file = ExternalLoadsFile;
V.AnalyzeTool.states_file = StatesFile;
V.AnalyzeTool.coordinates_file = MOTFile;
V.AnalyzeTool.speeds_file = SpeedsFile;

V.AnalyzeTool.lowpass_cutoff_frequency_for_coodinates = LowPassFilterForKinematics;

fileout = [data.Name '_Setup_MuscleAnalysis.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

muscle_analysis = V;
        