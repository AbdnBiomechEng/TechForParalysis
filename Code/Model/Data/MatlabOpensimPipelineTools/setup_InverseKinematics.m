function IK = setup_InverseKinematics(varargin)

% function IK = setup_InverseKinematics(data)
%
% This function will create a Setup_IK.XML file for use with OpenSIM (ver
% 2.31 and beyond) based on the data from the C3D file and the names of the 
% model files which will be used to track the kinematics using the scaled
% model.
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file - scaled)
%         'MarkerFile' - filename string of the TRC file containing
%           coordinates of the markers to be used in scaling (defaults to
%           the data.TRC_Filename if 'data' structure passed)
%         'IKTasksFile' - filename string of the ScaleTasksFile
%
%          OPTIONAL
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m
%         'InitialTime' - initial time for IK to run from (defaults to
%                   start of the MOT file time)
%         'FinalTime' - final time for IK to run until (defaults to
%                   time at end of the MOT file)
%         'ResultsDirectory' - name of directory where IK results are
%                   written (default - ./)
%         'InputDirectory' - name of directory where data arrises from
%                   (default - empty --> current directory used)
%         'ConstraintWeight' - A positive scalar that is used to weight the 
%                   importance of satisfying constraints.A weighting of 
%                   'Infinity' or if it is unassigned results in the
%                   constraints being strictly enforced (default -
%                   'infinity')
%         'Accuracy' - The accuracy of the solution in absolute terms. I.e. 
%                   the number of significant digits to which the solution 
%                   can be trusted (default = 0.00005)
%         'ReportErrors' - Flag (true or false) indicating whether or not 
%                   to report marker and coordinate errors from the inverse 
%                   kinematics solution (default - 'true')
%         'OutputFile' - Name of the motion file (.mot) to which the
%                   results should be written (default - empty [] -->
%                   the same name as TRC file will be used')
%
% Output - IK (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  V = setup_InverseKinematics('data', data,'ModelSetFile','c:\models\model123.osim' ...
%              'IKTasksFile', 'c:\models\model123_IKTaskFile.xml')
% 
% Written by Glen Lichtwark (The University of Queensland)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy

ModelFile = [];
MarkerFile = [];
IKTasksFile = [];
InitialTime = [];
FinalTime = [];

CoordinateFile = 'Unassigned';

ResultsDirectory = './';
InputDirectory = [];
ConstraintWeight = 'infinity';
Accuracy = 0.00001;

ReportErrors = 'true';

OutputFile = [];

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

% define the initial and final times for the scaling if these aren't
% prescribe and the MarkerFile if this isn't prescribed
if ~isempty(data)
    if isfield(data,'time')
        if isempty(InitialTime)
            InitialTime = data.time(1);
        end
        if isempty(FinalTime)
            FinalTime = data.time(end);
        end
        if isempty(MarkerFile)
            MarkerFile = data.TRC_Filename;
        end
        if isempty(OutputFile)
            OutputFile = [data.TRC_Filename(1:end-4) '_ik.mot'];
        end
    end
end

[~, Name, ~] = fileparts(ModelFile);

%setup root and InverseKinematicsTool
root = 'OpenSimDocument';
        
V.ATTRIBUTE.Version = '20302';
V.InverseKinematicsTool.ATTRIBUTE.name = Name;

% define results and input directories
V.InverseKinematicsTool.results_directory = ResultsDirectory;
if isempty(InputDirectory)
    V.InverseKinematicsTool.input_directory = cd;
else V.InverseKinematicsTool.input_directory = InputDirectory;
end

% define model file
if ~isempty(ModelFile)
    V.InverseKinematicsTool.model_file = ModelFile;
end

% define solver parameters
V.InverseKinematicsTool.constraint_weight = ConstraintWeight;
V.InverseKinematicsTool.accuracy = Accuracy;

% Define IKTaskSet file
if ~isempty(IKTasksFile)
    V.InverseKinematicsTool.IKTaskSet.ATTRIBUTE.file = IKTasksFile;
end

% Define marker (TRC) file and coordinate file (if present)
V.InverseKinematicsTool.marker_file = MarkerFile;

V.InverseKinematicsTool.coordinate_file = CoordinateFile;

% Define the time range
V.InverseKinematicsTool.time_range = num2str([InitialTime FinalTime]);

% Define outputs - report errors or output filename
V.InverseKinematicsTool.report_errors = ReportErrors;
if isempty(OutputFile)
     V.InverseKinematicsTool.output_motion_file = [MarkerFile(1:end-4) '_ik.mot'];
else V.InverseKinematicsTool.output_motion_file = OutputFile;
end

fileout = [data.Name '_Setup_InverseKinematics.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

IK = V;
        