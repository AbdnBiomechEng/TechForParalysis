function ID = setup_InverseDynamics(varargin)

% function ID = setup_ID(data)
%
% This function will create a Setup_ID.XML file for use with OpenSIM
% based on the data from the C3D file and the names of the model files
% which will be used to scale the model.
% 
% Input - 'ModelFile' - string which is the filename (including path) of the
%                OSIM file (model file)
%         'MOTFile' - filename string of the MOT file for ID
%         'GRFFile' - filename string of the XML file containing GRF info
%                   and pointing to the MOT file containing the GRF data
%         (if this is not specified, the model will use actuator set in
%         model)
%
%         OPTIONAL parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m and any other processing (e.g. IK)
%         'ResultsDirectory' - name of directory where IK results are
%                   written (default - ./)
%         'InitialTime' - initial time for ID to run from (defaults to
%                       start of the MOT file time)
%         'FinalTime' - final time for ID to run until (defaults to
%                       time at end of the MOT file)
%         'ForcesToExclude' - List of forces by individual or grouping name 
%                       (e.g. All, actuators, muscles, ...) to be excluded 
%                       when computing model dynamics (default - 'Muscles')
%         'LowPassFilterForKinematics' - low pass filter cutoff frequency for 
%                       coordinate file data (default -1 = none)
%         'OutputFile' - string which of STO output file (default -
%                       'c3dfilename_id.sto')
%
% Output - ID (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  ID = setup_InverseDynamics('data',data, 'ModelFile','c:\Model.osim','MOTFile','c:\tiral1.mot',...
%               'GRFFile','c:\trial1_grf.xml','LowPassFilterFreq',10); 
% 
% Written by Glen Lichtwark (The University of Queensland)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy

ModelFile = [];
MOTFile = [];
GRFFile = [];

ResultsDirectory = [];

InitialTime = [];
FinalTime = [];

ForcesToExclude = 'Muscles'; 

LowPassFilterForKinematics = -1;
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

% define the initial and final times for inverse dynamics from the data structure 
% if this is passed to function and these aren't prescribed otherwise
% also define Output file if this isn't prescribed
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

%setup root and InverseDynamicsTool
root = 'OpenSimDocument';
        
V.ATTRIBUTE.Version = '20302';
[~, Name, ~] = fileparts(ModelFile);
[~, trial, ~] = fileparts(MOTFile);
V.InverseDynamicsTool.ATTRIBUTE.name = Name;

% define results directories
if isempty(ResultsDirectory)
    V.InverseDynamicsTool.results_directory = ' ';
else V.InverseDynamicsTool.results_directory = ResultsDirectory;
end

% define the model file
if ~isempty(ModelFile)
    V.InverseDynamicsTool.model_file = ModelFile;
end

V.InverseDynamicsTool.time_range = num2str([InitialTime FinalTime],6);
V.InverseDynamicsTool.forces_to_exclude = ForcesToExclude;

% define the GRF file for ID - this will be created an xml file in the
% setup process
V.InverseDynamicsTool.external_loads_file = GRFFile;

% define the MOT file for ID
V.InverseDynamicsTool.coordinates_file = MOTFile;

% filter
V.InverseDynamicsTool.lowpass_cutoff_frequency_for_coordinates = LowPassFilterForKinematics;

if isempty(OutputFile)
    V.InverseDynamicsTool.output_gen_force_file = [trial '_id.sto'];
else V.InverseDynamicsTool.output_gen_force_file = OutputFile;
end

fileout = [data.Name '_Setup_InverseDynamics.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root, Pref);

ID = V;
        