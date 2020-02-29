function Scale = setup_scale(varargin)

% function Scale = setup_scale(data)
%
% This function will create a Setup_Scale.XML file for use with OpenSIM
% based on the data from the C3D file and the names of the model files
% which will be used to scale the model.
% 
% Input - 'ModelFile' - string which is the filename of the OSIM file
%               (model file)
%         'MarkerFile' - filename string of the TRC file containing
%           coordinates of the markers to be used in scaling (defaults to
%           the data.TRC_Filename if 'data' structure passed)
%         'ScaleTasksFile' - filename string of the ScaleTasksFile
%
%          OPTIONAL INPUT FILES
%         'MarkerSetFile' - filename string of the MarkerSetFile
%         'ScaleSetFile' - filename string of the ScaleSetFile
%         'MeasurementSetFile' - filename string of the MMeasurementSetFile
%         'CoordinateFile' - filaname string of the CoordinateFile
%
%         OPTIONAL parameters
%         'data' - structure containing data from C3D file after processing
%                with C3D2TRC.m
%         'Mass' - subject mass (default - 75 or taken from 'data' structure if provided)
%         'Height' - subject height (defaule - 1750 or taken from 'data' structure if provided)
%         'Name' - subject name which is used as output model name
%               (defaults to model file name unless specified here or in
%               'data' structure)
%         'InitialTime' - time to start scaling from MOT file (defaults 0
%               if not specified or taken from time colum in 'data' structure
%               if provided as input)
%         'FinalTime' - finishing time for scaling from MOT file (defaults 1
%               if not specified or taken from time colum in 'data' structure
%               if provided as input)
%         'PreserveMass' - preserve the mass distribution (default = true)
%         'MaxMarkerMovement' - Maximum amount of movement allowed in marker 
%               data when averaging frames of the static trial. (Default = 
%               -1 --> no limit).
%         'OutputModelFile' - Name of the model to output (defaults to
%               subject name + '_SCALED.osim' if data structure sent with C3D info) 
%
% Output - Scale (optional) - This is the structure which is used to make the
%                         XML output file
%
% E.g.  Scale = setup_scale('data',data,'ModelSetFile','c:\models\model123.osim' ...
%              'MeasurementsSetFile', 'c:\models\model123_meas_set_file.xml')
% 
% Written by Glen Lichtwark (The University of Queensland)
%
% Inspirations from Tim Dorn's Gait Extract Toolbox -writeXML.m (University of
% Melbourne)

% setup default input files to emtpy

ModelFile = [];
MarkerFile = [];
Mass = 75;
Height = 1750;
Name = [];
InitialTime = [];
FinalTime = [];

MarkerSetFile = [];
MeasurementSetFile = [];
ScaleSetFile = [];
ScaleTasksFile = [];
CoordinateFile = '';
PreserveMass = 'true';
MaxMarkerMovement = -1;

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

if isempty(ModelFile)
    error('No model file input. Please input a model file to scale');
end

% determine subject name, mass and height from data structure if present
% otherwise set Name to the current model file name
% also define the initial and final times for the scaling if these aren't
% prescribe and the MarkerFile if this isn't prescribed
if ~isempty(data)
    if isfield(data,'Mass')
        Mass = data.Mass;
        Height = data.Height;
        Name = data.Name;
        if isempty(InitialTime)
            InitialTime = data.time(1);
        end
        if isempty(FinalTime)
            FinalTime = data.time(end);
        end
        if isempty(MarkerFile)
            MarkerFile = data.TRC_Filename;
        end
    end
else isempty(Name)
    [~, Name, ~] = fileparts(ModelFile);
end

if isempty(MarkerFile)
    error('MarkerFile - TRC file containing marker coordinates - not prescribed)');
end

%setup root and InverseKinematicsTool
root = 'OpenSimDocument';
        
V.ATTRIBUTE.Version = '20302';

%scale tool
V.ScaleTool.ATTRIBUTE.name = Name;

% ScaleTool -> GenericModelMaker
V.ScaleTool.GenericModelMaker.ATTRIBUTE.name = '';
if ~isempty(ModelFile)
    V.ScaleTool.GenericModelMaker.model_file = ModelFile;
end
if ~isempty(MarkerSetFile)
    V.ScaleTool.GenericModelMaker.marker_set_file = MarkerSetFile;
end

V.ScaleTool.mass = Mass;
V.ScaleTool.height = Height;

% ScaleTool -> ModelScaler
V.ScaleTool.ModelScaler.ATTRIBUTE.name = '';

V.ScaleTool.ModelScaler.apply = 'true';
if ~isempty(MeasurementSetFile)
    V.ScaleTool.ModelScaler.scaling_order = 'measurements';
    V.ScaleTool.ModelScaler.MeasurementSet.ATTRIBUTE.file = MeasurementSetFile;
end
if ~isempty(ScaleSetFile)
    V.ScaleTool.ModelScaler.scaling_order = 'manualScale';
    V.ScaleTool.ModelScaler.ScaleSet.ATTRIBUTE.file = ScaleSetFile;
end
if ~isempty(MeasurementSetFile) && ~isempty(ScaleSetFile)
    V.ScaleTool.ModelScaler.scaling_order = 'measurements manualScale';
end
V.ScaleTool.ModelScaler.marker_file = MarkerFile;
V.ScaleTool.ModelScaler.time_range = num2str([InitialTime FinalTime]);
V.ScaleTool.ModelScaler.preserve_mass_distribution = PreserveMass;
% V.ModelScaler.output_model_file = [data.Name '_SCALED_NO_MARKERS.osim'];
V.ScaleTool.ModelScaler.output_scale_file = [Name '_ScaleSet_Applied.xml'];

% ScaleTool -> MarkerPlacer
if ~isempty(ScaleTasksFile)
    V.ScaleTool.MarkerPlacer.ATTRIBUTE.name = '';

    V.ScaleTool.MarkerPlacer.IKTaskSet.ATTRIBUTE.file = ScaleTasksFile;
    V.ScaleTool.MarkerPlacer.marker_file = MarkerFile;
    V.ScaleTool.MarkerPlacer.coordinate_file = CoordinateFile;
    V.ScaleTool.MarkerPlacer.time_range = num2str([InitialTime FinalTime]);
    if isempty(OutputModelFile)
        V.ScaleTool.MarkerPlacer.output_model_file = [Name, '_SCALED.osim'];
    else V.ScaleTool.MarkerPlacer.output_model_file = OutputModelFile;
    end
    V.ScaleTool.MarkerPlacer.output_motion_file = [Name, '_static_output.mot'];
    V.ScaleTool.MarkerPlacer.max_marker_movement = MaxMarkerMovement;
end

fileout = [data.Name '_Setup_Scale.xml'];

Pref.StructItem = false;

xml_write(fileout, V, root,Pref);

Scale = V;
        