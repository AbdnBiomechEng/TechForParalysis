%% Scaling (Batch Processing)

% Pull in the modeling and scaling classes straight from the OpenSim distribution
import org.opensim.modeling.*
import org.opensim.modeling.Scale
import org.opensim.modeling.ScaleSet
import org.opensim.modeling.ScaleTool
import org.opensim.modeling.ModelScaler

% Load Plugin needed for opening model
Model.LoadOpenSimLibrary('C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\ScapulothoracicJointPlugin40\WinX64\ScapulothoracicJointPlugin40_WinX64');

% Get the model
[modelFile,modelFilePath] = uigetfile('*.osim','Pick the the model file to be used.');

% Load the original model and initialize
model1 = Model(fullfile(modelFilePath, modelFile));
model1.initSystem;

% Go to the subject's folder where .trc files are located
trc_data_folder = 'C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\ExperimentalData\Markers\FP';
xml_data_folder = 'C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\Simulations\IK_RR';

trialsForScale=dir(fullfile(xml_data_folder, '*.xml'));
nTest=size(trialsForScale);

% Loop to scale all the models and carry out IK for each
for test=1:nTest

% Get the name of the file for this xml test
scaleFile = trialsForScale(test).name;

% Create name of test from .xml file name
name = regexprep(scaleFile,'.xml','');
fullpath = fullfile(xml_data_folder, scaleFile);

% Setup the ScaleTool for this test
[SetupForScale,SetupScalePath] = uigetfile('*.xml','Pick the scale.xml setup file for this subject/model.');
ScTool = ScaleTool([SetupScalePath SetupForScale]);
ScTool.ScaleTool.run()

% Tell Tool to use the loaded model
ikTool.setModel(model1);
trialsForIK = dir(fullfile(trc_data_folder, '*.trc'));
nTrials = size(trialsForIK);
    
for trial= 1:nTrials
    % Get the name of the file for this trial
    markerFile = trialsForIK(trial).name;
    
    % Create name of trial from .trc file name
    name = regexprep(markerFile,'.trc','');
    fullpath = fullfile(trc_data_folder, markerFile);
    
    % Get trc data to determine time range
    markerData = MarkerData(fullpath);
    
    % Get initial and final time 
    initial_time = markerData.getStartFrameTime();
    final_time = markerData.getLastFrameTime();
    
    % Setup the ikTool for this trial
    ikTool.setName(name);
    ikTool.setMarkerDataFileName(fullpath);
    ikTool.setStartTime(initial_time);
    ikTool.setEndTime(final_time);
    ikTool.setOutputMotionFileName(fullfile([cd '\' resultsfolder], [name '_ik.mot']));
    
    % Save the settings in a setup file
    outfile = ['Setup_IK_' name '.xml'];
    ikTool.print(fullfile(SetupPath, outfile));
    
    % Printing out current cycle
    fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
    
    % Run IK
    tic;
    ikTool.run();
    time(trial) = toc;

    % Move Marker Errors + Locations + Setup File into Output Folder
    [status, message] = movefile('*_ik_model_marker_locations.sto', [cd '\' resultsfolder]);
    if(status ~= 1)
        disp(['Locations file FAILED at relocating because ' message]);
    end
    
    [status, message] = movefile('*_ik_marker_errors.sto', [cd '\' resultsfolder]);
    if(status ~= 1)
        disp(['Errors file FAILED at relocating because ' message]);
    end
    
    [status, message] = movefile('Setup_IK_*.xml', [cd '\' resultsfolder]);
    if(status ~= 1)
        disp(['Setup file FAILED at relocating because ' message]);
    end
    
end
end

TotalTime = sum(time)/60;
disp(['Total time = ', num2str(TotalTime), ' minutes'])


