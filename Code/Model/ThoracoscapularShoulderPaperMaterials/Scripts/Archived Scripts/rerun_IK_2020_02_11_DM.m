%% INVERSE KINEMATICS (Batch Processing)

% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% Load Plugin needed for opening model
Model.LoadOpenSimLibrary('C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\ScapulothoracicJointPlugin40\WinX64\ScapulothoracicJointPlugin40_WinX64');

% Get the model
[modelFile,modelFilePath] = uigetfile('*.osim','Pick the the model file to be used.');

% Load the model and initialize
model = Model(fullfile(modelFilePath, modelFile));
model.initSystem();

% Go to the subject's folder where .trc files are located
trc_data_folder = 'C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\ExperimentalData\Markers\FP';

% Get and operate on the files
[SetupForIK,SetupPath] = uigetfile('*.xml','Pick the IK.xml setup file for this subject/model.');
ikTool = InverseKinematicsTool([SetupPath SetupForIK]);

% Specify results folder
resultsfolder = 'output_rerun_2020_02_11';
mkdir(resultsfolder);

% Tell Tool to use the loaded model
ikTool.setModel(model);
trialsForIK = dir(fullfile(trc_data_folder, '*.trc'));
nTrials = size(trialsForIK);

% Loop through IK for all the trials
for trial= 1:nTrials
    
    % Get the name of the file for this trial
    markerFile = trialsForIK(trial).name;
    
    % Create name of trial from .trc file name
    name = regexprep(markerFile,'.trc','');
    fullpath = fullfile(trc_data_folder, markerFile);
    
    % Get trc data to determine time range
    markerData = MarkerData(fullpath);
    
    % Get initial and intial time 
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

TotalTime = sum(time)/60;
disp(['Total time = ', num2str(TotalTime), ' minutes'])


