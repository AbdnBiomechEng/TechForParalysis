function setupAndRunIKBatch(trc_data_folder,results_folder,model_name)

% Based on setupAndRunIKBatchExample.m by Edith Arnold
% trc_data_folder is the folder that contains the .trc files
% results_folder is the folder where the results will be printed
% model_name is the .osim model used for IK

% Pull in the modeling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% Choose a generic setup file to work from
ikTool = InverseKinematicsTool('IK_setup.xml');

% Get the model
model = Model(model_name);
model.initSystem();

% Tell Tool to use the loaded model
ikTool.setModel(model);

trialsForIK = dir(fullfile(trc_data_folder, '*.trc'));

nTrials = size(trialsForIK);

% Loop through the trials
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
    ikTool.setOutputMotionFileName(fullfile(results_folder, [name '_ik.mot']));
    
    % Save the settings in a setup file
    outfile = ['Setup_IK_' name '.xml'];
    ikTool.print(outfile);
    
    fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
    % Run IK
    ikTool.run();

end
