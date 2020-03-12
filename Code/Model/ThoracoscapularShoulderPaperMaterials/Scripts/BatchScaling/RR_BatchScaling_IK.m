% Pull in the modeling and scaling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% Load Plugin needed for opening model
Model.LoadOpenSimLibrary('C:\Users\rsk02\Documents\GitHub\TechForParalysis\Code\Model\ThoracoscapularShoulderPaperMaterials\ScapulothoracicJointPlugin40\WinX64\ScapulothoracicJointPlugin40_WinX64');

% Get the generic model
[modelFile,modelFilePath] = uigetfile('*.osim','Pick the the model file to be used.');

%Load the original model and initialize
model_0 = Model(fullfile(modelFilePath, modelFile));
model_0.initSystem;
%% 
% Go to the subject's folder where .trc files are located
trc_data_folder = 'C:\Users\rsk02\Documents\GitHub\TechForParalysis\Code\Model\ThoracoscapularShoulderPaperMaterials\ExperimentalData\Markers';

% Get and operate on the files
[SetupForIK,SetupPath] = uigetfile('*.xml','Pick the IK.xml configuration file.');
ikTool = InverseKinematicsTool([SetupPath SetupForIK]);

%% 
[filename, filepath] = uigetfile('*.xml','Pick the scale.xml configuration file.');

% Select scale tool and apply it to the xml file
SCTool = ScaleTool([filepath filename]);

load('ScaleFactors.mat'); % load the data file for s.f.
i=1;
for ntest=4:17 % loop to go through all trials for xml
ij_c7(i)=ScaleFactor(1,ntest);                 % x thorax
ij_px(i)=ScaleFactor(2,ntest);                 % y thorax
ij_ac(i)=ScaleFactor(3,ntest);                 % z thorax
ijc_acc_x(i)=ScaleFactor(4,ntest);             % x clavicle
ijc_acc_y(i)=ScaleFactor(5,ntest);             % y clavicle
ijc_acc_z(i)=ScaleFactor(6,ntest);             % z clavicle
aa_ac(i)=ScaleFactor(7,ntest);                 % x scapula
ai_ts(i)=ScaleFactor(8,ntest);                 % y scapula
aa_ts(i)=ScaleFactor(9,ntest);                 % z scapula
EpL_EpM_x(i)=ScaleFactor(10,ntest);            % x humerus
gu_centelbow(i)=ScaleFactor(11,ntest);         % y humerus
EpL_EpM_z(i)=ScaleFactor(12,ntest);            % z humerus
rs_us_x(i)=ScaleFactor(13,ntest);              % x ulna & x radius
centelbow_centusrs(i)=ScaleFactor(14,ntest);   % y ulna, y radius
rs_us_z(i)=ScaleFactor(15,ntest);              % z ulna & z radius


% Getting to measures under model scalar which is under scale tool
ms = SCTool.getModelScaler().getScaleSet();

% change thorax s.f.
a_t = ms.get('thorax_xyz');
prop_t = a_t.updPropertyByIndex(0);
PropertyHelper.setValueDouble(ij_c7(i), prop_t, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ij_px(i), prop_t, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(ij_ac(i), prop_t, 2) % changing the z-axis s.f.

% change clavicle s.f.
a_c = ms.get('clavicle_xyz');
prop_c = a_c.updPropertyByIndex(0);
PropertyHelper.setValueDouble(ijc_acc_x(i), prop_c, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ijc_acc_y(i), prop_c, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(ijc_acc_z(i), prop_c, 2) % changing the z-axis s.f.

% change scapula s.f.
a_s = ms.get('scapula_xyz');
prop_s = a_s.updPropertyByIndex(0);
PropertyHelper.setValueDouble(aa_ac(i), prop_s, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ai_ts(i), prop_s, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(aa_ts(i), prop_s, 2) % changing the z-axis s.f.

% change humerus s.f.
a_h = ms.get('humerus_xyz');
prop_h = a_h.updPropertyByIndex(0);
PropertyHelper.setValueDouble(EpL_EpM_x(i), prop_h, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(gu_centelbow(i), prop_h, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(EpL_EpM_z(i), prop_h, 2) % changing the z-axis s.f.

% change ulna s.f.
a_u = ms.get('ulna_xyz');
prop_u = a_u.updPropertyByIndex(0);
PropertyHelper.setValueDouble(rs_us_x(i), prop_u, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(centelbow_centusrs(i), prop_u, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(rs_us_z(i), prop_u, 2) % changing the z-axis s.f.

% change radius s.f.
a_r = ms.get('radius_xyz');
prop_r = a_r.updPropertyByIndex(0);
PropertyHelper.setValueDouble(rs_us_x(i), prop_r, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(centelbow_centusrs(i), prop_r, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(rs_us_z(i), prop_r, 2) % changing the z-axis s.f.

% change hand s.f.
a_ha = ms.get('hand_xyz');
prop_ha = a_ha.updPropertyByIndex(0);
PropertyHelper.setValueDouble(1, prop_ha, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 2) % changing the z-axis s.f.

% set file name & path for the generic model
GM=SCTool.getGenericModelMaker();
GM_path=GM.updPropertyByIndex(0);
PropertyHelper.setValueString('C:\Users\rsk02\Documents\GitHub\TechForParalysis\Code\Model\ThoracoscapularShoulderPaperMaterials\Model\ThoracoscapularShoulderModel02_generic.osim', GM_path, 0)

% set file name & path for scaled .osim file
ScaledModelPath='C:\Users\rsk02\Documents\GitHub\TechForParalysis\Code\Model\ThoracoscapularShoulderPaperMaterials\Model\ScaledModels\';
NewNameOSIM=['Model',int2str(i),'.osim'];
PathOSIM=[ScaledModelPath,NewNameOSIM];
OM_1=SCTool.getModelScaler();
OM_2=OM_1.updPropertyByIndex(7);
PropertyHelper.setValueString(PathOSIM,OM_2,0);

% name & print xml scaling files
NewNameXML=['ScaleProtocol_',int2str(i),'.xml'];
SCTool.print(NewNameXML);

% run xml scaling file & as a result it will create .osim file
K=SCTool.ScaleTool.run();

% Specify IK results folder
resultsfolder = ['IK output_',int2str(i)];
mkdir(resultsfolder);

% IK for the experimental lab tests (i.e. elbow extension, flexion, etc)
% Tell Tool to use the loaded model
model = Model(fullfile(ScaledModelPath, NewNameOSIM));
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



i=i+1;
end


