
% Pull in the modeling and scaling classes straight from the OpenSim distribution
import org.opensim.modeling.*

% Load Plugin needed for opening model
Model.LoadOpenSimLibrary('C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\ScapulothoracicJointPlugin40\WinX64\ScapulothoracicJointPlugin40_WinX64');

% Get the original model
[modelFile,modelFilePath] = uigetfile('*.osim','Pick the the model file to be used.');

% Load the original model and initialize
model_0 = Model(fullfile(modelFilePath, modelFile));
model_0.initSystem;

%% 
[filename, filepath] = uigetfile('*.xml','Pick the scale.xml setup file.');

% Select scale tool and apply it to the xml file
%SCTool = ScaleTool(fullfile(filepath, filename));
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


%Getting to measures under model scalar which is under scale tool
ms = SCTool.getModelScaler().getScaleSet();
a_t = ms.get('thorax');
prop_t = a_t.updPropertyByIndex(0);
PropertyHelper.setValueDouble(ij_c7(i), prop_t, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ij_px(i), prop_t, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(ij_ac(i), prop_t, 2) % changing the z-axis s.f.

a_c = ms.get('clavicle');
prop_c = a_c.updPropertyByIndex(0);
PropertyHelper.setValueDouble(ijc_acc_x(i), prop_c, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ijc_acc_y(i), prop_c, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(ijc_acc_z(i), prop_c, 2) % changing the z-axis s.f.

a_s = ms.get('scapula');
prop_s = a_s.updPropertyByIndex(0);
PropertyHelper.setValueDouble(aa_ac(i), prop_s, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(ai_ts(i), prop_s, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(aa_ts(i), prop_s, 2) % changing the z-axis s.f.

a_h = ms.get('humerus');
prop_h = a_h.updPropertyByIndex(0);
PropertyHelper.setValueDouble(EpL_EpM_x(i), prop_h, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(gu_centelbow(i), prop_h, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(EpL_EpM_z(i), prop_h, 2) % changing the z-axis s.f.

a_u = ms.get('ulna');
prop_u = a_u.updPropertyByIndex(0);
PropertyHelper.setValueDouble(rs_us_x(i), prop_u, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(centelbow_centusrs(i), prop_u, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(rs_us_z(i), prop_u, 2) % changing the z-axis s.f.

a_r = ms.get('radius');
prop_r = a_r.updPropertyByIndex(0);
PropertyHelper.setValueDouble(rs_us_x(i), prop_r, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(centelbow_centusrs(i), prop_r, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(rs_us_z(i), prop_r, 2) % changing the z-axis s.f.

a_ha = ms.get('hand');
prop_ha = a_ha.updPropertyByIndex(0);
PropertyHelper.setValueDouble(1, prop_ha, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 2) % changing the z-axis s.f.

GM=SCTool.getGenericModelMaker();
GM_path=GM.updPropertyByIndex(0);
PropertyHelper.setValueString('C:\Users\rsk02\Desktop\IK routine\ThoracoscapularShoulderPaperMaterials\Model\ThoracoscapularShoulderModel02.osim', GM_path, 0)

% set file name & path in the .osim to 
OM_1=SCTool.getModelScaler();
OM_2=OM_1.updPropertyByIndex(7);
PropertyHelper.setValueString('C:\Users\rsk02\Desktop\IK routine\BatchScaling\modelTEST.osim',OM_2,0);

% set file name & path to output model after scaling and marker placement
% OM_1=SCTool.getMarkerPlacer();
% OM_2=OM_1.updPropertyByIndex(6);

NewNameXML=['ScaleProtocol_',int2str(i),'.xml'];
SCTool.print(NewNameXML);

K=SCTool.ScaleTool.run();


i=i+1;
end


