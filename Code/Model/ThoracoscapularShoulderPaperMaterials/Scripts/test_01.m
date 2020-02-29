% setting up the data - ScaleFactor consists of model & experimental length
% between two markers & the corresponding s.f.
ij_c7=ScaleFactor(1,:);                 % x thorax
ij_px=ScaleFactor(2,:);                 % y thorax
ij_ac=ScaleFactor(3,:);                 % z thorax
c7_t8=ScaleFactor(4,:);                 % y-2 thorax
ijc_acc=ScaleFactor(5,:);               % y clavicle
aa_ac=ScaleFactor(6,:);                 % x scapula
ai_ts=ScaleFactor(7,:);                 % y scapula
aa_ts=ScaleFactor(8,:);                 % z scapula
EpL_EpM=ScaleFactor(9,:);               % x & z humerus
gu_centelbow=ScaleFactor(10,:);         % y humerus
rs_us=ScaleFactor(11,:);                % x & z ulna, x & z radius
centelbow_centusrs=ScaleFactor(12,:);   % y ulna, y radius

% import OpenSim modeling library
import org.opensim.modeling.*
[filename, filepath] = uigetfile('*.xml');

% Select scale tool and apply it to the xml file
SCTool = ScaleTool(fullfile(filepath, filename));

%Getting to measures under model scalar which is under scale tool
ms = SCTool.getModelScaler().getScaleSet();
a_t = ms.get('thorax');
prop_t = a_t.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_t, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_t, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_t, 2) % changing the z-axis s.f.

a_c = ms.get('clavicle');
prop_c = a_c.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_c, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_c, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_c, 2) % changing the z-axis s.f.

a_s = ms.get('scapula');
prop_s = a_s.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_s, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_s, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_s, 2) % changing the z-axis s.f.

a_h = ms.get('humerus');
prop_h = a_h.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_h, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_h, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_h, 2) % changing the z-axis s.f.

a_u = ms.get('ulna');
prop_u = a_u.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_u, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_u, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_u, 2) % changing the z-axis s.f.

a_r = ms.get('radius');
prop_r = a_r.updPropertyByIndex(0);
PropertyHelper.setValueDouble(500, prop_r, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(600, prop_r, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(700, prop_r, 2) % changing the z-axis s.f.

a_ha = ms.get('hand');
prop_ha = a_ha.updPropertyByIndex(0);
PropertyHelper.setValueDouble(1, prop_ha, 0) % changing the x-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 1) % changing the y-axis s.f.
PropertyHelper.setValueDouble(1, prop_ha, 2) % changing the z-axis s.f.

SCTool.print('test.xml')
%% 

