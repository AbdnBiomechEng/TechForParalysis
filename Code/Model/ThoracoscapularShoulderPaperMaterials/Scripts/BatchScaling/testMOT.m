motData = dlmread('ABD_NL001_FP03_ik.mot', '\t', 11, 0);
time=motData(:,1);
ground_thorax_rot_x=motData(:,2);
ground_thorax_rot_y=motData(:,3);
ground_thorax_rot_z=motData(:,4);
ground_thorax_tx=motData(:,5);
ground_thorax_ty=motData(:,6);
ground_thorax_tz=motData(:,7);
clav_prot=motData(:,8);
clav_elev=motData(:,9);
scapula_abduction=motData(:,10);
scapula_elevation=motData(:,11);
scapula_upward_rot=motData(:,12);
scapula_winging=motData(:,13);
plane_elv=motData(:,14);
shoulder_elv=motData(:,15);
axial_rot=motData(:,16);
elbow_flexion=motData(:,17);
pro_sup=motData(:,18);

% plot IK result vs time 
figure (1); % fig1 showing thorax rotation in x, y & z
plot(time, ground_thorax_rot_x,'r', time, ground_thorax_rot_y,'g', time, ground_thorax_rot_z,'b');
xlabel('Time (s)');
ylabel('IK (degrees)');
legend('x-rot','y-rot','z-rot');
title('Thorax')



