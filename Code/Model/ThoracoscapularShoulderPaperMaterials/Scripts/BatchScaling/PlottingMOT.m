function []=PlottingMOT(x) % you need to define x='ABD_NL001_FP03_ik.mot' in
%the batch sclaing .m file

motData = dlmread(x, '\t', 11, 0);
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

%% plot IK result vs time

% fig1 showing thorax rotation & translation in x, y & z
figure(1);
subplot(2,1,1);
plot(time, ground_thorax_rot_x,'-k', time, ground_thorax_rot_y,':k', time, ground_thorax_rot_z,'--k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
legend('x-rot','y-rot','z-rot', 'Location','best');
title('Thorax Rotations & Translations','fontSize',12)
subplot(2,1,2);
plot(time, ground_thorax_tx,'-k', time, ground_thorax_ty,':k', time, ground_thorax_tz,'--k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Displacement (m)','FontWeight','bold');
legend('x-translation','y-translation','z-translation', 'Location','best');
savefig('Thorax_Rotations_Translations.fig');

% fig2 showing clavicle rotation & translation
figure(2);
plot(time, clav_prot,'-k',time, clav_elev,':k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
legend('Rotation','Elevation', 'Location','best');
title('Clavicle Rotations & Elevations','fontSize',12);
savefig('Clavicle_Rotations_Translations.fig');

% fig3 showing scapula rotations & translations
figure(3);
plot(time, scapula_abduction,'-k',time, scapula_elevation,'--k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
L1='Abduction';
L2='Elevation';
title('Scapula Rotations & Translations','fontSize',12);
hold on
plot(time, scapula_upward_rot,'color',[140 140 140]/255,'LineWidth',1);
L3='Upward-rot';
plot(time, scapula_winging,':k','LineWidth',1);
L4='Internal-rot';
legend(L1,L2,L3,L4,'Location','best');
hold off
savefig('Scapula_Rotations_Translations.fig');

% fig4 showing shoulder rotations and elevations.
%Plane elevation is any movements taken at the elevation plane, and
% the angle at the elevation plane is zero. Shoulder elevation is
% humeral or glenohumeral elevation(abduction). axial rotation is internal
% and external rotation of the humerus - best looked at when the elbow is
% slightly bent.
figure(4)
plot(time, plane_elv,'-k',time, shoulder_elv,':k', time, axial_rot,'--k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
title('Shoulder Rotations & Elevations','fontSize',12);
legend('Plane-elevation','Shoulder-elevation','Axial-rotation', 'Location','best');
savefig('Shoulder_Rotations_Elevations.fig');

% fig5 showing elbow flexion, pronation and supination at the writs.
figure(5)
subplot(2,1,1);
plot(time, elbow_flexion,'-k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
title('Elbow flexion, Ponation & Supination','fontSize',12);
legend('Flexion', 'Location','best');
subplot(2,1,2);
plot(time, pro_sup,'-k','LineWidth',1);
xlim([2.4 4.6]);
xlabel('Time (s)','FontWeight','bold');
ylabel('Joint angle (^o)','FontWeight','bold');
legend('Pronation & Supination', 'Location','best');
savefig('Elbow_flex_pronation_supination.fig');
end

