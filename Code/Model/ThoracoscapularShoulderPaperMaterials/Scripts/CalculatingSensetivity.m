% You might have to move this where the .mot files are stored
%time=0:0.01:2.2;
%Thorax
i=1;
for x=1:7
    K= num2str(x);
    Fname1=['S3_ABD_LOAD_ik_' K '.mot'];
    S3_ABD_LOAD{i} = dlmread(Fname1, '\t', 11, 0);    
    
    Fname2=['S3_ABD_NL_ik_' K '.mot'];
    S3_ABD_NL{i} = dlmread(Fname2, '\t', 11, 0);

    Fname3=['S3_FLEX_LOAD_ik_' K '.mot'];
    S3_FLEX_LOAD{i} = dlmread(Fname3, '\t', 11, 0);  
    
    Fname4=['S3_FLEX_NL_ik_' K '.mot'];
    S3_FLEX_NL{i} = dlmread(Fname4, '\t', 11, 0);
    
    Fname5=['S3_HTH_ik_' K '.mot'];
    S3_HTH{i} = dlmread(Fname5, '\t', 11, 0);
    
    Fname6=['S3_LAT_ROT_LOAD_ik_' K '.mot'];
    S3_LAT_ROT_LOAD{i} = dlmread(Fname6, '\t', 11, 0);
    
    Fname7=['S3_LAT_ROT_NL_ik_' K '.mot'];
    S3_LAT_ROT_NL{i} = dlmread(Fname7, '\t', 11, 0);
    i=i+1;
end

for j=1:7
Dif_ABD_L{j}=(S3_ABD_LOAD{1,7}-mean(S3_ABD_LOAD{1,7}))-(S3_ABD_LOAD{1,j}-mean(S3_ABD_LOAD{1,j}));
RMS_ABD_L{j}=sqrt(mean(Dif_ABD_L{1,j}.^2));

Dif_ABD_NL{j}=(S3_ABD_NL{1,7}-mean(S3_ABD_NL{1,7}))-(S3_ABD_NL{1,j}-mean(S3_ABD_NL{1,j}));
RMS_ABD_NL{j}=sqrt(mean(Dif_ABD_NL{1,j}.^2));

Dif_FLEX_L{j}=(S3_FLEX_LOAD{1,7}-mean(S3_FLEX_LOAD{1,7}))-(S3_FLEX_LOAD{1,j}-mean(S3_FLEX_LOAD{1,j}));
RMS_FLEX_L{j}=sqrt(mean(Dif_FLEX_L{1,j}.^2));

Dif_FLEX_NL{j}=(S3_FLEX_NL{1,7}(1:2100,:)-mean(S3_FLEX_NL{1,7}(1:2100,:)))-(S3_FLEX_NL{1,j}(1:2100,:)-mean(S3_FLEX_NL{1,j}(1:2100,:)));
RMS_FLEX_NL{j}=sqrt(mean(Dif_FLEX_NL{1,j}.^2));

Dif_HTH{j}=(S3_HTH{1,7}-mean(S3_HTH{1,7}))-(S3_HTH{1,j}-mean(S3_HTH{1,j}));
RMS_HTH{j}=sqrt(mean(Dif_HTH{1,j}.^2));

Dif_LAT_ROT_L{j}=(S3_LAT_ROT_LOAD{1,7}-mean(S3_LAT_ROT_LOAD{1,7}))-(S3_LAT_ROT_LOAD{1,j}-mean(S3_LAT_ROT_LOAD{1,j}));
RMS_LAT_ROT_L{j}=sqrt(mean(Dif_LAT_ROT_L{1,j}.^2));

Dif_LAT_ROT_NL{j}=(S3_LAT_ROT_NL{1,7}-mean(S3_LAT_ROT_NL{1,7}))-(S3_LAT_ROT_NL{1,j}-mean(S3_LAT_ROT_NL{1,j}));
RMS_LAT_ROT_NL{j}=sqrt(mean(Dif_LAT_ROT_NL{1,j}.^2));
end


RMS_ABD_L=reshape(cell2mat(RMS_ABD_L),18,7); 
RMS_ABD_NL=reshape(cell2mat(RMS_ABD_NL),18,7);
RMS_FLEX_L=reshape(cell2mat(RMS_FLEX_L),18,7);
RMS_FLEX_NL=reshape(cell2mat(RMS_FLEX_NL),18,7);
RMS_HTH=reshape(cell2mat(RMS_HTH),18,7);
RMS_LAT_ROT_L=reshape(cell2mat(RMS_LAT_ROT_L),18,7);
RMS_LAT_ROT_NL=reshape(cell2mat(RMS_LAT_ROT_NL),18,7);

%Thorax
g_thor_rotx=[RMS_ABD_L(2,1:6); RMS_ABD_NL(2,1:6); RMS_FLEX_L(2,1:6);...
    RMS_FLEX_NL(2,1:6); RMS_HTH(2,1:6); RMS_LAT_ROT_L(2,1:6); RMS_LAT_ROT_NL(2,1:6)];
g_thor_roty=[RMS_ABD_L(3,1:6); RMS_ABD_NL(3,1:6); RMS_FLEX_L(3,1:6);...
    RMS_FLEX_NL(3,1:6); RMS_HTH(3,1:6); RMS_LAT_ROT_L(3,1:6); RMS_LAT_ROT_NL(3,1:6)];
g_thor_rotz=[RMS_ABD_L(4,1:6); RMS_ABD_NL(4,1:6); RMS_FLEX_L(4,1:6);...
    RMS_FLEX_NL(4,1:6); RMS_HTH(4,1:6); RMS_LAT_ROT_L(4,1:6); RMS_LAT_ROT_NL(4,1:6)];
g_thor_tx=[RMS_ABD_L(5,1:6); RMS_ABD_NL(5,1:6); RMS_FLEX_L(5,1:6);...
    RMS_FLEX_NL(5,1:6); RMS_HTH(5,1:6); RMS_LAT_ROT_L(5,1:6); RMS_LAT_ROT_NL(5,1:6)];
g_thor_ty=[RMS_ABD_L(6,1:6); RMS_ABD_NL(6,1:6); RMS_FLEX_L(6,1:6);...
    RMS_FLEX_NL(6,1:6); RMS_HTH(6,1:6); RMS_LAT_ROT_L(6,1:6); RMS_LAT_ROT_NL(6,1:6)];
g_thor_tz=[RMS_ABD_L(7,1:6); RMS_ABD_NL(7,1:6); RMS_FLEX_L(7,1:6);...
    RMS_FLEX_NL(7,1:6); RMS_HTH(7,1:6); RMS_LAT_ROT_L(7,1:6); RMS_LAT_ROT_NL(7,1:6)];
g_thor_rotx_m=mean(g_thor_rotx);
g_thor_roty_m=mean(g_thor_roty);
g_thor_rotz_m=mean(g_thor_rotz);
g_thor_tx_m=mean(g_thor_tx);
g_thor_ty_m=mean(g_thor_ty);
g_thor_tz_m=mean(g_thor_tz);
g_thor_rotx_sd=std(g_thor_rotx);
g_thor_roty_sd=std(g_thor_roty);
g_thor_rotz_sd=std(g_thor_rotz);
g_thor_tx_sd=std(g_thor_tx);
g_thor_ty_sd=std(g_thor_ty);
g_thor_tz_sd=std(g_thor_tz);


%Clavicle
clav_prot=[RMS_ABD_L(8,1:6); RMS_ABD_NL(8,1:6); RMS_FLEX_L(8,1:6);...
    RMS_FLEX_NL(8,1:6); RMS_HTH(8,1:6); RMS_LAT_ROT_L(8,1:6); RMS_LAT_ROT_NL(8,1:6)];
clav_elev=[RMS_ABD_L(9,1:6); RMS_ABD_NL(9,1:6); RMS_FLEX_L(9,1:6);...
    RMS_FLEX_NL(9,1:6); RMS_HTH(9,1:6); RMS_LAT_ROT_L(9,1:6); RMS_LAT_ROT_NL(9,1:6)];
clav_prot_m=mean(clav_prot);
clav_elev_m=mean(clav_elev);
clav_prot_sd=std(clav_prot);
clav_elev_sd=std(clav_elev);

%Scapula
scapula_abduction=[RMS_ABD_L(10,1:6); RMS_ABD_NL(10,1:6); RMS_FLEX_L(10,1:6);...
    RMS_FLEX_NL(10,1:6); RMS_HTH(10,1:6); RMS_LAT_ROT_L(10,1:6); RMS_LAT_ROT_NL(10,1:6)];
scapula_elevation=[RMS_ABD_L(11,1:6); RMS_ABD_NL(11,1:6); RMS_FLEX_L(11,1:6);...
    RMS_FLEX_NL(11,1:6); RMS_HTH(11,1:6); RMS_LAT_ROT_L(11,1:6); RMS_LAT_ROT_NL(11,1:6)];
scapula_upward_rot=[RMS_ABD_L(12,1:6); RMS_ABD_NL(12,1:6); RMS_FLEX_L(12,1:6);...
    RMS_FLEX_NL(12,1:6); RMS_HTH(12,1:6); RMS_LAT_ROT_L(12,1:6); RMS_LAT_ROT_NL(12,1:6)];
scapula_winging=[RMS_ABD_L(13,1:6); RMS_ABD_NL(13,1:6); RMS_FLEX_L(13,1:6);...
    RMS_FLEX_NL(13,1:6); RMS_HTH(13,1:6); RMS_LAT_ROT_L(13,1:6); RMS_LAT_ROT_NL(13,1:6)];
scapula_abduction_m=mean(scapula_abduction);
scapula_elevation_m=mean(scapula_elevation);
scapula_upward_rot_m=mean(scapula_upward_rot);
scapula_winging_m=mean(scapula_winging);
scapula_abduction_sd=std(scapula_abduction);
scapula_elevation_sd=std(scapula_elevation);
scapula_upward_rot_sd=std(scapula_upward_rot);
scapula_winging_sd=std(scapula_winging);

%Shoulder
plane_elv=[RMS_ABD_L(14,1:6); RMS_ABD_NL(14,1:6); RMS_FLEX_L(14,1:6);...
    RMS_FLEX_NL(14,1:6); RMS_HTH(14,1:6); RMS_LAT_ROT_L(14,1:6); RMS_LAT_ROT_NL(14,1:6)];
shoulder_elv=[RMS_ABD_L(15,1:6); RMS_ABD_NL(15,1:6); RMS_FLEX_L(15,1:6);...
    RMS_FLEX_NL(15,1:6); RMS_HTH(15,1:6); RMS_LAT_ROT_L(15,1:6); RMS_LAT_ROT_NL(15,1:6)];
axial_rot=[RMS_ABD_L(16,1:6); RMS_ABD_NL(16,1:6); RMS_FLEX_L(16,1:6);...
    RMS_FLEX_NL(16,1:6); RMS_HTH(16,1:6); RMS_LAT_ROT_L(16,1:6); RMS_LAT_ROT_NL(16,1:6)];
plane_elv_m=mean(plane_elv);
shoulder_elv_m=mean(shoulder_elv);
axial_rot_m=mean(axial_rot);
plane_elv_sd=std(plane_elv);
shoulder_elv_sd=std(shoulder_elv);
axial_rot_sd=std(axial_rot);

%Elbow
elbow_flexion=[RMS_ABD_L(17,1:6); RMS_ABD_NL(17,1:6); RMS_FLEX_L(17,1:6);...
    RMS_FLEX_NL(17,1:6); RMS_HTH(17,1:6); RMS_LAT_ROT_L(17,1:6); RMS_LAT_ROT_NL(17,1:6)];
pro_sup=[RMS_ABD_L(18,1:6); RMS_ABD_NL(18,1:6); RMS_FLEX_L(18,1:6);...
    RMS_FLEX_NL(18,1:6); RMS_HTH(18,1:6); RMS_LAT_ROT_L(18,1:6); RMS_LAT_ROT_NL(18,1:6)];
elbow_flexion_m=mean(elbow_flexion);
pro_sup_m=mean(pro_sup);
elbow_flexion_sd=std(elbow_flexion);
pro_sup_sd=std(pro_sup);

M={'M1', 'M2', 'M3', 'M4','M5', 'M6'}; % 'M10', 'M11','M12','M13'};

%% 
% 1- Thorax
% Thorax_rot=figure;
% plot(1:6,g_thor_rotx_m(1,:),'or');
% hold on
% plot(1:6,g_thor_roty_m(1,:),'xb');
% plot(1:6,g_thor_rotz_m(1,:),'*k');
% grid on
% xticks([1:1:6]);
% xlim([0 7]);
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (^o)','FontWeight','bold')
% %yticks([0:0.1:1])
% %ylim([0 1])
% legend('Rotation-x','Rotation-y','Rotation-z', 'Location','best')
% title('Thorax Rotations','fontSize',12);
% %errorbar(g_thor_rotx_m(1,1:9),g_thor_rotx_sd(1:9))
% savefig('Thorax_rot.fig');
% saveas(Thorax_rot,'Thorax_rot.jpg')
% hold off

Thorax_rot_sub=figure;
subplot(3,1,1)
boxplot(g_thor_rotx(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Thorax Rotation-x')
%
subplot(3,1,2)
boxplot(g_thor_roty(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Thorax Rotation-y')
%
subplot(3,1,3)
boxplot(g_thor_rotz(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Thorax Rotation-z')
savefig('Thorax_rot_sub.fig')
saveas(Thorax_rot_sub,'Thorax_rot_sub.jpg')
%% 
% Thorax_t=figure;
% plot(1:9,g_thor_tx_m(1,1:9),'or');
% hold on
% plot(1:9,g_thor_ty_m(1,1:9),'xb');
% plot(1:9,g_thor_tz_m(1,1:9),'*k');
% grid on
% xticks([1:1:6]);
% xlim([0 7]);
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (m)','FontWeight','bold')
% yticks([0:0.0005:0.004])
% ylim([0 0.004])
% legend('Translation-x','Translation-y','Translation-z', 'Location','best')
% title('Thorax Translations','fontSize',12);
% savefig('Thorax_translation.fig');
% saveas(Thorax_t,'Thorax_t.jpg')
% hold off

Thorax_t_sub=figure;
subplot(3,1,1)
boxplot(g_thor_tx(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (m)','FontWeight','bold')
title('Thorax Translation-x')
%
subplot(3,1,2)
boxplot(g_thor_ty(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (m)','FontWeight','bold')
title('Thorax Translation-y')
%
subplot(3,1,3)
boxplot(g_thor_tz(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (m)','FontWeight','bold')
title('Thorax Translation-z')
savefig('Thorax_t_sub.fig')
saveas(Thorax_t_sub,'Thorax_t_sub.jpg')
%% 
% 2- Sternoclavicular joint
% Clavicle=figure;
% plot(1:9,clav_prot_m(1,1:9),'or');
% hold on
% plot(1:9,clav_elev_m(1,1:9),'xb');
% grid on
% xticks([1:1:6]);
% xlim([0 10])
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (^o)','FontWeight','bold')
% yticks([0:0.2:2])
% ylim([0 2])
% legend('Protraction & Retraction','Elevation', 'Location','best')
% title('Sternoclavicular Joint','fontSize',12);
% savefig('SternoclavicularJoint.fig');
% saveas(Clavicle,'Clavicle.jpg')
% hold off

Clavicle_sub=figure;
subplot(2,1,1)
boxplot(clav_prot(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Clavicle Protraction & Retraction')
%
subplot(2,1,2)
boxplot(clav_elev(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Clavicle Elevation')
savefig('Clavicle_sub.fig')
saveas(Clavicle_sub,'Clavicle_sub.jpg')
%% 
% 3- Scapulathoracic joint
% Scapula=figure;
% plot(1:9,scapula_abduction_m(1,1:9),'or');
% hold on
% plot(1:9,scapula_elevation_m(1,1:9),'xb');
% plot(1:9,scapula_upward_rot_m(1,1:9),'*k');
% plot(1:9,scapula_winging_m(1,1:9),'+m','LineWidth',1);
% grid on
% xticks([1:1:6]);
% xlim([0 7]);
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (^o)','FontWeight','bold')
% yticks([0:0.2:2.4])
% ylim([0 2.4])
% legend('Abduction','Elevation','Upward Rotation','Internal Rotation', 'Location','best')
% title('Scapulathoracic Joint','fontSize',12);
% savefig('ScapulathoracicJoint.fig');
% saveas(Scapula,'Scapula.jpg')
% hold off

Scapula_sub=figure;
subplot(4,1,1)
boxplot(scapula_abduction(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Scapula Abduction')
%
subplot(4,1,2)
boxplot(scapula_elevation(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Scapula Elevation')
%
subplot(4,1,3)
boxplot(scapula_upward_rot(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Scapula Upward Rotation')
%
subplot(4,1,4)
boxplot(scapula_winging(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Scapula Internal Rotation')
savefig('Scapula_sub.fig')
saveas(Scapula_sub,'Scapula_sub.jpg')
%% 
% 4- Glenohumeral joint
% Shoulder=figure;
% plot(1:9,plane_elv_m(1,1:9),'or');
% hold on
% plot(1:9,shoulder_elv_m(1,1:9),'xb');
% plot(1:9,axial_rot_m(1,1:9),'*k');
% grid on
% xticks([1:1:6]);
% xlim([0 7]);
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (^o)','FontWeight','bold')
% yticks([0:1:10])
% ylim([0 10])
% legend('Plane Elevation','Shoulder Elevation','Axial Rotation', 'Location','best')
% title('Glenohumeral Joint','fontSize',12);
% savefig('GlenohumeralJoint.fig');
% saveas(Shoulder,'Shoulder.jpg')
% hold off

Shoulder_sub=figure;
subplot(3,1,1)
boxplot(plane_elv(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Glenohumeral Joint Plane Elevation')
%
subplot(3,1,2)
boxplot(shoulder_elv(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Glenohumeral Joint Shoulder Elevation')
%
subplot(3,1,3)
boxplot(axial_rot(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Glenohumeral Joint Axial Rotation')
savefig('Shoulder_sub.fig')
saveas(Shoulder_sub,'Shoulder_sub.jpg')
%% 
% 5- Elbow joint
% Elbow=figure;
% plot(1:9,pro_sup_m(1,1:9),'or');
% hold on
% plot(1:9,elbow_flexion_m(1,1:9),'xb');
% grid on
% xticks([1:1:6]);
% xlim([0 7]);
% xticklabels(M)
% xtickangle(45)
% xlabel('Scaling Methods','FontWeight','bold')
% ylabel('RMS (^o)','FontWeight','bold')
% yticks([0:0.5:4])
% ylim([0 4])
% legend('Pronation & Supination','Elbow Flexion', 'Location','best')
% title('Elbow Joint','fontSize',12);
% savefig('ElbowJoint.fig');
% saveas(Elbow,'Elbow.jpg')
% hold off

Elbow_sub=figure;
subplot(2,1,1)
boxplot(pro_sup(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Pronation & Supination')
%
subplot(2,1,2)
boxplot(elbow_flexion(:,1:6))
grid on
xticks([1:1:6]);
xlim([0 7]);
xticklabels(M)
xtickangle(45)
%xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
title('Elbow Flexion')
savefig('Elbow_sub.fig');
saveas(Elbow_sub,'Elbow_sub.jpg')
