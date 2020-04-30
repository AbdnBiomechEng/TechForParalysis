% You might have to move this where the .mot files are stored
%time=0:0.01:2.2;
%Thorax
i=1;
for x=1:14
    K= num2str(x);
    Fname1=['S2_ABD_LOAD_ik_' K '.mot'];
    S2_ABD_LOAD{i} = dlmread(Fname1, '\t', 11, 0);    
    
    Fname2=['S2_ABD_NL_ik_' K '.mot'];
    S2_ABD_NL{i} = dlmread(Fname2, '\t', 11, 0);

    Fname3=['S2_FLEX_LOAD_ik_' K '.mot'];
    S2_FLEX_LOAD{i} = dlmread(Fname3, '\t', 11, 0);  
    
    Fname4=['S2_FLEX_NL_ik_' K '.mot'];
    S2_FLEX_NL{i} = dlmread(Fname4, '\t', 11, 0);
    
    Fname5=['S2_HTH_ik_' K '.mot'];
    S2_HTH{i} = dlmread(Fname5, '\t', 11, 0);
    
    Fname6=['S2_LAT_ROT_LOAD_ik_' K '.mot'];
    S2_LAT_ROT_LOAD{i} = dlmread(Fname6, '\t', 11, 0);
    
    Fname7=['S2_LAT_ROT_NL_ik_' K '.mot'];
    S2_LAT_ROT_NL{i} = dlmread(Fname7, '\t', 11, 0);
    i=i+1;
end

for j=1:14
Dif_ABD_L{j}=S2_ABD_LOAD{1,10}-S2_ABD_LOAD{1,j};
RMS_ABD_L{j}=sqrt(mean(Dif_ABD_L{1,j}.^2));

Dif_ABD_NL{j}=S2_ABD_NL{1,10}-S2_ABD_NL{1,j};
RMS_ABD_NL{j}=sqrt(mean(Dif_ABD_NL{1,j}.^2));

Dif_FLEX_L{j}=S2_FLEX_LOAD{1,10}-S2_FLEX_LOAD{1,j};
RMS_FLEX_L{j}=sqrt(mean(Dif_FLEX_L{1,j}.^2));

Dif_FLEX_NL{j}=S2_FLEX_NL{1,10}-S2_FLEX_NL{1,j};
RMS_FLEX_NL{j}=sqrt(mean(Dif_FLEX_NL{1,j}.^2));

Dif_HTH{j}=S2_HTH{1,10}-S2_HTH{1,j};
RMS_HTH{j}=sqrt(mean(Dif_HTH{1,j}.^2));

Dif_LAT_ROT_L{j}=S2_LAT_ROT_LOAD{1,10}-S2_LAT_ROT_LOAD{1,j};
RMS_LAT_ROT_L{j}=sqrt(mean(Dif_LAT_ROT_L{1,j}.^2));

Dif_LAT_ROT_NL{j}=S2_LAT_ROT_NL{1,10}-S2_LAT_ROT_NL{1,j};
RMS_LAT_ROT_NL{j}=sqrt(mean(Dif_LAT_ROT_NL{1,j}.^2));


end

RMS_ABD_L=reshape(cell2mat(RMS_ABD_L),18,14); 
RMS_ABD_NL=reshape(cell2mat(RMS_ABD_NL),18,14);
RMS_FLEX_L=reshape(cell2mat(RMS_FLEX_L),18,14);
RMS_FLEX_NL=reshape(cell2mat(RMS_FLEX_NL),18,14);
RMS_HTH=reshape(cell2mat(RMS_HTH),18,14);
RMS_LAT_ROT_L=reshape(cell2mat(RMS_LAT_ROT_L),18,14);
RMS_LAT_ROT_NL=reshape(cell2mat(RMS_LAT_ROT_NL),18,14);

%Thorax
g_thor_rotx=[RMS_ABD_L(2,:); RMS_ABD_NL(2,:); RMS_FLEX_L(2,:);...
    RMS_FLEX_NL(2,:); RMS_HTH(2,:); RMS_LAT_ROT_L(2,:); RMS_LAT_ROT_NL(2,:)];
g_thor_roty=[RMS_ABD_L(3,:); RMS_ABD_NL(3,:); RMS_FLEX_L(3,:);...
    RMS_FLEX_NL(3,:); RMS_HTH(3,:); RMS_LAT_ROT_L(3,:); RMS_LAT_ROT_NL(3,:)];
g_thor_rotz=[RMS_ABD_L(4,:); RMS_ABD_NL(4,:); RMS_FLEX_L(4,:);...
    RMS_FLEX_NL(4,:); RMS_HTH(4,:); RMS_LAT_ROT_L(4,:); RMS_LAT_ROT_NL(4,:)];
g_thor_tx=[RMS_ABD_L(5,:); RMS_ABD_NL(5,:); RMS_FLEX_L(5,:);...
    RMS_FLEX_NL(5,:); RMS_HTH(5,:); RMS_LAT_ROT_L(5,:); RMS_LAT_ROT_NL(5,:)];
g_thor_ty=[RMS_ABD_L(6,:); RMS_ABD_NL(6,:); RMS_FLEX_L(6,:);...
    RMS_FLEX_NL(6,:); RMS_HTH(6,:); RMS_LAT_ROT_L(6,:); RMS_LAT_ROT_NL(6,:)];
g_thor_tz=[RMS_ABD_L(7,:); RMS_ABD_NL(7,:); RMS_FLEX_L(7,:);...
    RMS_FLEX_NL(7,:); RMS_HTH(7,:); RMS_LAT_ROT_L(7,:); RMS_LAT_ROT_NL(7,:)];
g_thor_rotx_m=mean(g_thor_rotx);
g_thor_roty_m=mean(g_thor_roty);
g_thor_rotz_m=mean(g_thor_rotz);
g_thor_tx_m=mean(g_thor_tx);
g_thor_ty_m=mean(g_thor_ty);
g_thor_tz_m=mean(g_thor_tz);


%Clavicle
clav_prot=[RMS_ABD_L(8,:); RMS_ABD_NL(8,:); RMS_FLEX_L(8,:);...
    RMS_FLEX_NL(8,:); RMS_HTH(8,:); RMS_LAT_ROT_L(8,:); RMS_LAT_ROT_NL(8,:)];
clav_elev=[RMS_ABD_L(9,:); RMS_ABD_NL(9,:); RMS_FLEX_L(9,:);...
    RMS_FLEX_NL(9,:); RMS_HTH(9,:); RMS_LAT_ROT_L(9,:); RMS_LAT_ROT_NL(9,:)];
clav_prot_m=mean(clav_prot);
clav_elev_m=mean(clav_elev);

%Scapula
scapula_abduction=[RMS_ABD_L(10,:); RMS_ABD_NL(10,:); RMS_FLEX_L(10,:);...
    RMS_FLEX_NL(10,:); RMS_HTH(10,:); RMS_LAT_ROT_L(10,:); RMS_LAT_ROT_NL(10,:)];
scapula_elevation=[RMS_ABD_L(11,:); RMS_ABD_NL(11,:); RMS_FLEX_L(11,:);...
    RMS_FLEX_NL(11,:); RMS_HTH(11,:); RMS_LAT_ROT_L(11,:); RMS_LAT_ROT_NL(11,:)];
scapula_upward_rot=[RMS_ABD_L(12,:); RMS_ABD_NL(12,:); RMS_FLEX_L(12,:);...
    RMS_FLEX_NL(12,:); RMS_HTH(12,:); RMS_LAT_ROT_L(12,:); RMS_LAT_ROT_NL(12,:)];
scapula_winging=[RMS_ABD_L(13,:); RMS_ABD_NL(13,:); RMS_FLEX_L(13,:);...
    RMS_FLEX_NL(13,:); RMS_HTH(13,:); RMS_LAT_ROT_L(13,:); RMS_LAT_ROT_NL(13,:)];
scapula_abduction_m=mean(scapula_abduction);
scapula_elevation_m=mean(scapula_elevation);
scapula_upward_rot_m=mean(scapula_upward_rot);
scapula_winging_m=mean(scapula_winging);

%Shoulder
plane_elv=[RMS_ABD_L(14,:); RMS_ABD_NL(14,:); RMS_FLEX_L(14,:);...
    RMS_FLEX_NL(14,:); RMS_HTH(14,:); RMS_LAT_ROT_L(14,:); RMS_LAT_ROT_NL(14,:)];
shoulder_elv=[RMS_ABD_L(15,:); RMS_ABD_NL(15,:); RMS_FLEX_L(15,:);...
    RMS_FLEX_NL(15,:); RMS_HTH(15,:); RMS_LAT_ROT_L(15,:); RMS_LAT_ROT_NL(15,:)];
axial_rot=[RMS_ABD_L(16,:); RMS_ABD_NL(16,:); RMS_FLEX_L(16,:);...
    RMS_FLEX_NL(16,:); RMS_HTH(16,:); RMS_LAT_ROT_L(16,:); RMS_LAT_ROT_NL(16,:)];
plane_elv_m=mean(plane_elv);
shoulder_elv_m=mean(shoulder_elv);
axial_rot_m=mean(axial_rot);

%Elbow
elbow_flexion=[RMS_ABD_L(17,:); RMS_ABD_NL(17,:); RMS_FLEX_L(17,:);...
    RMS_FLEX_NL(17,:); RMS_HTH(17,:); RMS_LAT_ROT_L(17,:); RMS_LAT_ROT_NL(17,:)];
pro_sup=[RMS_ABD_L(18,:); RMS_ABD_NL(18,:); RMS_FLEX_L(18,:);...
    RMS_FLEX_NL(18,:); RMS_HTH(18,:); RMS_LAT_ROT_L(18,:); RMS_LAT_ROT_NL(18,:)];
elbow_flexion_m=mean(elbow_flexion);
pro_sup_m=mean(pro_sup);

M={'M1', 'M2', 'M3', 'M4','M5', 'M6', 'M7', 'M8','M9', 'M10', 'M11',...
    'M12','M13'};

Thorax_rot=figure;
plot(1:13,[g_thor_rotx_m(1,1:9) g_thor_rotx_m(1,11:14)],'or');
hold on
plot(1:13,[g_thor_roty_m(1,1:9) g_thor_roty_m(1,11:14)],'xb');
plot(1:13,[g_thor_rotz_m(1,1:9) g_thor_rotz_m(1,11:14)],'*k');
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
legend('Rotation-x','Rotation-y','Rotation-z', 'Location','best')
title('Thorax Rotations','fontSize',12);
%savefig('Thorax_rot.fig');
hold off

Thorax_t=figure;
plot(1:13,[g_thor_tx_m(1,1:9) g_thor_tx_m(1,11:14)],'or');
hold on
plot(1:13,[g_thor_ty_m(1,1:9) g_thor_ty_m(1,11:14)],'xb');
plot(1:13,[g_thor_tz_m(1,1:9) g_thor_tz_m(1,11:14)],'*k');
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (m)','FontWeight','bold')
legend('Translation-x','Translation-y','Translation-z', 'Location','best')
title('Thorax Translation','fontSize',12);
%savefig('Thorax_translation.fig');
hold off

Clavicle=figure;
plot(1:13,[clav_prot_m(1,1:9) clav_prot_m(1,11:14)],'or');
hold on
plot(1:13,[clav_elev_m(1,1:9) clav_elev_m(1,11:14)],'xb');
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
legend('Rotation','Elevation', 'Location','best')
title('Clavicle','fontSize',12);
%savefig('Clavicle.fig');
hold off

Scapula=figure;
plot(1:13,[scapula_abduction_m(1,1:9) scapula_abduction_m(1,11:14)],'or');
hold on
plot(1:13,[scapula_elevation_m(1,1:9) scapula_elevation_m(1,11:14)],'xb');
plot(1:13,[scapula_upward_rot_m(1,1:9) scapula_upward_rot_m(1,11:14)],'*k');
plot(1:13,[scapula_winging_m(1,1:9) scapula_winging_m(1,11:14)],'+m','LineWidth',1);
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
legend('Abduction','Elevation','Upward Rotation','Internal Rotation', 'Location','best')
title('Scapula','fontSize',12);
%savefig('Thorax_rot.fig');
hold off

Shoulder=figure;
plot(1:13,[plane_elv_m(1,1:9) plane_elv_m(1,11:14)],'or');
hold on
plot(1:13,[shoulder_elv_m(1,1:9) shoulder_elv_m(1,11:14)],'xb');
plot(1:13,[axial_rot_m(1,1:9) axial_rot_m(1,11:14)],'*k');
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
legend('Plane-elevation','Shoulder-elevation','Axial-rotation', 'Location','best')
title('Shoulder','fontSize',12);
%savefig('Shoulder.fig');
hold off

Elbow_Fig=figure;
plot(1:13,[pro_sup_m(1,1:9) pro_sup_m(1,11:14)],'or');
hold on
plot(1:13,[elbow_flexion_m(1,1:9) elbow_flexion_m(1,11:14)],'xb');
grid on
xticks([1:1:13]);
xticklabels(M)
xtickangle(45)
xlabel('Scaling Methods','FontWeight','bold')
ylabel('RMS (^o)','FontWeight','bold')
legend('Pronation & Supination','Elbow Flexion', 'Location','best')
title('Elbow','fontSize',12);
%savefig('Elbow.fig');
hold off



