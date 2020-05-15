%% combining RMS results for all seven activities from the two subjects
load('RMS_Subject2.mat')
load('RMS_Subject3.mat')
n2=1;
n3=1;
for i=1:12
    if rem(i,2)==0 %the second column is subject 3 - in magenta
        RMS_axial_rot(:,i)=axial_rot_3(:,n2);
        RMS_clav_elev(:,i)=clav_elev_3(:,n2);
        RMS_clav_prot(:,i)=clav_prot_3(:,n2);
        RMS_elbow_flex(:,i)=elbow_flexion_3(:,n2);
        RMS_g_thor_rotx(:,i)=g_thor_rotx_3(:,n2);
        RMS_g_thor_roty(:,i)=g_thor_roty_3(:,n2);
        RMS_g_thor_rotz(:,i)=g_thor_rotz_3(:,n2);
        RMS_g_thor_tx(:,i)=g_thor_tx_3(:,n2);
        RMS_g_thor_ty(:,i)=g_thor_ty_3(:,n2);
        RMS_g_thor_tz(:,i)=g_thor_tz_3(:,n2);
        RMS_plane_elv(:,i)=plane_elv_3(:,n2);
        RMS_pro_sup(:,i)=pro_sup_3(:,n2);
        RMS_scap_abd(:,i)=scapula_abduction_3(:,n2);
        RMS_scap_elev(:,i)=scapula_elevation_3(:,n2);
        RMS_scap_upward_rot(:,i)=scapula_upward_rot_3(:,n2);
        RMS_scap_wing(:,i)=scapula_winging_3(:,n2);
        RMS_shoulder_elev(:,i)=shoulder_elv_3(:,n2);        
        n2=n2+1;
        
    else %the first column is subject 2 - in blue   
        RMS_axial_rot(:,i)=axial_rot_2(:,n3);
        RMS_clav_elev(:,i)=clav_elev_2(:,n3);
        RMS_clav_prot(:,i)=clav_prot_2(:,n3);
        RMS_elbow_flex(:,i)=elbow_flexion_2(:,n3);
        RMS_g_thor_rotx(:,i)=g_thor_rotx_2(:,n3);
        RMS_g_thor_roty(:,i)=g_thor_roty_2(:,n3);
        RMS_g_thor_rotz(:,i)=g_thor_rotz_2(:,n3);
        RMS_g_thor_tx(:,i)=g_thor_tx_2(:,n3);
        RMS_g_thor_ty(:,i)=g_thor_ty_2(:,n3);
        RMS_g_thor_tz(:,i)=g_thor_tz_2(:,n3);
        RMS_plane_elv(:,i)=plane_elv_2(:,n3);
        RMS_pro_sup(:,i)=pro_sup_2(:,n3);
        RMS_scap_abd(:,i)=scapula_abduction_2(:,n3);
        RMS_scap_elev(:,i)=scapula_elevation_2(:,n3);
        RMS_scap_upward_rot(:,i)=scapula_upward_rot_2(:,n3);
        RMS_scap_wing(:,i)=scapula_winging_2(:,n3);
        RMS_shoulder_elev(:,i)=shoulder_elv_2(:,n3);        
        n3=n3+1;
        
    end
end
%% plotting the combined data as a box plot


% 1- Thorax 1
thorax_rot=figure;
subplot(3,1,1)
formatedboxplot(RMS_g_thor_rotx)
title('Thorax Rotation x')
subplot(3,1,2)
formatedboxplot(RMS_g_thor_roty)
title('Thorax Rotation y')
subplot(3,1,3)
formatedboxplot(RMS_g_thor_rotz)
title('Thorax Rotation z')

% 2- Thorax 2
thorax_trans=figure;
subplot(3,1,1)
formatedboxplot(RMS_g_thor_tx)
title('Thorax Translation x')
subplot(3,1,2)
formatedboxplot(RMS_g_thor_ty)
title('Thorax Translation y')
subplot(3,1,3)
formatedboxplot(RMS_g_thor_tz)
title('Thorax Translation z')

% 3- Sternoclavicular joint
Clavicle=figure;
subplot(2,1,1)
formatedboxplot(RMS_clav_prot)
title('Clavicle Protraction & Retraction')
subplot(2,1,2)
formatedboxplot(RMS_clav_elev)
title('Clavicle Elevation')

% 4- Scapulathoracic joint 1 of 2
Scapula1=figure;
subplot(2,1,1)
formatedboxplot(RMS_scap_abd)
title('Scapula Abduction')
subplot(2,1,2)
formatedboxplot(RMS_scap_elev)
title('Scapula Elevation')

% 5- Scapulathoracic joint 2 of 2
Scapula2=figure;
subplot(2,1,1)
formatedboxplot(RMS_scap_upward_rot)
title('Scapula Upward Rotation')
subplot(2,1,2)
formatedboxplot(RMS_scap_wing)
title('Scapula Internal Rotation')

% 6- Glenohumeral joint
Shoulder=figure;
subplot(3,1,1)
formatedboxplot(RMS_plane_elv)
title('Glenohumeral Joint Plane Elevation')
subplot(3,1,2)
formatedboxplot(RMS_shoulder_elev)
title('Glenohumeral Joint Shoulder Elevation')
subplot(3,1,3)
formatedboxplot(RMS_axial_rot)
title('Glenohumeral Joint Axial Rotation')

% 7- Elbow joint
Elbow=figure;
subplot(2,1,1)
formatedboxplot(RMS_pro_sup)
title('Pronation & Supination')
subplot(2,1,2)
formatedboxplot(RMS_elbow_flex)
title('Elbow Flexion')

savefig(thorax_rot,'thorax_rot.fig')
savefig(thorax_trans,'thorax_trans.fig')
savefig(Clavicle,'Clavicle.fig')
savefig(Scapula1,'Scapula.fig')
savefig(Scapula2,'Scapula.fig')
savefig(Shoulder,'Elbow.fig')
savefig(Elbow,'Elbow.fig')

saveas(thorax_rot,'thorax_rot.jpg')
saveas(thorax_trans,'thorax_trans.jpg')
saveas(Clavicle,'Clavicle.jpg')
saveas(Scapula1,'Scapula.jpg')
saveas(Scapula2,'Scapula.jpg')
saveas(Shoulder,'Elbow.jpg')
saveas(Elbow,'Elbow.jpg')
