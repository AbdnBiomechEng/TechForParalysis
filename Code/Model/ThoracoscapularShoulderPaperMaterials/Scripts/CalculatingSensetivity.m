% You might have to move this where the .mot files are stored
time=0:0.01:2.2;
%Thorax
i=1;
for x=0:13
    K= num2str(x);
    Fname=['motData' K];
    load(Fname)
    
    %Thorax
    ground_thorax_rot_x(:,i)=motData(241:461,2);
    ground_thorax_rot_y(:,i)=motData(241:461,3);
    ground_thorax_rot_z(:,i)=motData(241:461,4);
    ground_thorax_tx(:,i)=motData(241:461,5);
    ground_thorax_ty(:,i)=motData(241:461,6);
    ground_thorax_tz(:,i)=motData(241:461,7);
    
    %Clavicle
    clav_prot(:,i)=motData(241:461,8);
    clav_elev(:,i)=motData(241:461,9);

    %Scapula
    scapula_abduction(:,i)=motData(241:461,10);
    scapula_elevation(:,i)=motData(241:461,11);
    scapula_upward_rot(:,i)=motData(241:461,12);
    scapula_winging(:,i)=motData(241:461,13);

    %Shoulder
    plane_elv(:,i)=motData(241:461,14);
    shoulder_elv(:,i)=motData(241:461,15);
    axial_rot(:,i)=motData(241:461,16);

    %Elbow
    elbow_flexion(:,i)=motData(241:461,17);

    %Wrist
    pro_sup(:,i)=motData(241:461,18);
    
    
    i=i+1;
end

%Thorax
Dif_rot_x=(ground_thorax_rot_x(:,1)-ground_thorax_rot_x);
Dif_rot_y=(ground_thorax_rot_y(:,1)-ground_thorax_rot_y);
Dif_rot_z=(ground_thorax_rot_z(:,1)-ground_thorax_rot_z);
RMS_Dif_rot_x=sqrt(mean(Dif_rot_x.^2));
RMS_Dif_rot_y=sqrt(mean(Dif_rot_y.^2));
RMS_Dif_rot_z=sqrt(mean(Dif_rot_z.^2));

Dif_t_x=(ground_thorax_rot_x(:,1)-ground_thorax_rot_x);
Dif_t_y=(ground_thorax_rot_y(:,1)-ground_thorax_rot_y);
Dif_t_z=(ground_thorax_rot_z(:,1)-ground_thorax_rot_z);
RMS_Dif_t_x=sqrt(mean(Dif_t_x.^2));
RMS_Dif_t_y=sqrt(mean(Dif_t_y.^2));
RMS_Dif_t_z=sqrt(mean(Dif_t_z.^2));
RMS_Thorax=[RMS_Dif_rot_x; RMS_Dif_rot_y; RMS_Dif_rot_z;...
    RMS_Dif_t_x; RMS_Dif_t_y; RMS_Dif_t_z];

%Clavicle
Dif_clav_prot =clav_prot(:,1)-clav_prot;
Dif_clav_elev =clav_elev(:,1)-clav_elev;
RMS_Dif_clav_prot=sqrt(mean(Dif_clav_prot.^2));
RMS_Dif_clav_elev=sqrt(mean(Dif_clav_elev.^2));
RMS_Clavicle=[RMS_Dif_clav_prot; RMS_Dif_clav_elev];

%Scapula
Dif_scapula_abduction=(scapula_abduction(:,1)-scapula_abduction);
Dif_scapula_elevation=(scapula_elevation(:,1)-scapula_elevation);
Dif_scapula_upward_rot=(scapula_upward_rot(:,1)-scapula_upward_rot);
Dif_scapula_winging=(scapula_winging(:,1)-scapula_winging);
RMS_Dif_scapula_abduction=sqrt(mean(Dif_scapula_abduction.^2));
RMS_Dif_scapula_elevation=sqrt(mean(Dif_scapula_elevation.^2));
RMS_Dif_scapula_upward_rot=sqrt(mean(Dif_scapula_upward_rot.^2));
RMS_Dif_scapula_winging=sqrt(mean(Dif_scapula_winging.^2));
RMS_Scapula=[RMS_Dif_scapula_abduction; RMS_Dif_scapula_elevation;...
    RMS_Dif_scapula_upward_rot; RMS_Dif_scapula_winging];

%Shoulder
Dif_plane_elv=(plane_elv(:,1)-plane_elv);
Dif_shoulder_elv=(shoulder_elv(:,1)-shoulder_elv);
Dif_axial_rot=(axial_rot(:,1)-axial_rot);
RMS_Dif_plane_elv=sqrt(mean(Dif_plane_elv.^2));
RMS_Dif_shoulder_elv=sqrt(mean(Dif_shoulder_elv.^2));
RMS_Dif_axial_rot=sqrt(mean(Dif_axial_rot.^2));
RMS_Shoulder=[RMS_Dif_plane_elv; RMS_Dif_shoulder_elv; ...
    RMS_Dif_axial_rot];

%Elbow
Dif_elbow_flexion=(elbow_flexion(:,1)-elbow_flexion);
RMS_Dif_elbow_flexion=sqrt(mean(Dif_elbow_flexion.^2));

%Wrist
Dif_pro_sup=(pro_sup(:,1)-pro_sup);
RMS_Dif_pro_sup=sqrt(mean(Dif_pro_sup.^2));


