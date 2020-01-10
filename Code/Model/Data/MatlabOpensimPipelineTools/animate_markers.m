function data = animate_markers(data,a_data,fp_info,step)
% function data = animate_markers(data,a_data,fp_info,step)
% 
% This function will animate vicon c3d data that has been processed and
% display the markers in 3D space relative to the global coordiate and
% forceplates.
%
% Input - data - structure format containing data that has been loaded with
%                loadc3dfile.m and then sorted appropriately using the
%                sort_c3d.m function
%         a_data - structure containing the analog force plate data from 
%                the loadc3dfile.m function
%         fp_info - structure containing the coordinates of the force plate
%                corners, from the loadc3dfile.m function
%         step - number of frames to step over for animation (optional -
%                default = 10). This just speeds up animation.
%
% Output - data (optional) - outputs the data format which has additional
%                calculations force plate calculations (e.g. COP, Mz),  
%                downsampled to same frequency as the marker data
%
% Note, this function will animate specific data sets (PlugInGait,
%       OpenSimModel and OpenSimFullBodyModel) if the are the model names
%       in the data structure.
%                
% Authors: Rod Barrett and Glen Lichtwark (Griffith University and The 
%          University of Queensland)
%

if nargin < 4
    step = 5;
end

if nargin < 3
    fp_info = [];
end

if nargin < 2
    a_data = [];
end

switch data.Model
    case 'PlugInGait'
        if isfield(data.Markers,'LFHD')
        HeadX = [data.Markers.LFHD(:,1) data.Markers.RFHD(:,1) data.Markers.RBHD(:,1) data.Markers.LBHD(:,1) data.Markers.LFHD(:,1)];
        HeadY = [data.Markers.LFHD(:,2) data.Markers.RFHD(:,2) data.Markers.RBHD(:,2) data.Markers.LBHD(:,2) data.Markers.LFHD(:,2)];
        HeadZ = [data.Markers.LFHD(:,3) data.Markers.RFHD(:,3) data.Markers.RBHD(:,3) data.Markers.LBHD(:,3) data.Markers.LFHD(:,3)];
        else HeadX = zeros(size(data.Markers.LANK,1),1)*NaN;
            HeadY = zeros(size(data.Markers.LANK,1),1)*NaN;
            HeadZ = zeros(size(data.Markers.LANK,1),1)*NaN;
        end

        if isfield(data.Markers,'C7')
        TorsoX = [data.Markers.C7(:,1) data.Markers.T10(:,1) data.Markers.STRN(:,1) data.Markers.CLAV(:,1) data.Markers.C7(:,1)];
        TorsoY = [data.Markers.C7(:,2) data.Markers.T10(:,2) data.Markers.STRN(:,2) data.Markers.CLAV(:,2) data.Markers.C7(:,2)];
        TorsoZ = [data.Markers.C7(:,3) data.Markers.T10(:,3) data.Markers.STRN(:,3) data.Markers.CLAV(:,3) data.Markers.C7(:,3)];
        else TorsoX = zeros(size(data.Markers.LANK,1),1)*NaN;
            TorsoY = zeros(size(data.Markers.LANK,1),1)*NaN;
            TorsoZ = zeros(size(data.Markers.LANK,1),1)*NaN;
        end

        if isfield(data.Markers,'LASI')
        PelvisX = [data.Markers.LASI(:,1) data.Markers.RASI(:,1) data.Markers.RPSI(:,1) data.Markers.LPSI(:,1) data.Markers.LASI(:,1)];
        PelvisY = [data.Markers.LASI(:,2) data.Markers.RASI(:,2) data.Markers.RPSI(:,2) data.Markers.LPSI(:,2) data.Markers.LASI(:,2)];
        PelvisZ = [data.Markers.LASI(:,3) data.Markers.RASI(:,3) data.Markers.RPSI(:,3) data.Markers.LPSI(:,3) data.Markers.LASI(:,3)];
        else error('Need to have at least the leg markers to animate');
        end

        LeftLegX = [data.Markers.LASI(:,1) data.Markers.LTHI(:,1) data.Markers.LKNE(:,1) data.Markers.LTIB(:,1)...
            data.Markers.LANK(:,1) data.Markers.LTOE(:,1) data.Markers.LHEE(:,1) data.Markers.LANK(:,1)...
            data.Markers.LKNE(:,1) data.Markers.LASI(:,1)];
        LeftLegY = [data.Markers.LASI(:,2) data.Markers.LTHI(:,2) data.Markers.LKNE(:,2) data.Markers.LTIB(:,2)...
            data.Markers.LANK(:,2) data.Markers.LTOE(:,2) data.Markers.LHEE(:,2) data.Markers.LANK(:,2)...
            data.Markers.LKNE(:,2) data.Markers.LASI(:,2)];
        LeftLegZ = [data.Markers.LASI(:,3) data.Markers.LTHI(:,3) data.Markers.LKNE(:,3) data.Markers.LTIB(:,3)...
            data.Markers.LANK(:,3) data.Markers.LTOE(:,3) data.Markers.LHEE(:,3) data.Markers.LANK(:,3)...
            data.Markers.LKNE(:,3) data.Markers.LASI(:,3)];

        RightLegX = [data.Markers.RASI(:,1) data.Markers.RTHI(:,1) data.Markers.RKNE(:,1) data.Markers.RTIB(:,1)...
            data.Markers.RANK(:,1) data.Markers.RTOE(:,1) data.Markers.RHEE(:,1) data.Markers.RANK(:,1)...
            data.Markers.RKNE(:,1) data.Markers.RASI(:,1)];
        RightLegY = [data.Markers.RASI(:,2) data.Markers.RTHI(:,2) data.Markers.RKNE(:,2) data.Markers.RTIB(:,2)...
            data.Markers.RANK(:,2) data.Markers.RTOE(:,2) data.Markers.RHEE(:,2) data.Markers.RANK(:,2)...
            data.Markers.RKNE(:,2) data.Markers.RASI(:,2)];
        RightLegZ = [data.Markers.RASI(:,3) data.Markers.RTHI(:,3) data.Markers.RKNE(:,3) data.Markers.RTIB(:,3)...
            data.Markers.RANK(:,3) data.Markers.RTOE(:,3) data.Markers.RHEE(:,3) data.Markers.RANK(:,3)...
            data.Markers.RKNE(:,3) data.Markers.RASI(:,3)];

        if isfield(data,'Calc_Markers')
            if isfield(data.Calc_Markers,'LTIO')
            PantsX = [data.Calc_Markers.LTIO(:,1) data.Calc_Markers.LFEO(:,1) data.Calc_Markers.LFEP(:,1)...
                data.Calc_Markers.RFEP(:,1) data.Calc_Markers.RFEO(:,1) data.Calc_Markers.RTIO(:,1)];
            PantsY = [data.Calc_Markers.LTIO(:,2) data.Calc_Markers.LFEO(:,2) data.Calc_Markers.LFEP(:,2)...
                data.Calc_Markers.RFEP(:,2) data.Calc_Markers.RFEO(:,2) data.Calc_Markers.RTIO(:,2)];
            PantsZ = [data.Calc_Markers.LTIO(:,3) data.Calc_Markers.LFEO(:,3) data.Calc_Markers.LFEP(:,3)...
                data.Calc_Markers.RFEP(:,3) data.Calc_Markers.RFEO(:,3) data.Calc_Markers.RTIO(:,3)];

            else PantsX = zeros(size(data.Markers.LANK,1),1)*NaN;
            PantsY = zeros(size(data.Markers.LANK,1),1)*NaN;
            PantsZ = zeros(size(data.Markers.LANK,1),1)*NaN;
            end
            else PantsX = zeros(size(data.Markers.LANK,1),1)*NaN;
            PantsY = zeros(size(data.Markers.LANK,1),1)*NaN;
            PantsZ = zeros(size(data.Markers.LANK,1),1)*NaN;
        end

        if isfield(data,'ForcePlate')

            if length(data.ForcePlate(1).COP) > length(PelvisX)
                % if the forceplate and kinematic frequencies are different
                XI = 1:a_data.Rate/data.Rate:length(a_data.COPx1); 
            else XI = 1:length(PelvisX);
            end

            for i = 1:fp_info.Number
                Force(i).FPX = [fp_info.(['Corner1_' num2str(i)]).X fp_info.(['Corner2_' num2str(i)]).X fp_info.(['Corner3_' num2str(i)]).X fp_info.(['Corner4_' num2str(i)]).X fp_info.(['Corner1_' num2str(i)]).X];
                Force(i).FPY = [fp_info.(['Corner1_' num2str(i)]).Y fp_info.(['Corner2_' num2str(i)]).Y fp_info.(['Corner3_' num2str(i)]).Y fp_info.(['Corner4_' num2str(i)]).Y fp_info.(['Corner1_' num2str(i)]).Y];
                Force(i).FPZ = [fp_info.(['Corner1_' num2str(i)]).Z fp_info.(['Corner2_' num2str(i)]).Z fp_info.(['Corner3_' num2str(i)]).Z fp_info.(['Corner4_' num2str(i)]).Z fp_info.(['Corner1_' num2str(i)]).Z];
                Force(i).GRFX = interp1([data.ForcePlate(i).COP(:,1) data.ForcePlate(i).COP(:,1)+data.ForcePlate(i).F(:,1)],XI);
                Force(i).GRFY = interp1([data.ForcePlate(i).COP(:,2) data.ForcePlate(i).COP(:,2)+data.ForcePlate(i).F(:,2)],XI);
                Force(i).GRFZ = interp1([data.ForcePlate(i).COP(:,3) data.ForcePlate(i).COP(:,3)+data.ForcePlate(i).F(:,3)],XI);
            end

        end

        figure;
        set(gcf, 'doublebuffer', 'on');
        clear j

        for i = data.Start_Frame:step:data.End_Frame
            cla
            hold on
            for f = 1:fp_info.Number
                plot3(Force(f).FPX, Force(f).FPY, Force(f).FPZ, 'k', ...
                    Force(f).GRFX(i,:), Force(f).GRFY(i,:), Force(f).GRFZ(i,:),'m-');
            end
            plot3(TorsoX(i,:),TorsoY(i,:),TorsoZ(i,:),'-xk',...
                HeadX(i,:),HeadY(i,:),HeadZ(i,:),'-xk',...
                PelvisX(i,:),PelvisY(i,:),PelvisZ(i,:),'-xb',...
                LeftLegX(i,:),LeftLegY(i,:),LeftLegZ(i,:),'-xr',...
                RightLegX(i,:),RightLegY(i,:),RightLegZ(i,:),'-xg', ...
                PantsX(i,:), PantsY(i,:), PantsZ(i,:), 'k.-')
            axis equal
            axis([-1000 1000 -2000 2000 -1 2000])    
            line([0 300], [0 0], [0 0], 'Color', 'b', 'LineWidth', 3)
            line([0 0], [0 300], [0 0], 'Color', 'g', 'LineWidth', 3)
            line([0 0], [0 0], [0 300], 'Color', 'r', 'LineWidth', 3)
            xlabel('X'), ylabel('Y'), zlabel('Z')
            view(-60,30)
            drawnow
            %grid on
            hold off
        end
    case 'OpenSimModel'
        if isfield(data.Markers,'LFHD')
            HeadX = [data.Markers.LFHD(:,1) data.Markers.RFHD(:,1) data.Markers.RBHD(:,1) data.Markers.LBHD(:,1) data.Markers.LFHD(:,1)];
            HeadY = [data.Markers.LFHD(:,2) data.Markers.RFHD(:,2) data.Markers.RBHD(:,2) data.Markers.LBHD(:,2) data.Markers.LFHD(:,2)];
            HeadZ = [data.Markers.LFHD(:,3) data.Markers.RFHD(:,3) data.Markers.RBHD(:,3) data.Markers.LBHD(:,3) data.Markers.LFHD(:,3)];
            else HeadX = zeros(size(data.Markers.LANK,1),1)*NaN;
                HeadY = zeros(size(data.Markers.LANK,1),1)*NaN;
                HeadZ = zeros(size(data.Markers.LANK,1),1)*NaN;
            end

            if isfield(data.Markers,'C7')
            TorsoX = [data.Markers.C7(:,1) data.Markers.T10(:,1) data.Markers.STRN(:,1) data.Markers.CLAV(:,1) data.Markers.C7(:,1)];
            TorsoY = [data.Markers.C7(:,2) data.Markers.T10(:,2) data.Markers.STRN(:,2) data.Markers.CLAV(:,2) data.Markers.C7(:,2)];
            TorsoZ = [data.Markers.C7(:,3) data.Markers.T10(:,3) data.Markers.STRN(:,3) data.Markers.CLAV(:,3) data.Markers.C7(:,3)];
            else TorsoX = zeros(size(data.Markers.LANK,1),1)*NaN;
                TorsoY = zeros(size(data.Markers.LANK,1),1)*NaN;
                TorsoZ = zeros(size(data.Markers.LANK,1),1)*NaN;
            end

            if isfield(data.Markers,'LASI')
            PelvisX = [data.Markers.LASI(:,1) data.Markers.RASI(:,1) data.Markers.RPSI(:,1) data.Markers.LPSI(:,1) data.Markers.LASI(:,1)];
            PelvisY = [data.Markers.LASI(:,2) data.Markers.RASI(:,2) data.Markers.RPSI(:,2) data.Markers.LPSI(:,2) data.Markers.LASI(:,2)];
            PelvisZ = [data.Markers.LASI(:,3) data.Markers.RASI(:,3) data.Markers.RPSI(:,3) data.Markers.LPSI(:,3) data.Markers.LASI(:,3)];
            else error('Need to have at least the leg markers to animate');
            end

            if isfield(data.Markers,'LHJC')
            LeftLegX = [data.Markers.LHJC(:,1) data.Markers.LKJC(:,1) data.Markers.LAJC(:,1)];
            LeftLegY = [data.Markers.LHJC(:,2) data.Markers.LKJC(:,2) data.Markers.LAJC(:,2)];
            LeftLegZ = [data.Markers.LHJC(:,3) data.Markers.LKJC(:,3) data.Markers.LAJC(:,3)];

            RightLegX = [data.Markers.RHJC(:,1) data.Markers.RKJC(:,1) data.Markers.RAJC(:,1)];
            RightLegY = [data.Markers.RHJC(:,2) data.Markers.RKJC(:,2) data.Markers.RAJC(:,2)];
            RightLegZ = [data.Markers.RHJC(:,3) data.Markers.RKJC(:,3) data.Markers.RAJC(:,3)];
            end

            if isfield(data.Markers,'LHJC_func')
            LeftLegX = [data.Markers.LHJC_func(:,1) data.Markers.LKJC_func(:,1) data.Markers.LAJC(:,1)];
            LeftLegY = [data.Markers.LHJC_func(:,2) data.Markers.LKJC_func(:,2) data.Markers.LAJC(:,2)];
            LeftLegZ = [data.Markers.LHJC_func(:,3) data.Markers.LKJC_func(:,3) data.Markers.LAJC(:,3)];

            RightLegX = [data.Markers.RHJC_func(:,1) data.Markers.RKJC_func(:,1) data.Markers.RAJC(:,1)];
            RightLegY = [data.Markers.RHJC_func(:,2) data.Markers.RKJC_func(:,2) data.Markers.RAJC(:,2)];
            RightLegZ = [data.Markers.RHJC_func(:,3) data.Markers.RKJC_func(:,3) data.Markers.RAJC(:,3)];
            end

            if isfield(data.Markers,'LACR1')
            LeftArmX = [data.Markers.LACR1(:,1) data.Markers.LUA1(:,1) data.Markers.LElbow(:,1) ...
                data.Markers.LFA1(:,1) data.Markers.LWRR(:,1) data.Markers.LWRU(:,1) data.Markers.LCAR(:,1)];
            LeftArmY = [data.Markers.LACR1(:,2) data.Markers.LUA1(:,2) data.Markers.LElbow(:,2) ...
                data.Markers.LFA1(:,2) data.Markers.LWRR(:,2) data.Markers.LWRU(:,2) data.Markers.LCAR(:,2)];
            LeftArmZ = [data.Markers.LACR1(:,3) data.Markers.LUA1(:,3) data.Markers.LElbow(:,3) ...
                data.Markers.LFA1(:,3) data.Markers.LWRR(:,3) data.Markers.LWRU(:,3) data.Markers.LCAR(:,3)];

            RightArmX = [data.Markers.RACR1(:,1) data.Markers.RUA1(:,1) data.Markers.RElbow(:,1) ...
                data.Markers.RFA1(:,1) data.Markers.RWRR(:,1) data.Markers.RWRU(:,1) data.Markers.RCAR(:,1)];
            RightArmY = [data.Markers.RACR1(:,2) data.Markers.RUA1(:,2) data.Markers.RElbow(:,2) ...
                data.Markers.RFA1(:,2) data.Markers.RWRR(:,2) data.Markers.RWRU(:,2) data.Markers.RCAR(:,2)];
            RightArmZ = [data.Markers.RACR1(:,3) data.Markers.RUA1(:,3) data.Markers.RElbow(:,3) ...
                data.Markers.RFA1(:,3) data.Markers.RWRR(:,3) data.Markers.RWRU(:,3) data.Markers.RCAR(:,3)];
            end

            if isfield(data.Markers,'LMT1')
            LeftFootX = [data.Markers.LCAL(:,1) data.Markers.LMT1(:,1) data.Markers.LToe(:,1) ...
                data.Markers.LMT5(:,1) data.Markers.LCAL(:,1)];
            LeftFootY = [data.Markers.LCAL(:,2) data.Markers.LMT1(:,2) data.Markers.LToe(:,2) ...
                data.Markers.LMT5(:,2) data.Markers.LCAL(:,2)];
            LeftFootZ = [data.Markers.LCAL(:,3) data.Markers.LMT1(:,3) data.Markers.LToe(:,3) ...
                data.Markers.LMT5(:,3) data.Markers.LCAL(:,3)];

            RightFootX = [data.Markers.RCAL(:,1) data.Markers.RMT1(:,1) data.Markers.RToe(:,1) ...
                data.Markers.RMT5(:,1) data.Markers.RCAL(:,1)];
            RightFootY = [data.Markers.RCAL(:,2) data.Markers.RMT1(:,2) data.Markers.RToe(:,2) ...
                data.Markers.RMT5(:,2) data.Markers.RCAL(:,2)];
            RightFootZ = [data.Markers.RCAL(:,3) data.Markers.RMT1(:,3) data.Markers.RToe(:,3) ...
                data.Markers.RMT5(:,3) data.Markers.RCAL(:,3)];
            end

            MarkersX = [data.Markers.RTH1(:,1) data.Markers.RTH2(:,1) data.Markers.RTH3(:,1) data.Markers.RTH4(:,1) ...
                data.Markers.LTH1(:,1) data.Markers.LTH2(:,1) data.Markers.LTH3(:,1) data.Markers.LTH4(:,1) ...
                data.Markers.LTB1(:,1) data.Markers.LTB2(:,1) data.Markers.LTB3(:,1) data.Markers.LTB4(:,1) ...
                data.Markers.RTB1(:,1) data.Markers.RTB2(:,1) data.Markers.RTB3(:,1) data.Markers.RTB4(:,1) ...
                data.Markers.RBack(:,1)];
            MarkersY = [data.Markers.RTH1(:,2) data.Markers.RTH2(:,2) data.Markers.RTH3(:,2) data.Markers.RTH4(:,2) ...
                data.Markers.LTH1(:,2) data.Markers.LTH2(:,2) data.Markers.LTH3(:,2) data.Markers.LTH4(:,2) ...
                data.Markers.LTB1(:,2) data.Markers.LTB2(:,2) data.Markers.LTB3(:,2) data.Markers.LTB4(:,2) ...
                data.Markers.RTB1(:,2) data.Markers.RTB2(:,2) data.Markers.RTB3(:,2) data.Markers.RTB4(:,2) ...
                data.Markers.RBack(:,2)];
            MarkersZ = [data.Markers.RTH1(:,3) data.Markers.RTH2(:,3) data.Markers.RTH3(:,3) data.Markers.RTH4(:,3) ...
                data.Markers.LTH1(:,3) data.Markers.LTH2(:,3) data.Markers.LTH3(:,3) data.Markers.LTH4(:,3) ...
                data.Markers.LTB1(:,3) data.Markers.LTB2(:,3) data.Markers.LTB3(:,3) data.Markers.LTB4(:,3) ...
                data.Markers.RTB1(:,3) data.Markers.RTB2(:,3) data.Markers.RTB3(:,3) data.Markers.RTB4(:,3) ...
                data.Markers.RBack(:,3)];

            if isfield(data,'ForcePlate')
                if length(data.ForcePlate(1).COP) > length(MarkersX)
                    % if the forceplate and kinematic frequencies are different
                    XI = 1:a_data.Rate/data.Rate:length(a_data.COPx1); 
                else XI = 1:length(MarkersX);
                end

                for i = 1:fp_info.Number
                    Force(i).FPX = [fp_info.(['Corner1_' num2str(i)]).X fp_info.(['Corner2_' num2str(i)]).X fp_info.(['Corner3_' num2str(i)]).X fp_info.(['Corner4_' num2str(i)]).X fp_info.(['Corner1_' num2str(i)]).X];
                    Force(i).FPY = [fp_info.(['Corner1_' num2str(i)]).Y fp_info.(['Corner2_' num2str(i)]).Y fp_info.(['Corner3_' num2str(i)]).Y fp_info.(['Corner4_' num2str(i)]).Y fp_info.(['Corner1_' num2str(i)]).Y];
                    Force(i).FPZ = [fp_info.(['Corner1_' num2str(i)]).Z fp_info.(['Corner2_' num2str(i)]).Z fp_info.(['Corner3_' num2str(i)]).Z fp_info.(['Corner4_' num2str(i)]).Z fp_info.(['Corner1_' num2str(i)]).Z];
                    Force(i).GRFX = interp1([data.ForcePlate(i).COP(:,1) data.ForcePlate(i).COP(:,1)+data.ForcePlate(i).F(:,1)],XI);
                    Force(i).GRFY = interp1([data.ForcePlate(i).COP(:,2) data.ForcePlate(i).COP(:,2)+data.ForcePlate(i).F(:,2)],XI);
                    Force(i).GRFZ = interp1([data.ForcePlate(i).COP(:,3) data.ForcePlate(i).COP(:,3)+data.ForcePlate(i).F(:,3)],XI);
                end
            end

            figure;
            set(gcf, 'doublebuffer', 'on');

            for i = data.Start_Frame:step:data.End_Frame
                cla
                hold on
                for j = 1:fp_info.Number
                    plot3(Force(j).FPX, Force(j).FPY, Force(j).FPZ, 'k', ...
                        Force(j).GRFX(i,:), Force(j).GRFY(i,:), Force(j).GRFZ(i,:),'m-');
                end
                plot3(TorsoX(i,:),TorsoY(i,:),TorsoZ(i,:),'-xk',...
                    HeadX(i,:),HeadY(i,:),HeadZ(i,:),'-xk',...
                    PelvisX(i,:),PelvisY(i,:),PelvisZ(i,:),'-xb',...
                    LeftFootX(i,:),LeftFootY(i,:),LeftFootZ(i,:),'-xr',...
                    RightFootX(i,:),RightFootY(i,:),RightFootZ(i,:),'-xg', ...
                    LeftLegX(i,:),LeftLegY(i,:),LeftLegZ(i,:),'-xr', ...
                    RightLegX(i,:),RightLegY(i,:),RightLegZ(i,:),'-xg', ...    
                    LeftArmX(i,:),LeftArmY(i,:),LeftArmZ(i,:),'-xr', ...
                    RightArmX(i,:),RightArmY(i,:),RightArmZ(i,:),'-xg', ...
                    MarkersX(i,:),MarkersY(i,:),MarkersZ(i,:),'xk');
                hold off
                axis equal
                axis([-1000 1000 -2000 2000 -100 2000])    
                line([0 300], [0 0], [0 0], 'Color', 'b', 'LineWidth', 3)
                line([0 0], [0 300], [0 0], 'Color', 'g', 'LineWidth', 3)
                line([0 0], [0 0], [0 300], 'Color', 'r', 'LineWidth', 3)
                xlabel('X'), ylabel('Y'), zlabel('Z')
                view(-60,30)
                drawnow
                %grid on
            end
    otherwise 

        marker_names = fields(data.Markers);

        for i = 1:length(marker_names)
           Markers.X(:,i) = data.Markers.(char(marker_names(i)))(:,1);
           Markers.Y(:,i) = data.Markers.(char(marker_names(i)))(:,2);
           Markers.Z(:,i) = data.Markers.(char(marker_names(i)))(:,3);
        end

        if ~isempty(a_data) && ~isempty(fp_info)
            if isfield(data,'ForcePlate')
                if length(data.ForcePlate(1).COP) > length(Markers.X)
                    % if the forceplate and kinematic frequencies are different
                    XI = 1:a_data.Rate/data.Rate:length(a_data.COPx1); 
                else XI = 1:length(Markers.X);
                end

                for i = 1:fp_info.Number
                    Force(i).FPX = [fp_info.(['Corner1_' num2str(i)]).X fp_info.(['Corner2_' num2str(i)]).X fp_info.(['Corner3_' num2str(i)]).X fp_info.(['Corner4_' num2str(i)]).X fp_info.(['Corner1_' num2str(i)]).X];
                    Force(i).FPY = [fp_info.(['Corner1_' num2str(i)]).Y fp_info.(['Corner2_' num2str(i)]).Y fp_info.(['Corner3_' num2str(i)]).Y fp_info.(['Corner4_' num2str(i)]).Y fp_info.(['Corner1_' num2str(i)]).Y];
                    Force(i).FPZ = [fp_info.(['Corner1_' num2str(i)]).Z fp_info.(['Corner2_' num2str(i)]).Z fp_info.(['Corner3_' num2str(i)]).Z fp_info.(['Corner4_' num2str(i)]).Z fp_info.(['Corner1_' num2str(i)]).Z];
                    Force(i).GRFX = interp1([data.ForcePlate(i).COP(:,1) data.ForcePlate(i).COP(:,1)+data.ForcePlate(i).F(:,1)],XI);
                    Force(i).GRFY = interp1([data.ForcePlate(i).COP(:,2) data.ForcePlate(i).COP(:,2)+data.ForcePlate(i).F(:,2)],XI);
                    Force(i).GRFZ = interp1([data.ForcePlate(i).COP(:,3) data.ForcePlate(i).COP(:,3)+data.ForcePlate(i).F(:,3)],XI);
                end
            end
                end
        
        figure;
        set(gcf, 'doublebuffer', 'on');

        for i = data.Start_Frame:step:data.End_Frame
            if ~isempty(fp_info)
            for j = 1:fp_info.Number
                plot3(Force(j).FPX, Force(j).FPY, Force(j).FPZ, 'k', ...
                    Force(j).GRFX(i,:), Force(j).GRFY(i,:), Force(j).GRFZ(i,:),'m-');hold on;
            end
            end
            plot3(Markers.X(i,:), Markers.Y(i,:), Markers.Z(i,:), 'b.');  hold off;
            axis equal  
            axis([-1000 1000 -2000 2000 -100 2000])   
            %axis([min(min(Markers.X)) max(max(Markers.X)) min(min(Markers.Y)) max(max(Markers.Y)) min(min(Markers.Z)) max(max(Markers.Z))])      
            line([0 300], [0 0], [0 0], 'Color', 'b', 'LineWidth', 3)
            line([0 0], [0 300], [0 0], 'Color', 'g', 'LineWidth', 3)
            line([0 0], [0 0], [0 300], 'Color', 'r', 'LineWidth', 3)
            xlabel('X'), ylabel('Y'), zlabel('Z')
            view(-30,30)
            drawnow
        end
end
        
close(gcf)