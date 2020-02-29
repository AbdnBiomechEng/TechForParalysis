function btk_animate_markers(marker_data, fp_data, step_size)
% function btk_animate_markers(marker_data, fp_data, step)
% 
% This function will animate vicon c3d data that has been processed and
% display the markers in 3D space relative to the global coordiate and
% forceplates.
%
% Input - marker_data - structure format containing data that has been 
%                loaded with BTK_LOADC3D function
%         grw_data - structure containing the analog force plate data from 
%                the loadc3dfile.m function
%         fp_data - structure containing the coordinates of the force plate
%                corners, from the loadc3dfile.m function
%         step_size - number of frames to step over for animation (optional -
%                default = 10). This just speeds up animation.
%               
% Author: Glen Lichtwark (The University of Queensland)
% Updated: 06/12/2011

if nargin < 3
    step_size = 5;
end

if nargin < 2
    fp_data = [];
end

if nargin < 1
    error('No data passed to function');
end

if ~isempty(fp_data)
    for i = 1:length(fp_data.FP_data)
        Force(i).FPX = [fp_data.FP_data(i).corners(1,:)'; fp_data.FP_data(i).corners(1,1)];
        Force(i).FPY = [fp_data.FP_data(i).corners(2,:)'; fp_data.FP_data(i).corners(2,1)];
        Force(i).FPZ = [fp_data.FP_data(i).corners(3,:)'; fp_data.FP_data(i).corners(3,1)];
        Force(i).GRFX = interp1(fp_data.Time, [fp_data.GRF_data(i).P(:,1) fp_data.GRF_data(i).P(:,1)+fp_data.GRF_data(i).F(:,1)], marker_data.Time);
        Force(i).GRFY = interp1(fp_data.Time, [fp_data.GRF_data(i).P(:,2) fp_data.GRF_data(i).P(:,2)+fp_data.GRF_data(i).F(:,2)], marker_data.Time);
        Force(i).GRFZ = interp1(fp_data.Time, [fp_data.GRF_data(i).P(:,3) fp_data.GRF_data(i).P(:,3)+fp_data.GRF_data(i).F(:,3)], marker_data.Time);
    end
end

marker_names = fieldnames(marker_data.Markers);
% tf = strncmp('C',marker_names,1);
% marker_names(tf) = [];

for i = 1:length(marker_names)
   Markers.X(:,i) = marker_data.Markers.(marker_names{i})(:,1);
   Markers.Y(:,i) = marker_data.Markers.(marker_names{i})(:,2);
   Markers.Z(:,i) = marker_data.Markers.(marker_names{i})(:,3);
end

figure;
set(gcf, 'doublebuffer', 'on');

for i = 1:step_size:marker_data.Info.NumFrames
    if ~isempty(fp_data)
        for j = 1:length(Force)
            plot3(Force(j).FPX, Force(j).FPY, Force(j).FPZ, 'k', ...
                Force(j).GRFX(i,:), Force(j).GRFY(i,:), Force(j).GRFZ(i,:),'m-');hold on;
        end
    end
    plot3(Markers.X(i,:), Markers.Y(i,:), Markers.Z(i,:), 'b.');  hold off;
    axis equal  
    axis([min(min(Markers.X)) max(max(Markers.X)) min(min(Markers.Y)) max(max(Markers.Y)) min(min(Markers.Z)) max(max(Markers.Z))])      
    line([0 300], [0 0], [0 0], 'Color', 'b', 'LineWidth', 3)
    line([0 0], [0 300], [0 0], 'Color', 'g', 'LineWidth', 3)
    line([0 0], [0 0], [0 300], 'Color', 'r', 'LineWidth', 3)
    xlabel('X'), ylabel('Y'), zlabel('Z')
%     view(-90,10)
    drawnow
    pause(0.075)
    %grid on
end

close(gcf)