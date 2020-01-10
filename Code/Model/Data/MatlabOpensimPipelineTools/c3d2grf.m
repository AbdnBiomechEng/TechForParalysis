function data = c3d2grf(data,a_data,fp_info,subsample)
% This function will use the information from the forceplates and the
% analog force plate data to create a GRF vector in the global coordinate
% system, sampled at the same frequency as the kinematic marker data. All
% COP's, forces and GRF moments are expressed relative to the global
% coordinate system (not the force plate system)
%
% Input - data - structure format containing data that has been loaded with
%                loadc3dfile.m 
%         a_data - structure containing the analog force plate data from 
%                the loadc3dfile.m function
%         fp_info - structure containing the coordinates of the force plate
%                corners and origin coordinates, from the loadc3dfile.m function
%         subsample - 'on' or 'off', whether to subsample forceplate data at same
%                frequency as kinematic data (default = 'on')
%       
% Output - data - structure format containing orginal data plus the new
%                force plate data in a field called Forceplate
%
% Author: Glen Lichtwark (Griffith University)
%

if nargin < 4
    subsample = 'off';
end

if fp_info.Number > 0 
    
    if strcmp(subsample,'on')
        XI = a_data.Rate/data.Rate:a_data.Rate/data.Rate:length(a_data.COPx1);
        data.FP_Rate = data.Rate;
    else XI = 1:length(a_data.COPx1);
        data.FP_Rate = a_data.Rate;
    end      

    for i = 1:fp_info.Number
        % first calculate the orientation of the forceplate relative to the
        % global coordinate system. This is done from knowledge of the position
        % of the corners of the forceplate
        plate_centre = round(mean([fp_info.(['Corner1_' num2str(i)]).X fp_info.(['Corner2_' num2str(i)]).X ...
            fp_info.(['Corner3_' num2str(i)]).X fp_info.(['Corner4_' num2str(i)]).X;
            fp_info.(['Corner1_' num2str(i)]).Y fp_info.(['Corner2_' num2str(i)]).Y ...
            fp_info.(['Corner3_' num2str(i)]).Y fp_info.(['Corner4_' num2str(i)]).Y],2));
        alpha = atan2(round(fp_info.(['Corner1_' num2str(i)]).X) - round(fp_info.(['Corner4_' num2str(i)]).X),...
            round(fp_info.(['Corner1_' num2str(i)]).Y) - round(fp_info.(['Corner4_' num2str(i)]).Y));

        % now caluclate the COP position relative to the global coordinate
        % system
        data.ForcePlate(i).COP(:,1) = interp1(-a_data.(['COPx' num2str(i)])*cos(alpha)+... %need to flip the x-axis system to match that of the global coordinate system 
            a_data.(['COPy' num2str(i)])*sin(alpha)+plate_centre(1),XI);
        data.ForcePlate(i).COP(:,2) = interp1(a_data.(['COPy' num2str(i)])*cos(alpha)+...
            a_data.(['COPx' num2str(i)])*sin(alpha)+plate_centre(2),XI);
        data.ForcePlate(i).COP(:,3) = -interp1(a_data.(['COPz' num2str(i)]),XI)'; %need to flip the z-axis system to match that of the global coordinate system

        % now calculate the forces and moments about the global coordinate
        % system
        data.ForcePlate(i).F(:,1) = interp1(-a_data.(['Fx' num2str(i)])*cos(alpha) + ... %need to flip the x-axis system to match that of the global coordinate system 
            a_data.(['Fy' num2str(i)])*sin(alpha),XI);
        data.ForcePlate(i).F(:,2) = interp1(a_data.(['Fy' num2str(i)])*cos(alpha) - ...
            -a_data.(['Fx' num2str(i)])*sin(alpha),XI);
        data.ForcePlate(i).F(:,3) = -interp1(a_data.(['Fz' num2str(i)]),XI)'; %need to flip the z-axis system to match that of the global coordinate system

        % first calculate the moment about the COP in the FoP frame of
        % reference
        c = -fp_info.(['Origin_' num2str(i)])(3)/1000;

        % Mx and My should be 0 by definition
        Mx = (interp1(a_data.(['Mx' num2str(i)]),XI)'/1000 - (c*interp1(a_data.(['Fy' num2str(i)]),XI)') ... %...
                - (interp1(a_data.(['COPy' num2str(i)])/1000,XI)'.*interp1(a_data.(['Fz' num2str(i)]),XI)'));
        My = (interp1(a_data.(['My' num2str(i)]),XI)'/1000 + (c*interp1(a_data.(['Fx' num2str(i)]),XI)') ...
                + (interp1(a_data.(['COPx' num2str(i)])/1000,XI)'.*interp1(a_data.(['Fz' num2str(i)]),XI)')); %...
        % Mz is the free torque
        Mz = (interp1(a_data.(['Mz' num2str(i)]),XI)'/1000 ...
                + (interp1(a_data.(['COPy' num2str(i)])/1000,XI)'.*interp1(a_data.(['Fx' num2str(i)]),XI)') ...
                - (interp1(a_data.(['COPx' num2str(i)])/1000,XI)'.*interp1(a_data.(['Fy' num2str(i)]),XI)'));

        % Now rotate moments to the global coordinate system like the forces    
        data.ForcePlate(i).M(:,1) = -Mx*cos(alpha) + My*sin(alpha);
        data.ForcePlate(i).M(:,2) = My*cos(alpha) - (-Mx*sin(alpha));
        data.ForcePlate(i).M(:,3) = -Mz;
    end 
    
    else
        disp('There is no force plate data in this file');   
end
