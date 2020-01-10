function [data] = grf_assign_treadmill(data,markers,segments,FP_thresh,filt_cut)

if nargin < 5
    filt_cut = 20;
end

if nargin < 4
    FP_thresh = 25;
end

if nargin < 3
    segments = {'calcn_r';'calcn_l'};
end

if nargin < 2
    markers = {'RHEEL';'LHEEL'};
end

if nargin < 1
    load data.mat
end

% loop through both plates and determine which foot is contact with which
% plate at all instants in time
for f = 1:length(data.fp_data.GRF_data)
    clear a b c d on off Ind_on
    % define the vertical ground reaction force of plate
    F = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.fp_data.GRF_data(f).F(:,3));
    % dfine when a foot is on or off each plate - use a force threshold
    % defined aboved
    a = find(F>FP_thresh);
    b = find(flipud(F)>FP_thresh);
    on = [a(1); flipud(abs((b(diff(b)>15))-length(F)+1))];
    off = [a(diff(a)>15); (length(F)-b(1)+1)];
    if off(1) < on(1)
        off(1) = [];
    end
    
    % use the defined marker attached to foot to determine relationship
    % between this marker and the COP for the plate
    
    if strcmp(data.marker_data.Info.units.ALLMARKERS,'m')
        p_sc = 1000;
    else p_sc = 1;
    end
    
    % define right foot marker
    D.(markers{1}) = interp1(data.marker_data.Time,data.marker_data.Markers.(markers{1}),data.fp_data.Time,'linear','extrap')*p_sc;
    % define left foot marker
    D.(markers{2}) = interp1(data.marker_data.Time,data.marker_data.Markers.(markers{2}),data.fp_data.Time,'linear','extrap')*p_sc;
    
    % define COP - filter this to make calculation easier (noise becomes a
    % problem)
    COP = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.fp_data.GRF_data(f).P);
    
    % initilise some arrays - the left and right foot forces for this plate
    % and associate moment and COP as well as an indice which will detemine
    for i = 1:length(markers)
        data.GRF.(markers{i})(f).F = zeros(size(data.fp_data.GRF_data(f).F));
        data.GRF.(markers{i})(f).M = zeros(size(data.fp_data.GRF_data(f).M));
        data.GRF.(markers{i})(f).P = zeros(size(data.fp_data.GRF_data(f).P));
    end
    
    % go through each on and off event period and determine which marker is
    % closest to the COP of the force vector for that time. If it is
    % closest to right foot, then assign to right foot and vice versa.
    for i = 1:min([length(on) length(off)])
        clear c d
        % calcualte the distance between right and left markers and COP
        C1 = mean(dist_markers2D(D.(markers{1})(on(i):off(i),1:2),COP(on(i):off(i),1:2)));%left
        C2 = mean(dist_markers2D(D.(markers{2})(on(i):off(i),1:2),COP(on(i):off(i),1:2)));%right
        % if C1 is smaller it must be left and C2 smaller must be right        
        if C2 < C1
            if C2 < 300
                % store in a force structure for each foot and each plate (f)
                % will have its own force assigned to that foot.
                data.GRF.(markers{2})(f).F(on(i):off(i),:) = data.fp_data.GRF_data(f).F(on(i):off(i),:);
                data.GRF.(markers{2})(f).M(on(i):off(i),:) = data.fp_data.GRF_data(f).M(on(i):off(i),:);
                data.GRF.(markers{2})(f).P(on(i):off(i),:) = data.fp_data.GRF_data(f).P(on(i):off(i),:);
            end
        else
            if C1 < 300
                data.GRF.(markers{1})(f).F(on(i):off(i),:) = data.fp_data.GRF_data(f).F(on(i):off(i),:);
                data.GRF.(markers{1})(f).M(on(i):off(i),:) = data.fp_data.GRF_data(f).M(on(i):off(i),:);
                data.GRF.(markers{1})(f).P(on(i):off(i),:) = data.fp_data.GRF_data(f).P(on(i):off(i),:);
            end
        end
        
    end
    
    % subplot(2,1,1), plot([data.GRF.FR(f).F(:,3) data.GRF.FL(f).F(:,3)]); hold on;
   
end

% the summed (or total) force vector for each foot is just the sum of each plate  
data.GRF.(segments{1}).F = data.GRF.(markers{1})(1).F + data.GRF.(markers{1})(2).F;
data.GRF.(segments{2}).F = data.GRF.(markers{2})(1).F + data.GRF.(markers{2})(2).F;

% the summed (or total) free moment for each foot is just the sum of each plate  
data.GRF.(segments{1}).M = data.GRF.(markers{1})(1).M + data.GRF.(markers{1})(2).M;
data.GRF.(segments{2}).M = data.GRF.(markers{2})(1).M + data.GRF.(markers{2})(2).M;

% calcualte the centre of pressure as the weighted avaerage between the COP
% in the X and Y direction - Z = 0. 
data.GRF.(segments{1}).P(:,1) = (data.GRF.(markers{1})(1).P(:,1).*(data.GRF.(markers{1})(1).F(:,3)./data.GRF.(segments{1}).F(:,3)))+(data.GRF.(markers{1})(2).P(:,1).*(data.GRF.(markers{1})(2).F(:,3)./data.GRF.(segments{1}).F(:,3)));
data.GRF.(segments{1}).P(:,2) = (data.GRF.(markers{1})(1).P(:,2).*(data.GRF.(markers{1})(1).F(:,3)./data.GRF.(segments{1}).F(:,3)))+(data.GRF.(markers{1})(2).P(:,2).*(data.GRF.(markers{1})(2).F(:,3)./data.GRF.(segments{1}).F(:,3)));
data.GRF.(segments{1}).P(:,3) = zeros(size(data.GRF.(segments{1}).P(:,1)));
% determine when the COP is NaN (0/0) which is when foot is not on plate
Rnan = isnan(data.GRF.(segments{1}).P(:,1));

data.GRF.(segments{2}).P(:,1) = (data.GRF.(markers{2})(1).P(:,1).*(data.GRF.(markers{2})(1).F(:,3)./data.GRF.(segments{2}).F(:,3)))+(data.GRF.(markers{2})(2).P(:,1).*(data.GRF.(markers{2})(2).F(:,3)./data.GRF.(segments{2}).F(:,3)));
data.GRF.(segments{2}).P(:,2) = (data.GRF.(markers{2})(1).P(:,2).*(data.GRF.(markers{2})(1).F(:,3)./data.GRF.(segments{2}).F(:,3)))+(data.GRF.(markers{2})(2).P(:,2).*(data.GRF.(markers{2})(2).F(:,3)./data.GRF.(segments{2}).F(:,3)));
data.GRF.(segments{2}).P(:,3) = zeros(size(data.GRF.(segments{2}).P(:,1)));
% determine when the COP is NaN (0/0) which is when foot is not on plate
Lnan = isnan(data.GRF.(segments{2}).P(:,1));

% calculate the free moments
%data.GRF.(segments{1}).M = zeros(size(data.GRF.(segments{1}).P));
%data.GRF.(segments{1}).M(:,3) = data.GRF.FR(1).M(:,3)+data.GRF.FR(2).M(:,3)-data.GRF.FR(1).F(:,1).*data.

% filter the total (summed) force and moment data
data.GRF.(segments{1}).F = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.GRF.(segments{1}).F);
data.GRF.(segments{2}).F = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.GRF.(segments{2}).F);

data.GRF.(segments{1}).M = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.GRF.(segments{1}).M);
data.GRF.(segments{2}).M = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,data.GRF.(segments{2}).M);

% filtering makes some zero values for vertical force, make this zero
data.GRF.(segments{1}).F(data.GRF.(segments{1}).F(:,3)<0,3) = 0;
data.GRF.(segments{2}).F(data.GRF.(segments{2}).F(:,3)<0,3) = 0;

% filter the total COP - interpolate NaNs before filtering and then
% reassign NaNs to a value of zero
data.GRF.(segments{1}).P = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,interp1(find(~Rnan),data.GRF.(segments{1}).P(~Rnan,:),(1:length(data.GRF.(segments{1}).P(:,1)))','linear','extrap'));
data.GRF.(segments{1}).P(Rnan,:) = 0;
data.GRF.(segments{2}).P = matfiltfilt(1/data.fp_data.Info(1).frequency,filt_cut,2,interp1(find(~Lnan),data.GRF.(segments{2}).P(~Lnan,:),(1:length(data.GRF.(segments{2}).P(:,1)))','linear','extrap'));
data.GRF.(segments{2}).P(Lnan,:) = 0;

%figure(1); plot([data.GRF.(segments{1}).F(:,3) data.GRF.(segments{2}).F(:,3)]); pause;



