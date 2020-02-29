function data = assign_forces(data,assign_markers,assign_bodies, thresh)
% function data = assign_forces(data,assign_markers,assign_bodies, thresh)
%
% Function to assign any recorded forces to a specific body based on the
% position of a nominated marker attached to the body.
%
% INPUT -   data - structure containing fields from from previously loaded
%               C3D file using btk_loadc3d.m as well as a filename string
%           assign_markers - cell array of marker names to be used as
%               guides to match to a force vector COP (e.g. heel marker)
%           assign_bodies - cell array with the matching body name that any
%               matching forces can be assigned to
%           thresh - the treshold mean distance from the marker to the COP
%               that is used to assess a positive assignment (in meters - 
%               defaults to 0.2m)
%
% OUTPUT -  data - structure containing the relevant data with assignments           
%
% Written by Glen Lichtwark (University of Queensland)
% Updated September 2012

if nargin<4
    thresh = 0.2;
end

if ~isfield(data,'marker_data')
    error(['Ensure that the data structure contains a field called "marker_data" '...
    'that contains marker fields and coordinates - see btk_loadc3d']);
end

if ~isfield(data,'fp_data')
    error(['Ensure that the data structure contains a field called "fp_data" '...
    'that contains force plate data - see btk_loadc3d']);
end

if iscell(assign_markers) 
    if iscell(assign_bodies)
        if length(assign_markers) ~= length(assign_bodies)
            error('Cell arrays for assigned markers and bodies must be the same length')
        else N = length(assign_markers);
        end
    else error('Assigned marker list and assigned bodies must both be cell arrays of same length with paired matchings')
    end
else error('The body and marker assignment lists must be cell arrays');
end

% examine the force signals for each forcelpate and determine whether any
% of the markers are close to the COP to indicate that this force should be
% assigned to the bodies that the marker attaches to

% define the ratio of force sampling frequency to marker sampling frequency
F = data.fp_data.Info(1).frequency/data.marker_data.Info.frequency; %assume same sampling frequency on all plates!!!

% determine if units are in mm or m
if strcmp(data.marker_data.Info.units.ALLMARKERS,'mm')
    p_sc = 1000;
else p_sc = 1;
end

for i = 1:length(data.fp_data.GRF_data)
    
    data.AssignForce.ExForce{i} = ['ExternalForce_' num2str(i)];
    
    a = find(data.fp_data.GRF_data(i).F(:,2) > 15);
    aa = round(a(1)/F:a(end)/F);
   
    if ~isempty(a)
        D = zeros(1:N);
        for j = 1:N
            D(j) = nanmean(dist_markers(data.fp_data.GRF_data(i).P(aa*F,:),...
                [data.marker_data.Markers.(assign_markers{j})(aa,1) ...
                zeros(size(data.marker_data.Markers.(assign_markers{j})(aa,1))) ...
                data.marker_data.Markers.(assign_markers{j})(aa,3)]/p_sc));
        end
        [~,I] = min(D);
        if D(I)<thresh
            data.AssignForce.ApBodies{i} = assign_bodies{I};
        else data.AssignForce.ApBodies{i} = 'ground'; 
        end
    else data.AssignForce.ApBodies{i} = 'ground';       
    end
    
end

