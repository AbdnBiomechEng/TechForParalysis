function data = c3d2trc(varargin)
% function c3d2trc(file) OR
% function c3d2trc([data, a_data, fp_info, MP_info])
%
% Function to convert data from a C3D file into the TRC and MOT file
% formats for OpenSim
%
% INPUT -   file - the C3D file path that you wish to load (leave blank to
%               choose from a dialog box) OR
%           (data, a_data, fp_info, MP_info, filename) - data from previously loaded
%           C3D file using loadc3dfile.m as well as a filename string
%
% OUTPUT -  data - (optional) structure containing the relevant data from
%                 the c3dfile
%           Creates the TRC file and _grf.MOT file for OpenSim
%
% Written by Glen Lichtwark (University of Queensland) & Rod Barrett
% (Griffith University)

%% load data
if length(varargin) < 2
% load C3d file
    if nargin == 1
        file = varargin{1};
        if isempty(fileparts(file))
            pname = cd;
            pname = [pname '\'];
            fname = file;
        else [pname, name, ext] = fileparts(file);
            fname = [name ext];
        end
    else [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
    end

    [data, a_data, fp_info, MP_info] = loadc3dfile([pname, fname]);

    elseif length(varargin) == 5
        data = varargin{1};
        a_data = varargin{2};
        fp_info = varargin{3};
        MP_info = varargin{4};
        fname = varargin{5};
        pname = cd;
        pname = [pname '\'];
    else
        error('Wrong number of inputs... needs to be either filename or 4 structure arrays and filename string ([data, a_data, fp_info, MP_info filename])');
end

%%
% sort the C3D file so we know what is Marker data and what is calculated
% data --> if this is already sorted it won't matter
if ~isfield(data,'Markers')
data = sort_c3d(data);
end

data.C3D_Filename = [pname fname]; 

%%
% try and find the mass and height of the subject from the MP file and
% store into the data structure
if isfield(MP_info,'Bodymass')
    data.Mass = MP_info.Bodymass;
else data.Mass = 75; % default to 75kg
end

if isfield(MP_info,'Height')
    data.Height = MP_info.Height;
else data.Height = 1750; % default to 1750mm
end

%%
% calculate the force data in the global coordinate system
data = c3d2grf(data,a_data,fp_info,'off');

%%
% define some parameters 
nrows = data.End_Frame - data.Start_Frame + 1;
nmarkers = length(fieldnames(data.Markers));

nframe = 1:nrows;
time = (nframe/data.Rate)-1/data.Rate;
data.time = time;

%%
% animate the data - there are some custom models which have their own
% animation drawings otherwise just plot marker positions
% if strcmp(data.Model,'PlugInGait')
% animate_PiG(data,a_data,fp_info,5);
% elseif strcmp(data.Model,'OpenSimFullBodyModel') || strcmp(data.Model,'OpenSimModel')
%     animate_UWA(data,a_data,fp_info,5);
% else animate_markers(data,a_data,fp_info,5);
% end
animate_markers(data,a_data,fp_info,5);
%%
% we need to reorder the lab coordinate system to match that of the OpenSim
% system --> SKIP THIS STEP IF LAB COORDINATE SYSTEM IS SAME AS MODEL
% SYSTEM
markers = fieldnames(data.Markers); % get markers names
% go through each marker field and re-order from X Y Z to Y Z X
for i = 1:nmarkers
   data.Markers.(markers{i}) =  [data.Markers.(markers{i})(:,2)...
       data.Markers.(markers{i})(:,3) data.Markers.(markers{i})(:,1)];
end

%%
% now we need to make the headers for the column headings for the TRC file
% which are made up of the marker names and the XYZ for each marker

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = 'Frame#\tTime\t';
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe; time];

% now loop through each maker name and make marker name with 3 tabs for the
% first line and the X Y Z columns with the marker numnber on the second
% line all separated by tab delimeters
% each of the data columns (3 per marker) will be in floating format with a
% tab delimiter - also add to the data matrix
for i = 1:nmarkers
    dataheader1 = [dataheader1 markers{i} '\t\t\t'];    
    dataheader2 = [dataheader2 'X' num2str(i) '\t' 'Y' num2str(i) '\t'...
        'Z' num2str(i) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
    % add 3 rows of data for the X Y Z coordinates of the current marker
    % first check for NaN's and fill with a linear interpolant - warn the
    % user of the gaps
    clear m
    m = find(isnan(data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),1))>0);
    if ~isempty(m);
        clear t d
        disp(['Warning -' markers{i} ' data missing in parts. Frames ' num2str(m(1)) '-'  num2str(m(end))])
        t = time;
        t(m) = [];
        d = data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:);
        d(m,:) = [];
        data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:) = interp1(t,d,time,'linear','extrap');
    end
    data_out = [data_out; data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:)'];
end
dataheader1 = [dataheader1 '\n'];
dataheader2 = [dataheader2 '\n'];
format_text = [format_text '\n'];

disp('Writing trc file...') 

%Output marker data to an OpenSim TRC file

newfilename = strrep(fname,'c3d','trc');

data.TRC_Filename = [pname newfilename];

%open the file
fid_1 = fopen([pname newfilename],'w');

% first write the header data
fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n',newfilename);
fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', data.Rate, data.Rate, nrows, nmarkers, data.units, data.Rate,data.Start_Frame,data.End_Frame); 
fprintf(fid_1, dataheader1);
fprintf(fid_1, dataheader2);

% then write the output marker data
fprintf(fid_1, format_text,data_out);

% close the file
fclose(fid_1);

disp('Done.')

%%
% Write motion file containing GRFs

disp('Writing grf.mot file...')

if fp_info.Number > 0 
F = data.FP_Rate/data.Rate;

fp_time = 0:1/data.FP_Rate:((F*data.End_Frame)-(F*data.Start_Frame))/data.FP_Rate;

% initialise force data matrix with the time array and column header
force_data_out = fp_time';
force_header = 'time\t';
force_format = '%20.6f\t';

% go through each marker field and re-order from X Y Z to Y Z X and place
% into data array and add data to the force data matrix --> also need to
% divide by 1000 to convert to m from mm

% if there are more than two force plates then we will also need to combine
% into two data sets, one for left and one for right (or one if only one
% forceplate is hit)

if isfield(data.Markers,'RHEEL')
    RHEE = 'RHEEL';
    LHEE = 'LHEEL';
else RHEE = 'RHEE';
    LHEE = 'LHEE';
end

if isfield(data.Markers,'RCAL')
    RHEE = 'RCAL';
    LHEE = 'LCAL';
end

for i = 1:fp_info.Number
    % reoder data so lab coordinate system to match that of the OpenSim
    % system
   data.ForcePlate(i).COP =  [data.ForcePlate(i).COP(:,2)/1000 ...
       data.ForcePlate(i).COP(:,3)/1000 data.ForcePlate(i).COP(:,1)/1000];
   data.ForcePlate(i).F =  [data.ForcePlate(i).F(:,2) ...
       data.ForcePlate(i).F(:,3) data.ForcePlate(i).F(:,1)];
   data.ForcePlate(i).M =  [data.ForcePlate(i).M(:,2) ...
       data.ForcePlate(i).M(:,3) data.ForcePlate(i).M(:,1)];
   
   % do some cleaning of the COP before and after contact
   b = find(abs(diff(data.ForcePlate(i).COP(:,3)))>0);
   if ~isempty(b)
       for j = 1:3
            data.ForcePlate(i).COP(1:b(1),j) = data.ForcePlate(i).COP(b(1)+1,j);
            data.ForcePlate(i).COP(b(end):end,j) = data.ForcePlate(i).COP(b(end)-1,j);
       end
   end
   
   % define the period which we are analysing
   K = (F*data.Start_Frame):1:(F*data.End_Frame);
   
   % add the force, COP and moment data for current plate to the force matrix 
   force_data_out = [force_data_out data.ForcePlate(i).F(K,:) data.ForcePlate(i).COP(K,:) data.ForcePlate(i).M(K,:)];
   % define the header and formats
   force_header = [force_header num2str(i) '_ground_force_vx\t' num2str(i) '_ground_force_vy\t' num2str(i) '_ground_force_vz\t'...
       num2str(i) '_ground_force_px\t' num2str(i) '_ground_force_py\t' num2str(i) '_ground_force_pz\t' ...
       num2str(i) '_ground_torque_x\t' num2str(i) '_ground_torque_y\t' num2str(i) '_ground_torque_z\t'];
   force_format = [force_format '%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t'];
   
   % find which foot is on plate (if any) and store to the ForcePlate
   % structure object. Always make right foot strikes the first set of
   % columns in the MOT file and the left foot strikes the second set
   clear a
   a = find(data.ForcePlate(i).F(F*data.Start_Frame:F*data.End_Frame,2) > 0.01*max(data.ForcePlate(i).F(F*data.Start_Frame:F*data.End_Frame,2)));
   aa = a + (F*data.Start_Frame)-1;
   if ~isempty(a)
        if nanmean(dist_markers(data.ForcePlate(i).COP(aa,:),[data.Markers.(RHEE)(round(aa/F),1) zeros(size(data.Markers.(RHEE)(round(aa/F),1))) data.Markers.(RHEE)(round(aa/F),3)]/1000)) < ...
            nanmean(dist_markers(data.ForcePlate(i).COP(aa,:),[data.Markers.(LHEE)(round(aa/F),1) zeros(size(data.Markers.(LHEE)(round(aa/F),1))) data.Markers.(LHEE)(round(aa/F),3)]/1000))
            data.ForcePlate(i).Leg = {'Right'};         
        else data.ForcePlate(i).Leg = {'Left'};
        end
   else data.ForcePlate(i).Leg = {'None'};
   end
end

force_header = [force_header(1:end-2) '\n'];
force_format = [force_format(1:end-2) '\n'];

% assign a value of zero to any NaNs
force_data_out(logical(isnan(force_data_out))) = 0;

newfilename = [fname(1:end-4) '_grf.mot'];

data.GRF_Filename = [pname newfilename];

fid_2 = fopen([pname newfilename],'w');

% write the header information
fprintf(fid_2,'name %s\n',newfilename);
fprintf(fid_2,'datacolumns %d\n', size(force_data_out,2));  % total # of datacolumns
fprintf(fid_2,'datarows %d\n',length(fp_time)); % number of datarows
fprintf(fid_2,'range %f %f\n',time(1),time(end)); % range of time data
fprintf(fid_2,'endheader\n');
fprintf(fid_2,force_header);

% write the data
fprintf(fid_2,force_format,force_data_out');

fclose(fid_2);

disp('Done.')
else disp('No force plate information available.')
end
