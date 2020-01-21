function build_trc_file(filename,path,markers)
% Build_trc_file creates a .trc file that can be used as input to the
% OpenSim model OSM.osim
% For a different version run through Nexus, use Build_trc_file_Vicon.m
%
% Input
%     filename = the complete path and name if the .c3d file that we want
%                to analyze. The output .trc file has the same path and name.
%     markers (optional)  = cell array containing the names of markers
%                           to be included in the .trc file. If missing, all markers in
%                           OSM.osim are used (see line 16 below)
%
% Dimitra Blana, April 2016

[pathstr,name,~] = fileparts(filename);


if nargin<3,
    markers = {'c7','t8','ij','px','aa','ac','ai','ts','sc','gu','EpL','EpM','US','RS'};
end

if nargin<2
    path = pathstr;
end

if ~iscell(markers)
    error('Marker names should be in a cell array');
end

TRC_filename = fullfile(path,[name '.trc']);
data = btk_loadc3d(filename);

framerate = data.marker_data.Info.frequency;
startFrame = data.marker_data.First_Frame;
endFrame = data.marker_data.Last_Frame;

time = data.marker_data.Time;
nframe = double(startFrame):double(endFrame);

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = 'Frame#\tTime\t';
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe; time'];

% Get initial position of IJ to translate all data so that at the initial
% frame, IJ = (0,0,0)
IJ_init = data.marker_data.Markers.ij(1,:);
IJ_init_mat = repmat(IJ_init,length(time),1);


for imark = 1:length(markers)

    % now loop through each maker name and make marker name with 3 tabs for the
    % first line and the X Y Z columns with the marker numnber on the second
    % line all separated by tab delimeters
    dataheader1 = [dataheader1 markers{imark} '\t\t\t'];
    dataheader2 = [dataheader2 'X' num2str(imark) '\t' 'Y' num2str(imark) '\t'...
        'Z' num2str(imark) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];

    if (isfield(data.marker_data.Markers,markers{imark}))
        % get the input data
        XYZ = data.marker_data.Markers.(markers{imark}) - IJ_init_mat;
        % add 3 rows of data for the X Y Z coordinates of the current marker
        % DB: Vicon coordinate frame is x left, y back, z up. Convert to x
        % right, y up, z back for OpenSim model
        % RR: Vicon coordinate frame is x back, y left, z up. Convert to x
        % left, y up, z back for OpenSim model
        %data_out = [data_out; -XYZ(:,1)'/1000; XYZ(:,3)'/1000; XYZ(:,2)'/1000];
        data_out = [data_out; XYZ(:,2)'/1000; XYZ(:,3)'/1000; -XYZ(:,1)'/1000];        
        %data_out = [data_out; AxelRot(XYZ(:,2)'/1000, 90,[1 0 0],[]); AxelRot(XYZ(:,3)'/1000,0,[0 1 0],[]); AxelRot(-XYZ(:,1)'/1000,180,[0 0 1],[])];

    elseif strcmp(markers{imark},'GH') && isfield(data.marker_data.Markers,'GHhum')
        % GH in some cases is saved as GHhum
        % get the input data
        XYZ = data.marker_data.Markers.GHhum - IJ_init_mat;
        % add 3 rows of data for the X Y Z coordinates of the current marker
        % DB: Vicon coordinate frame is x left, y back, z up. Convert to x
        % right, y up, z back for OpenSim model
        % RR: Vicon coordinate frame is x back, y left, z up. Convert to x
        % left, y up, z back for OpenSim model
        %data_out = [data_out; -XYZ(:,1)'/1000; XYZ(:,3)'/1000; XYZ(:,2)'/1000];
        data_out = [data_out; XYZ(:,2)'/1000; XYZ(:,3)'/1000; -XYZ(:,1)'/1000];
        %data_out = [data_out; AxelRot(XYZ(:,2)'/1000, 90,[1 0 0],[]); AxelRot(XYZ(:,3)'/1000,0,[0 1 0],[]); AxelRot(-XYZ(:,1)'/1000,180,[0 0 1],[])];

    else
        error(['Invalid marker name ', markers{imark}]);
    end
end

dataheader1 = [dataheader1 '\n'];
dataheader2 = [dataheader2 '\n'];
format_text = [format_text '\n'];

%open the file
fid_1 = fopen(TRC_filename,'w');

% first write the header data
fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n',[name '.trc']);
fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', framerate, framerate, length(time), length(markers), 'm', framerate, startFrame, length(time));
fprintf(fid_1, dataheader1);
fprintf(fid_1, dataheader2);

% then write the output marker data
fprintf(fid_1, format_text,data_out);

% close the file
fclose(fid_1);

disp(['Written trc file ' TRC_filename]);
