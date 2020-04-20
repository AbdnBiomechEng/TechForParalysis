% Modify_trc_file reads a .trc file that was previously created in a gait
% lab and adds new marker names and coordinates. The new trc file can then
% used in an OpenSim model OSM.osim
% 
% Be careful with order of marker names. The names should be in the same
% order as the original trc file. The coordinate orders should match the
% arrangement of the marker names.
%
% Rosti Readioff, April 2020

% call the trc file that you want to modify.
[filename, filepath] = uigetfile('*.trc','Pick the trial.trc file.');
TRC_filename = [filepath,filename];
trcData = dlmread(TRC_filename, '\t', 6, 0); % read marker coordinates.
trcData_noTimeFrame=trcData(:,3:end); % all coordinates without the frame & time columns.
time = trcData(:,2); % time 
nframe = trcData(:,1); % frame

% Reading the text in the first four lines in the trc file.
opts = delimitedTextImportOptions("NumVariables", 189); 
% Specify range and delimiter
opts.DataLines = [1, 4];
opts.Delimiter = "\t";
HeaderTRC = readtable(TRC_filename, opts);
HeaderTRC = table2array(HeaderTRC); 
framerate=str2double(HeaderTRC(3,1)); % framerate
startFrame=str2double(HeaderTRC(3,7)); % start of frame
noMarkers=str2double(HeaderTRC(3,4)); % no. of markers in the original trc

% markers=HeaderTRC(4,:);
% markers(cellfun('isempty',markers)) = [] ; % deleting empty cells

% an if condition to check for the two different types of trc files
% generated in the lab, which are loaded trials and no loaded trials.
if noMarkers==56 % This is for no loaded trials.
    markers={'ORLAU Subject 2:ST1','ORLAU Subject 2:ST2',...
        'ORLAU Subject 2:ST3','ORLAU Subject 2:RACR1',...
        'ORLAU Subject 2:RACR2','ORLAU Subject 2:RACR3',...
        'ORLAU Subject 2:RUA1','ORLAU Subject 2:RUA2',...
        'ORLAU Subject 2:RUA3','ORLAU Subject 2:RUA4',...
        'ORLAU Subject 2:RFA1','ORLAU Subject 2:RFA2',...
        'ORLAU Subject 2:RFA3','ORLAU Subject 2:RFA4',...
        'ORLAU Subject 2:RH1','ORLAU Subject 2:RH2',...
        'ORLAU Subject 2:RH3','RACR4','RUA1','RUA2','RUA3','RUA4',...
        'RFA1','RFA2','RFA3','RFA4','RH1','RH2','RH3','RACR1','RACR2',...
        'RACR3','ST1','ST2','ST3','EpL','EpM','rs','us','RMC3',...
        'RMCP2','RMCP3','RMCP5','aa','ai','ts','pc','acc','gu','c7',...
        't8','ij','px','centelbow_humerus','centelbow_ulna','centusrs'};

% Below are the new markers added to the trial - please change this if you
% wanted to add more markers to the trial.
AddedMarkers={'ac', 'ijc', 'ij_proj_c7', 'ij_proj_px', 'ij_proj_ac',...
    'ai_proj_aa_x', 'ai_proj_aa_z', 'ai_proj_ts'};
markers=[markers AddedMarkers]; % combines old and new markers.

% Get initial position of IJ to translate all data so that at the initial
% frame, IJ = (0,0,0).
IJ_init =trcData_noTimeFrame(1,154:156);
IJ_init_mat = repmat(IJ_init,length(time),size(trcData_noTimeFrame,2)/3);

% all data are balanced relative to IJ with IJ at initial frame being
% (0,0,0) and divided by 1000 to convert from mm to metres.
trcData_balanced=(trcData_noTimeFrame-IJ_init_mat)/1000;

% coordinates for the new added markers - please change this if you
% added more markers to the trial.
ac=trcData_balanced(:,142:144);
ijc=trcData_balanced(:,154:156);
ij_proj_c7=[trcData_balanced(:,148) trcData_balanced(:,155:156)]; %x-c7, the rest ij
ij_proj_px=[trcData_balanced(:,154) trcData_balanced(:,158) trcData_balanced(:,156)]; %y-px, the rest ij
ij_proj_ac=[trcData_balanced(:,154) trcData_balanced(:,155) trcData_balanced(:,144)]; %z-ac, the rest ij
ai_proj_aa_x=[trcData_balanced(:,130) trcData_balanced(:,134) trcData_balanced(:,135)]; %x-aa, the rest ai
ai_proj_ts=[trcData_balanced(:,133) trcData_balanced(:,137) trcData_balanced(:,135)]; %y-ts, the rest ai
ai_proj_aa_z=[trcData_balanced(:,133) trcData_balanced(:,134) trcData_balanced(:,132)]; %z-aa, the rest ai

elseif noMarkers==62 % This is for loaded trials.
markers={'ORLAU Subject 2:ST1', 'ORLAU Subject 2:ST2',...
        'ORLAU Subject 2:ST3', 'ORLAU Subject 2:RACR1',...
        'ORLAU Subject 2:RACR2', 'ORLAU Subject 2:RACR3',...
        'ORLAU Subject 2:RUA1', 'ORLAU Subject 2:RUA2',...
        'ORLAU Subject 2:RUA3', 'ORLAU Subject 2:RUA4',...
        'ORLAU Subject 2:RFA1', 'ORLAU Subject 2:RFA2',...
        'ORLAU Subject 2:RFA3',	'ORLAU Subject 2:RFA4',...
        'ORLAU Subject 2:RH1', 'ORLAU Subject 2:RH2','ORLAU Subject 2:RH3',...
        'L1', 'L2', 'RACR4', 'L4', 'L5', 'L3', 'L6', 'RUA1', 'RUA2',...
        'RUA3', 'RUA4', 'RFA1', 'RFA2', 'RFA3', 'RFA4', 'RH1', 'RH2', ...
        'RH3','RACR1', 'RACR2', 'RACR3', 'ST1', 'ST2', 'ST3', 'EpL',...
        'EpM', 'rs', 'us', 'RMC3', 'RMCP2', 'RMCP3', 'RMCP5', 'aa',...
        'ai', 'ts', 'pc', 'acc', 'gu', 'c7', 't8', 'ij', 'px',...
        'centelbow_humerus', 'centelbow_ulna', 'centusrs'};
    
% Below are the new markers added to the trial - please change this if you
% wanted to add more markers to the trial.
AddedMarkers={'ac', 'ijc', 'ij_proj_c7', 'ij_proj_px', 'ij_proj_ac',...
    'ai_proj_aa_x', 'ai_proj_aa_z', 'ai_proj_ts'};
markers=[markers AddedMarkers]; % combines old and new markers.

% Get initial position of IJ to translate all data so that at the initial
% frame, IJ = (0,0,0).
IJ_init =trcData_noTimeFrame(1,172:174);
IJ_init_mat = repmat(IJ_init,length(time),size(trcData_noTimeFrame,2)/3);

% all data are balanced relative to IJ with IJ at initial frame being
% (0,0,0) and divided by 1000 to convert from mm to metres.
trcData_balanced=(trcData_noTimeFrame-IJ_init_mat)/1000;

% coordinates for the new added markers - please change this if you
% added more markers to the trial.
ac=trcData_balanced(:,160:162);
ijc=trcData_balanced(:,172:174);
ij_proj_c7=[trcData_balanced(:,166) trcData_balanced(:,173:174)]; %x-c7, the rest ij
ij_proj_px=[trcData_balanced(:,172) trcData_balanced(:,175) trcData_balanced(:,174)]; %y-px, the rest ij
ij_proj_ac=[trcData_balanced(:,172) trcData_balanced(:,173) trcData_balanced(:,162)]; %z-ac, the rest ij
ai_proj_aa_x=[trcData_balanced(:,148) trcData_balanced(:,152) trcData_balanced(:,153)]; %x-aa, the rest ai
ai_proj_ts=[trcData_balanced(:,151) trcData_balanced(:,155) trcData_balanced(:,153)]; %y-ts, the rest ai
ai_proj_aa_z=[trcData_balanced(:,151) trcData_balanced(:,152) trcData_balanced(:,150)]; %z-aa, the rest ai

else % if number of markers doesn't match any of the two then you get an error.
    error('More markers than expected');
end

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = ['Frame#\tTime\t'];
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';

% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe time];

% combine above with all the original xyz coordinates and coordinates newly
% generated for the new markers.
data_out = [data_out trcData_balanced ac ijc ij_proj_c7 ij_proj_px...
    ij_proj_ac ai_proj_aa_x ai_proj_ts ai_proj_aa_z]';

% Loop through each maker name and make marker name with 3 tabs for the
% first line and the X Y Z columns with the marker numnber on the second
% line all separated by tab delimeters.
for imark = 1:length(markers)
    dataheader1 = [dataheader1 markers{imark} '\t\t\t'];
    dataheader2 = [dataheader2 'X' num2str(imark) '\t' 'Y' num2str(imark) '\t'...
        'Z' num2str(imark) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
end

% Refine data headers and text format
dataheader1 = [dataheader1 '\n'];
dataheader2 = [dataheader2 '\n'];
format_text = [format_text '\n'];

%open the file
fid_1 = fopen(TRC_filename,'w');

% first write the header data
fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n',filename);
fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', framerate, framerate, length(time), length(markers), 'm', framerate, startFrame, length(time));
fprintf(fid_1, dataheader1);
fprintf(fid_1, dataheader2);

% then write the output marker data
fprintf(fid_1, format_text,data_out);

% close the file
fclose(fid_1);

disp(['Written trc file ' TRC_filename]);
