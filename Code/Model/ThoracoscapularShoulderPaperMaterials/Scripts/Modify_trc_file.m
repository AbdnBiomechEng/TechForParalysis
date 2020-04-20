[filename, filepath] = uigetfile('*.trc','Pick the trial.trc file.');
TRC_filename = [filepath,filename];
trcData = dlmread(TRC_filename, '\t', 6, 0);
trcData_noTimeFrame=trcData(:,3:end);
time = trcData(:,2);
nframe = trcData(:,1);

opts = delimitedTextImportOptions("NumVariables", 189);
% Specify range and delimiter
opts.DataLines = [1, 4];
opts.Delimiter = "\t";
HeaderTRC = readtable(TRC_filename, opts);
HeaderTRC = table2array(HeaderTRC);
framerate=str2double(HeaderTRC(3,1));
startFrame=str2double(HeaderTRC(3,7));
noMarkers=str2double(HeaderTRC(3,4));
% markers=HeaderTRC(4,:);
% markers(cellfun('isempty',markers)) = [] ; % deleting empty cells

if noMarkers==56
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

AddedMarkers={'ac', 'ijc', 'ij_proj_c7', 'ij_proj_px', 'ij_proj_ac',...
    'ai_proj_aa_x', 'ai_proj_aa_z', 'ai_proj_ts'};
markers=[markers AddedMarkers];

% Get initial position of IJ to translate all data so that at the initial
% frame, IJ = (0,0,0)
IJ_init =trcData_noTimeFrame(1,154:156);
IJ_init_mat = repmat(IJ_init,length(time),size(trcData_noTimeFrame,2)/3);

% all data are balanced relative to IJ with IJ at initial frame being
% (0,0,0) and divided by 1000 to convert from mm to metres.
trcData_balanced=(trcData_noTimeFrame-IJ_init_mat)/1000;
ac=trcData_balanced(:,142:144);
ijc=trcData_balanced(:,154:156);
ij_proj_c7=[trcData_balanced(:,148) trcData_balanced(:,155:156)]; %x-c7, the rest ij
ij_proj_px=[trcData_balanced(:,154) trcData_balanced(:,158) trcData_balanced(:,156)]; %y-px, the rest ij
ij_proj_ac=[trcData_balanced(:,154) trcData_balanced(:,155) trcData_balanced(:,144)]; %z-ac, the rest ij
ai_proj_aa_x=[trcData_balanced(:,130) trcData_balanced(:,134) trcData_balanced(:,135)]; %x-aa, the rest ai
ai_proj_ts=[trcData_balanced(:,133) trcData_balanced(:,137) trcData_balanced(:,135)]; %y-ts, the rest ai
ai_proj_aa_z=[trcData_balanced(:,133) trcData_balanced(:,134) trcData_balanced(:,132)]; %z-aa, the rest ai

elseif noMarkers==62
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

AddedMarkers={'ac', 'ijc', 'ij_proj_c7', 'ij_proj_px', 'ij_proj_ac',...
    'ai_proj_aa_x', 'ai_proj_aa_z', 'ai_proj_ts'};
markers=[markers AddedMarkers];

% Get initial position of IJ to translate all data so that at the initial
% frame, IJ = (0,0,0)
IJ_init =trcData_noTimeFrame(1,172:174);
IJ_init_mat = repmat(IJ_init,length(time),size(trcData_noTimeFrame,2)/3);

% all data are balanced relative to IJ with IJ at initial frame being
% (0,0,0) and divided by 1000 to convert from mm to metres.
trcData_balanced=(trcData_noTimeFrame-IJ_init_mat)/1000;
ac=trcData_balanced(:,160:162);
ijc=trcData_balanced(:,172:174);
ij_proj_c7=[trcData_balanced(:,166) trcData_balanced(:,173:174)]; %x-c7, the rest ij
ij_proj_px=[trcData_balanced(:,172) trcData_balanced(:,175) trcData_balanced(:,174)]; %y-px, the rest ij
ij_proj_ac=[trcData_balanced(:,172) trcData_balanced(:,173) trcData_balanced(:,162)]; %z-ac, the rest ij
ai_proj_aa_x=[trcData_balanced(:,148) trcData_balanced(:,152) trcData_balanced(:,153)]; %x-aa, the rest ai
ai_proj_ts=[trcData_balanced(:,151) trcData_balanced(:,155) trcData_balanced(:,153)]; %y-ts, the rest ai
ai_proj_aa_z=[trcData_balanced(:,151) trcData_balanced(:,152) trcData_balanced(:,150)]; %z-aa, the rest ai

else
    error('More markers than expected');
end

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = ['Frame#\tTime\t'];
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe time];

data_out = [data_out trcData_balanced ac ijc ij_proj_c7 ij_proj_px...
    ij_proj_ac ai_proj_aa_x ai_proj_ts ai_proj_aa_z]';

%data_out = [data_out trcData_balanced]';
 
for imark = 1:length(markers)
    % now loop through each maker name and make marker name with 3 tabs for the
    % first line and the X Y Z columns with the marker numnber on the second
    % line all separated by tab delimeters
    dataheader1 = [dataheader1 markers{imark} '\t\t\t'];
    dataheader2 = [dataheader2 'X' num2str(imark) '\t' 'Y' num2str(imark) '\t'...
        'Z' num2str(imark) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
end

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
