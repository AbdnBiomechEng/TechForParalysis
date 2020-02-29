function fileout = mot_grf(ik_mot_file, grf_mot_file)
% fileout = ik_grf(mot_file, grf_mot_file)
% 
% This function will combine the GRF data from a _GRF.MOT file with that
% from a motion file (MOT file e.g. that produced by IK analysis) so that the ground
% reaction vector can be visualised in the Opensim animation
%
% Input - mot_file - the motion file (*.MOT) containing the generalised coordinates
%         grf_mot_file - the motion file (*_.MOT) containing the ground reaction
%         force vector information (generally made using c3d2trc.m)
%         
%         output - the filename of the output MOT file
%
% Authors: Rod Barrett and Glen Lichtwark (Griffith University)
%

if nargin == 1
    error('Need to input 2 files (MOT and _GRF.MOT) or none so that you can choose the files');
end
if nargin < 1
[fname1, pname1] = uigetfile('*ik.mot', 'Select motion file (typically IK) -');
[fname2, pname2] = uigetfile('*grf.mot', 'Select motion file containing GRF data-');
else [pname1, name1, ext1] = fileparts(ik_mot_file);
    [pname2, name2, ext2] = fileparts(grf_mot_file);
    fname1 = [name1 ext1];
    fname2 = [name2 ext2];
    pname1 = [pname1 '\'];
    pname2 = [pname2 '\'];
end

filename1 = [pname1 fname1];
filename2 = [pname2 fname2];

% load files
[header1 data1] = hdrload(filename1);
[header2 data2] = hdrload(filename2);

if ~isempty(data1<2)
   [data, result]= readtext(filename1, '\t', '', '', 'empty2NaN') ;
   a = find(sum(result.numberMask,2)<size(result.emptyMask,2));
   b = find(sum(result.numberMask,2)==size(result.emptyMask,2));
   data1 = reshape([data{b,:}],length(b),size(result.emptyMask,2));
   header1 = [];
   for i = 1:length(data(a(end),:))
       header1 = [header1 data{a(end),i} ' '];
   end
end

data = [data1 interp1(data2(:,1),data2(:,2:end),data1(:,1),'spline','extrap')];
[nrows, ncols] = size(data);

% write the new header
fileout = [filename1(1:end-4) '_grf.mot'];
[~, name, ~] = fileparts(fileout);

fid = fopen(fileout,'w');

fprintf(fid,'name %s\n',name);
fprintf(fid,'datacolumns %d\n', ncols);  % total number of datacolumns
fprintf(fid,'datarows %d\n',nrows); % number of datarows
fprintf(fid,'range %f %f\n',data(1,1),data(end,1)); % range of time data
fprintf(fid,'endheader\n');

% find the spaces in the last row of both headers, which represents the
% columnn headings and place tabs in between rather than spaces
a = find(isspace(header1(end,:))==1);
b = find(isspace(header2(end,:))==1);

% determine if there are number in front of the GRF column headers, if so
% we will take these out so this can be used in the visualisation
if strcmp(header2(end,b(1)+1),'1')
    st = 3;
else st = 1;
end

column_headings = header1(end,1:a(1)-1);
for i = 1:length(a)-1
    column_headings = [column_headings '\t' header1(end,a(i)+1:a(i+1)-1)];
end
column_headings = [column_headings '\t' header1(end,a(end)+1:end) '\t' header2(end,b(1)+st:b(2)-1)];
for i = 2:length(b)-1
    column_headings = [column_headings '\t' header2(end,b(i)+st:b(i+1)-1)];
end
column_headings = [column_headings '\t' header2(end,b(end)+1:end) '\n'];
% write the column headings to file and the data
fprintf(fid,column_headings);

data_format = [];
for i = 1:ncols
    data_format = [data_format '%20.6f\t'];
end
data_format = [data_format '\n'];

fprintf(fid,data_format,data');

fclose(fid);
