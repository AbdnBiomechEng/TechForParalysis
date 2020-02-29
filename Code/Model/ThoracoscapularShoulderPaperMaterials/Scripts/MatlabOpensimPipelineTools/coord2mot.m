function fileout = coord2mot(D,fileout)
% function fileout = coord2mot(D)
% Function to write MOT file (fileout) given a set of 'coordinates' in structure D.
% The format for D is as follows -
% D.indegress = 'yes' (In degrees specifies if coordinate measures are in
% degrees (defaults to 'yes' if this is not specified)
% D.time = time points (n x 1 array)
% D.coord1 = coordinate 1 (must be the same length as time array)
% D.coord2 = coordinate 2 
% etc
% The structure field names with array data must be the same as the 
% coordinates in the model. 

% first checkt to see if degrees are specified or not - if not set to true
% ('yes')
if isfield(D,'indegrees')
   indegrees = D.indegrees;
   D = rmfield(D,'indegrees'); % remove this field from structure
else indegrees = 'yes';
end

% now check to make sure there is a field called time with time data in it
if isfield(D,'time')
    [m,n] = size(D.time);
    if n~= 1
        error('Time field must be an n x 1 array')
    end
    data_matrix = D.time;
    data_format = '%15.6f\t';
    header_columns = 'time\t';
    D = rmfield(D,'time'); % remove from structure
end

% find remaining fields and add to matrix to be written in file
fnames = fieldnames(D);

for i = 1:length(fnames)
   data_matrix = [data_matrix D.(fnames{i})]; 
   data_format = [data_format '%15.6f\t'];
   % catch activation or fibre length values in the case of states files
   if ~isempty(strfind(fnames{i},'activation'))
       newname = fnames{i};
       newname(strfind(fnames{i},'activation')-1) = '.';
       fnames{i} = newname;
   end
   if ~isempty(strfind(fnames{i},'fiber_length'))
       newname = fnames{i};
       newname(strfind(fnames{i},'fiber_length')-1) = '.';
       fnames{i} = newname;
   end
   if strcmp(fnames{i}(1),'N') && ~isempty(str2num(fnames{i}(2)))
       fnames{i}(1:end-1) = fnames{i}(2:end);
       fnames{i}(end) = [];
   end
   header_columns = [header_columns fnames{i} '\t']; 
end

fid = fopen(fileout,'w');

% write the header information
fprintf(fid,'Coordinates\nversion=1\n');
fprintf(fid,'nRows=%d\n', m);  % total # of rows
fprintf(fid,'nColumns=%d\n',length(fnames)+1); % number of columns
fprintf(fid,'indegrees=%s\n',indegrees); % range of time data
fprintf(fid,'\nUnits are S.I. units (second, meters, Newtons, ...)\nAngles are in degrees.\n\n');
fprintf(fid,'endheader\n');
fprintf(fid,[header_columns(1:end-1) 'n']);

% write the data
fprintf(fid,[data_format(1:end-1) 'n'],data_matrix');

% close the file
fclose(fid);