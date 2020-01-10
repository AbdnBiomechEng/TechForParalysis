function out = load_sto_file(filename)

% function data = load_sto_file(filename,delimiters)
%
% This function loads either STO or MOT files and stores each column of
% data in structure array with the field name which is the columning
% heading in the file. It discards any other information from the header of 
% the file.
%
% Input: filename - the STO or MOT filename
%
% Author: Glen Lichtwark 
% Last Modified: 17/11/2008

if nargin < 1
    [fname, pname] = uigetfile('*.*', 'File to load - ');
    file = [pname fname];
else file = filename;    
end

[file_data,s_data]= readtext(file, '\t', '', '', 'empty2NaN');

% search the numerical data (in s_data.numberMask) to find when the block 
% of data starts

a = find(abs(diff(sum(s_data.numberMask,2)))>0);
[m,n] = size(file_data);

% create an array with all of the data
num_dat = [file_data{a(end)+1:end,1:sum(s_data.numberMask(a(end)+1,:),2)}];

% reshape to put back into columns and rows
data = reshape(num_dat,m-a(end),sum(s_data.numberMask(a(end)+1,:),2));

% now find the column headings if there are any
if sum(s_data.stringMask(a(end),:)) == sum(s_data.numberMask(a(end)+1,:))
    data_label = file_data(a(end),:);
    b = a(end)-1;
else b = a(end);
end

% go through the data labels and find any that are duplicates (this occurs
% in the ground reaction force data where each forceplate has the same
% column headings) and add a number to distinguish the duplicates.

for i = 1:length(data_label)
    tf = strcmp(data_label(i),data_label);
    c = find(tf>0);
    if length(c) > 1
        for j = 1:length(c)
            data_label(c(j)) = cellstr([data_label{c(j)} num2str(j)]);
        end
    end
end

% now create the output structure with the field names from the data labels
% and the corresponding data from the columns of the data array
for i = 1:length(data_label)
    f_name = data_label{i};
    % find any spaces and replace with underscore
    e = findstr(f_name, ' ');
    if ~isempty(e)
        f_name(e) = '_';
    end
    e = findstr(f_name, '.');
    if ~isempty(e)
        f_name(e) = '_';
    end
    if ~isempty(str2num(f_name(1)))
        f_name = ['N' f_name];
    end
    out.(f_name) = data(:,i);
end

