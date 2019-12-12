function traces = readtrc(filename)
%READTRC Reads traces in a .trc file from a Rohde & Schwarz device

fid = fopen(filename);
Rows = textscan(fid, '%s','delimiter', '\n'); % creates a temporary array with the rows
fclose(fid);
    
% Looks for the start of each trace...
TraceStarts=strfind(Rows{1,1},'Scan'); 
TracesIdx = find(~cellfun('isempty', TraceStarts)); %.. and stores the indexes.
    
% Detect numbers in the rows, and discard 'VOID' elements
Columns= cellfun(@(x) textscan(x,'%f','delimiter','\t','CollectOutput',1,'treatAsEmpty',{'VOID'}), Rows{1,1});
Columns= cellfun(@transpose, Columns, 'UniformOutput', 0);
    
% Stores individual traces in a cell
for k=1:size(TracesIdx)-1
    traces{k,1}=cell2mat(Columns(TracesIdx(k)+4:TracesIdx(k+1)-1));
end
traces{k+1,1}=cell2mat(Columns(TracesIdx(size(TracesIdx))+4:size(Rows{1,1})));

% Store the header in the last cell element:
traces{k+2,1}=Rows{1,1}(1:16,1);

end

