function write_sto_file(data, filename)
% function to write a sto file based on an struture array where each field
% is an array to be added --> note each array must be the same length

if nargin < 2
    filename = 'output.sto';
end

fname = fieldnames(data);
[~,filestring,ext] = fileparts(filename);

%open the file
fid_1 = fopen(filename,'w');

D = data.time;
hd = ['time' '\t'];
fm = '%6.6f\t';

b = find(~strcmp('time',fname));

for i = 1:length(fname)-1
    if isstruct(data.(fname{b(i)}))
        f2names = fieldnames(data.(fname{b(i)}));
        for j = 1:length(f2names)
            D = [D data.(fname{b(i)}).(f2names{j})];
            hd = [hd fname{b(i)} '.' f2names{j} '\t'];
        end
    else D = [D data.(fname{b(i)})];
        hd = [hd fname{b(i)} '\t'];
    end
    fm = [fm '%6.6f\t'];
end

hd = [hd(1:end-1) 'n'];
fm = [fm(1:end-1) 'n'];

% first write the header data
fprintf(fid_1,'%s\n',filestring);
fprintf(fid_1,'version=1\n');
fprintf(fid_1,'nRows=%d\n', size(D,1)); 
fprintf(fid_1,'nColumns=%d\n', size(D,2)); 
fprintf(fid_1,'inDegrees=no\n'); 
fprintf(fid_1,'endheader\n'); 
fprintf(fid_1, hd);
% then write the output marker data
fprintf(fid_1,fm,D');

%close the file
fclose(fid_1);