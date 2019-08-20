function make_osimm(filename,dof_names,angles,time)
% Creates motion file for OpenSim
%
% Inputs:
% filename: the name of the motion file, without the extension
% angles: an num_dofs x n or n x num_dofs matrix of angles in degrees, 
% n the number of frames
% time (optional): a 1xn or nx1 vector of time values. If this is not
% provided, the timestep is assumed to be 0.01s.
%
% Dimitra Blana, February 2012

[nrows,ncolumns]=size(angles);
if ncolumns~=length(dof_names)
    angles=angles';
    nrows = ncolumns;
end

if nargin <4
    time = 0.01:0.01:0.01*nrows;
end

if size(time,2)~=1, time = time'; end
if size(time,1)~=nrows
    errordlg('The time vector does not have the same length as the angle data.','Dimension error');
    return;
end
data = [time angles*180/pi];  

% create motion file
% the header of the motion file is:
%
% <motion name>
% nRows=x
% nColumns=y
% endheader
% time dofnames
%

dofnames = 'time';
dofstr = '%f';
for idof = 1:length(dof_names)
    dofnames= [dofnames '  ' dof_names{idof}];
    dofstr =  [dofstr '  %f'];
end
dofstr = [dofstr '\n'];

fid = fopen([filename '.mot'],'wt');
fprintf(fid,'%s\n',filename);
fprintf(fid,'%s%i\n','nRows=',nrows);
fprintf(fid,'%s\n',['nColumns=' num2str(length(dof_names)+1)]);
fprintf(fid,'%s\n','endheader');
fprintf(fid,'%s\n',dofnames);
fprintf(fid,dofstr,data');
fclose(fid);

