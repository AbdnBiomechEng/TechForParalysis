function XML = grf2xml(data,varargin)
% FUNCTION XML = GRF2XML(data, fileout)
%
% Function to generate an appopriate XML file for inverse dynamics etc in
% the Opensim software. This file is generated from the information
% calculated in C3DTRC and relates to the MOT file generated. 
%
% INPUT -   REQUIRED -
%           ExternalLoadNames - cell array containing names of the external
%               loads to be used
%           AppliedToBodies - cell array containing names of bodies that
%               the forces are to be applied to
%           GRFFile - the name of the GRF file (a MOT or STO file) which
%               contains the force data to be applied to the model
%
%           OPTIONAL - 
%           MOTFile - motion file (.mot) or storage file (.sto) containing the
%               model kinematics used to transform a point expressed in ground to the
%               body of force application.If the point is not expressed in ground, the
%               point is not transformed
%           LowPassFilterForKinematics - cutoff frequency for low pass
%               filter frequency (default -1, no filter)
%           ForceExpressedinBody - Cell array containing the names of the 
%               body the force is expressed in (default is {'ground'} for 
%               each external load) - needs to have same number of cells as
%               ExternalLoadNames
%           PointExpressedinBody - Cell array containing the names of the 
%               body the force point is expressed in (default is {'ground'} for 
%               each external load) - needs to have same number of cells as
%               ExternalLoadNames
%           ForceIdentifier - Identifier (string) to locate the force to be 
%               applied in the data source (default - n_ground_force_v
%               where n is the order number of the external load) If default
%               is not used, then the same number of cells as
%               ExternalLoadNames needs to be provided with each cell
%               containing the string of the column header.
%               
%           PointIdentifier - Identifier (string) to locate the force to be 
%               applied in the data source (default - n_ground_force_p
%               where n is the order number of the external load) If default
%               is not used, then the same number of cells as
%               ExternalLoadNames needs to be provided with each cell
%               containing the string of the column header.
%           TorqueIdentifier - Identifier (string) to locate the force to be 
%               applied in the data source (default - n_ground_torque
%               where n is the order number of the external load) If default
%               is not used, then the same number of cells as
%               ExternalLoadNames needs to be provided with each cell
%               containing the string of the column header.
%            DataSourceName - Name of the data source (Storage) that will 
%               supply the force data (defaults to 'Unassigned')
%           
%           OutputFile - the name of the XML file that is output (if empty
%               this defaults to the same name as the MOT file.
%
% OUTPUT -  XML - (optional) structure containing the relevant structure
%           for generating the XML file
%
%   e.g. grf2xml('ExternalLoadNames',{'ExternalForce_1','ExternalForce_2'},...
%            'AppliedtoBodies',{'calcn_r','calcn_l'},'GRFFile','Trial01_grf.mot',...
%            'MOTFile',Trial01_ik.mot);
%
%           Creates the a _grf.MOT file for OpenSim
%
% Written by Glen Lichtwark (University of Queensland)

ExternalLoadNames = [];
AppliedToBodies = [];
GRFFile = [];

MOTFile = [];
LowPassFilterForKinematics = -1;
ForceExpressedinBody = [];
PointExpressedinBody = [];
ForceIdentifier = [];
PointIdentifier = [];
TorqueIdentifier = [];
DataSourceName = [];

OutputFile = [];

if ~isempty(varargin)
    if rem(length(varargin),2)
        error('Incorrect input arguments - must specify property and input')
    end
    for i = 1:2:length(varargin)
       n = varargin{i+1};
       eval([varargin{i} '= n;']); 
    end    
end

%setup root of the XML and the ExternalLoads definition
root = 'OpenSimDocument';
        
V.ATTRIBUTE.Version = '20302';
[~, filename, ~] = fileparts(data.TRC_Filename);
V.ExternalLoads.ATTRIBUTE.name = [data.Name filename];

if ~iscell(ExternalLoadNames)
    error('Please use cell array containing the names for each external load to be applied');
end

for i = 1:length(ExternalLoadNames)
    V.ExternalLoads.objects.ExternalForce(i).ATTRIBUTE.name = ExternalLoadNames{i};%['ExternalForce_' num2str(i)];
    V.ExternalLoads.objects.ExternalForce(i).isDisabled = 'false';
    V.ExternalLoads.objects.ExternalForce(i).applied_to_body = AppliedToBodies{i};
    if isempty(ForceExpressedinBody)
        V.ExternalLoads.objects.ExternalForce(i).force_expressed_in_body = 'ground';
    else V.ExternalLoads.objects.ExternalForce(i).force_expressed_in_body = ForceExpressedinBody{i};
    end
    if isempty(PointExpressedinBody)
        V.ExternalLoads.objects.ExternalForce(i).point_expressed_in_body = 'ground';
    else V.ExternalLoads.objects.ExternalForce(i).point_expressed_in_body = ForceExpressedinBody{i};
    end
    if isempty(ForceIdentifier)
        V.ExternalLoads.objects.ExternalForce(i).force_identifier = [num2str(i) '_ground_force_v'];
    else V.ExternalLoads.objects.ExternalForce(i).force_identifier = ForceIdentifier{i};
    end
    if isempty(PointIdentifier)
        V.ExternalLoads.objects.ExternalForce(i).point_identifier = [num2str(i) '_ground_force_p'];
    else V.ExternalLoads.objects.ExternalForce(i).point_identifier = PointIdentifier{i};
    end
    if isempty(TorqueIdentifier)
        V.ExternalLoads.objects.ExternalForce(i).torque_identifier = [num2str(i) '_ground_torque_'];
    else V.ExternalLoads.objects.ExternalForce(i).torque_identifier = TorqueIdentifier{i};
    end
    if isempty(DataSourceName)
        V.ExternalLoads.objects.ExternalForce(i).data_source_name = 'Unassigned';
    else V.ExternalLoads.objects.ExternalForce(i).data_source_name = DataSourceName;
    end
end

if isempty(GRFFile)
    error('Need to assign a GRF motion file');
else V.ExternalLoads.datafile=GRFFile;
end
if isempty(MOTFile)
    V.ExternalLoads.external_loads_model_kinematics_file= ' ';
else V.ExternalLoads.external_loads_model_kinematics_file=MOTFile;
end

V.ExternalLoads.lowpass_cutoff_frequency_for_load_kinematics = LowPassFilterForKinematics;

if isempty(OutputFile)
    OutputFile = [data.C3D_Filename(1:end-4) '_grf.xml'];
end

Pref.StructItem = false;

xml_write(OutputFile, V, root, Pref);

XML = V;

% % we need to rename the ExternalForce so they aren't numbered
% T = readtext(OutputFile);
% fid = fopen(OutputFile,'w');
% for i = 1:length(T)
%     if ~isempty(strfind(T{i},'ExternalForce'))
%         a = strfind(T{i},'ExternalForce');
%         T{i}(a(1)+13) = [];
%     end
%     fprintf(fid,'%s\n',T{i});
% end
% fclose(fid);

end
