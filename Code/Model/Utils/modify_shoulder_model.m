% Script that takes the full DAS3 model (saved in model_struct.mat) and
% simplifies it by keeping a subset of the original muscle elements
%
% For example, the coracobrachialis is represented by 3 elements in the
% full DAS3 model. The simplified model only keeps the 2nd element. This
% element is given a maximum isometric force and mass equal to the sum of the
% forces and masses of all three elements
% 
% The simplified model is saved in simplified_model_struct.mat
%
% The opensim version of the simplified model is das3_simplified.osim
% (created manually, not via this script)

% Load the original model structure, that includes all the muscles
full_model = load('simplified_arm_model_struct.mat');
% If you are updating an existing model, start with that
%simplified_model = load('simplified_model_struct.mat');
% otherwise start with the original model
simplified_model = full_model;

% Get muscle names, fmax and mass
elemnames = cell(full_model.model.nMus,1);
maxforces = zeros(full_model.model.nMus,1);
for imus=1:full_model.model.nMus
    elemnames{imus} = full_model.model.muscles{imus}.osim_name;
    maxforces(imus) = full_model.model.muscles{imus}.fmax;
end

% Create a table of element names and max forces
mus_forces = table(maxforces,'RowNames',elemnames);

% Get muscle names of simplified model
simple_elemnames = cell(simplified_model.model.nMus,1);
for imus=1:simplified_model.model.nMus
    simple_elemnames{imus} = simplified_model.model.muscles{imus}.osim_name;
end

newmodel = simplified_model.model;
newmodel.muscles = cell(1,1);
imus = 1;

% Now find max forces and mass based on the elements we've decided to keep
% - trap_scap_3: 1->5 
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap_1","trap_scap_2","trap_scap_3","trap_scap_4","trap_scap_5"],1).maxforces);
imus = imus+1;

% - trap_scap7: 6->7  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap_7"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap_6","trap_scap_7"],1).maxforces);
imus = imus+1;

% - trap_scap9: 8->9  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap_9"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap_8","trap_scap_9"],1).maxforces);
imus = imus+1;

% - trap_scap11: 10->11  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap11"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap10","trap_scap11"],1).maxforces);
imus = imus+1;

% - trap_clav_1: 1->2  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_clav_1"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_clav_1","trap_clav_2"],1).maxforces);
imus = imus+1;

% - lev_scap_1: 1->2 
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"lev_scap_1"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["lev_scap_1","lev_scap_2"],1).maxforces);
imus = imus+1;

% - pect_min_2: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pect_min_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pect_min_1","pect_min_2","pect_min_3","pect_min_4"],1).maxforces);
imus = imus+1;

% - rhomboid_3: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"rhomboid_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["rhomboid_1","rhomboid_2","rhomboid_3","rhomboid_4","rhomboid_5"],1).maxforces);
imus = imus+1;

% - serr_ant_2:1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_1","serr_ant_2","serr_ant_3"],1).maxforces);
imus = imus+1;

% - serr_ant_5:4->7  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_5"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_4","serr_ant_5","serr_ant_6","serr_ant_7"],1).maxforces);
imus = imus+1;

% - serr_ant_8: 8->12  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_8"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_8","serr_ant_9","serr_ant10","serr_ant11","serr_ant12"],1).maxforces);
imus = imus+1;

newmodel.muscles(imus:imus+25) = simplified_model.model.muscles(37:62);

newmodel.nMus = imus+25;

model = newmodel;
save('new_simplified_model_struct','model'); 
