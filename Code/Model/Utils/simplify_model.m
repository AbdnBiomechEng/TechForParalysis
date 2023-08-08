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
full_model = load('model_struct.mat');
% If you are updating an existing model, start with that
simplified_model = load('simplified_model_struct.mat');
% otherwise start with the original model
simplified_model = full_model;

% Get muscle names, fmax and mass
elemnames = cell(full_model.model.nMus,1);
maxforces = zeros(full_model.model.nMus,1);
mass = zeros(full_model.model.nMus,1);
for imus=1:full_model.model.nMus
    elemnames{imus} = full_model.model.muscles{imus}.osim_name;
    maxforces(imus) = full_model.model.muscles{imus}.fmax;
    mass(imus) = full_model.model.muscles{imus}.mass;
end

% Create a table of element names and max forces
mus_forces = table(maxforces,'RowNames',elemnames);

% Create a table of element names and mass
mus_mass = table(mass,'RowNames',elemnames);

% Get muscle names of simplified model
simple_elemnames = cell(simplified_model.model.nMus,1);
for imus=1:simplified_model.model.nMus
    simple_elemnames{imus} = simplified_model.model.muscles{imus}.osim_name;
end

newmodel = simplified_model.model;
newmodel.muscles = cell(1,1);
imus = 1;

% Now find max forces and mass based on the elements we've decided to keep
% - trap_scap_4: 1->6 
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap_4"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap_1","trap_scap_2","trap_scap_3","trap_scap_4","trap_scap_5","trap_scap_6"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["trap_scap_1","trap_scap_2","trap_scap_3","trap_scap_4","trap_scap_5","trap_scap_6"],1).mass);
imus = imus+1;

% - trap_scap10: 7->11  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_scap10"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_scap_7","trap_scap_8","trap_scap_9","trap_scap10","trap_scap11"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["trap_scap_7","trap_scap_8","trap_scap_9","trap_scap10","trap_scap11"],1).mass);
imus = imus+1;

% - trap_clav_1: 1->2  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"trap_clav_1"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["trap_clav_1","trap_clav_2"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["trap_clav_1","trap_clav_2"],1).mass);
imus = imus+1;

% - lev_scap_1: 1->2 
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"lev_scap_1"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["lev_scap_1","lev_scap_2"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["lev_scap_1","lev_scap_2"],1).mass);
imus = imus+1;

% - pect_min_2: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pect_min_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pect_min_1","pect_min_2","pect_min_3","pect_min_4"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["pect_min_1","pect_min_2","pect_min_3","pect_min_4"],1).mass);
imus = imus+1;

% - rhomboid_3: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"rhomboid_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["rhomboid_1","rhomboid_2","rhomboid_3","rhomboid_4","rhomboid_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["rhomboid_1","rhomboid_2","rhomboid_3","rhomboid_4","rhomboid_5"],1).mass);
imus = imus+1;

% - serr_ant_2:1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_1","serr_ant_2","serr_ant_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["serr_ant_1","serr_ant_2","serr_ant_3"],1).mass);
imus = imus+1;

% - serr_ant_5:4->7  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_5"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_4","serr_ant_5","serr_ant_6","serr_ant_7"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["serr_ant_4","serr_ant_5","serr_ant_6","serr_ant_7"],1).mass);
imus = imus+1;

% - serr_ant_8: 8->12  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"serr_ant_8"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["serr_ant_8","serr_ant_9","serr_ant10","serr_ant11","serr_ant12"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["serr_ant_8","serr_ant_9","serr_ant10","serr_ant11","serr_ant12"],1).mass);
imus = imus+1;

% - delt_scap_3: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"delt_scap_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["delt_scap_1","delt_scap_2","delt_scap_3","delt_scap_4","delt_scap_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["delt_scap_1","delt_scap_2","delt_scap_3","delt_scap_4","delt_scap_5"],1).mass);
imus = imus+1;

% - delt_scap10: 6->11  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"delt_scap10"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["delt_scap_6","delt_scap_7","delt_scap_8","delt_scap_9","delt_scap10","delt_scap11"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["delt_scap_6","delt_scap_7","delt_scap_8","delt_scap_9","delt_scap10","delt_scap11"],1).mass);
imus = imus+1;

% - delt_clav_2: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"delt_clav_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["delt_clav_1","delt_clav_2","delt_clav_3","delt_clav_4"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["delt_clav_1","delt_clav_2","delt_clav_3","delt_clav_4"],1).mass);
imus = imus+1;

% - coracobr_2: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"coracobr_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["coracobr_1","coracobr_2","coracobr_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["coracobr_1","coracobr_2","coracobr_3"],1).mass);
imus = imus+1;

% - infra_4: 1->6  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"infra_4"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["infra_1","infra_2","infra_3","infra_4","infra_5","infra_6"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["infra_1","infra_2","infra_3","infra_4","infra_5","infra_6"],1).mass);
imus = imus+1;

% - ter_min_2: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"ter_min_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["ter_min_1","ter_min_2","ter_min_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["ter_min_1","ter_min_2","ter_min_3"],1).mass);
imus = imus+1;

% - ter_maj_3: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"ter_maj_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["ter_maj_1","ter_maj_2","ter_maj_3","ter_maj_4"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["ter_maj_1","ter_maj_2","ter_maj_3","ter_maj_4"],1).mass);
imus = imus+1;

% - supra_3: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"supra_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["supra_1","supra_2","supra_3","supra_4"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["supra_1","supra_2","supra_3","supra_4"],1).mass);
imus = imus+1;

% - subscap_6: 1->11  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"subscap_6"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["subscap_1","subscap_2","subscap_3","subscap_4","subscap_5","subscap_6",...
    "subscap_7","subscap_8","subscap_9","subscap10","subscap11"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["subscap_1","subscap_2","subscap_3","subscap_4","subscap_5","subscap_6",...
    "subscap_7","subscap_8","subscap_9","subscap10","subscap11"],1).mass);
imus = imus+1;

% - bic_l: 1->1  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"bic_l"));
imus = imus+1;

% - bic_b_1: 1->2  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"bic_b_1"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["bic_b_1","bic_b_2"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["bic_b_1","bic_b_2"],1).mass);
imus = imus+1;

% - tric_long_2: 1->4  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"tric_long_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["tric_long_1","tric_long_2","tric_long_3","tric_long_4"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["tric_long_1","tric_long_2","tric_long_3","tric_long_4"],1).mass);
imus = imus+1;

% - lat_dorsi_2: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"lat_dorsi_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["lat_dorsi_1","lat_dorsi_2","lat_dorsi_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["lat_dorsi_1","lat_dorsi_2","lat_dorsi_3"],1).mass);
imus = imus+1;

% - lat_dorsi_4: 4->6  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"lat_dorsi_4"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["lat_dorsi_4","lat_dorsi_5","lat_dorsi_6"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["lat_dorsi_4","lat_dorsi_5","lat_dorsi_6"],1).mass);
imus = imus+1;

% - pect_maj_t_3: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pect_maj_t_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pect_maj_t_1","pect_maj_t_2","pect_maj_t_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["pect_maj_t_1","pect_maj_t_2","pect_maj_t_3"],1).mass);
imus = imus+1;

% - pect_maj_t_5: 4->6  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pect_maj_t_5"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pect_maj_t_4","pect_maj_t_5","pect_maj_t_6"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["pect_maj_t_4","pect_maj_t_5","pect_maj_t_6"],1).mass);
imus = imus+1;

% - pect_maj_c_2: 1->2  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pect_maj_c_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pect_maj_c_1","pect_maj_c_2"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["pect_maj_c_1","pect_maj_c_2"],1).mass);
imus = imus+1;

% - tric_med_4: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"tric_med_4"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["tric_med_1","tric_med_2","tric_med_3","tric_med_4","tric_med_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["tric_med_1","tric_med_2","tric_med_3","tric_med_4","tric_med_5"],1).mass);
imus = imus+1;

% - brachialis_4: 1->7  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"brachialis_4"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["brachialis_1","brachialis_2","brachialis_3","brachialis_4","brachialis_5","brachialis_6","brachialis_7"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["brachialis_1","brachialis_2","brachialis_3","brachialis_4","brachialis_5","brachialis_6","brachialis_7"],1).mass);
imus = imus+1;

% - brachiorad_2: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"brachiorad_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["brachiorad_1","brachiorad_2","brachiorad_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["brachiorad_1","brachiorad_2","brachiorad_3"],1).mass);
imus = imus+1;

% - pron_teres_1: 1->1  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pron_teres_1"));
imus = imus+1;

% - pron_teres_2: 1->1  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pron_teres_2"));
imus = imus+1;

% - supinator_3: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"supinator_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["supinator_1","supinator_2","supinator_3","supinator_4","supinator_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["supinator_1","supinator_2","supinator_3","supinator_4","supinator_5"],1).mass);
imus = imus+1;

% - pron_quad_2: 1->3  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"pron_quad_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["pron_quad_1","pron_quad_2","pron_quad_3"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["pron_quad_1","pron_quad_2","pron_quad_3"],1).mass);
imus = imus+1;

% - tric_lat_3: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"tric_lat_3"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["tric_lat_1","tric_lat_2","tric_lat_3","tric_lat_4","tric_lat_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["tric_lat_1","tric_lat_2","tric_lat_3","tric_lat_4","tric_lat_5"],1).mass);
imus = imus+1;

% - anconeus_2: 1->5  
newmodel.muscles(imus) = simplified_model.model.muscles(strcmp(simple_elemnames,"anconeus_2"));
newmodel.muscles{imus}.fmax = sum(mus_forces(["anconeus_1","anconeus_2","anconeus_3","anconeus_4","anconeus_5"],1).maxforces);
newmodel.muscles{imus}.mass = sum(mus_mass(["anconeus_1","anconeus_2","anconeus_3","anconeus_4","anconeus_5"],1).mass);

newmodel.nMus = imus;

model = newmodel;
save('simplified_model_struct','model'); 
