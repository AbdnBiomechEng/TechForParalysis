% Generate scaled models for each participant

osimfile = "1.model_scaled_bones.osim";

% First, get the approximate muscle lenghts and GH force vectors using
% das3_polynomials.m. 

mydir = "tempA";
musclepolyfile = "musclepolyA";
GHpolyfile = "GHpolyA";

% Run this in the DAS3 repository:
%das3_polynomials(osimfile,mydir,musclepolyfile,GHpolyfile);
% Then move musclepolyfile and GHpolyfile into them participant folder.

% Create matlab structure using the musclepoly and GH poly files by calling
% das3_readosim.m

model = das3_readosim(osimfile,musclepolyfile,GHpolyfile);
save("model_struct_A", "model");

% If not all muscles are needed, keep a subset of muscles:
new_model_file = "model_struct_A";
include_subset_muscles("model_struct_A",new_model_file);

% This model is a scaled model with no change in maximum muscle force and
% no muscle weakness. It is, in a sence, the "healthy" version of this
% participant.

% To generate an "injured" version, we will need to also scale
% the maximum forces:

osimfile_injured = "2.model_scaled_bonesAndMuscles.osim";
model = das3_readosim(osimfile_injured,musclepolyfile,GHpolyfile);

% and set the maximum voluntary excitations based on testing:
max_controls_t = readtable("2.max_control_list.csv");
max_controls = array2table(table2array(max_controls_t(:,2))');
max_controls.Properties.VariableNames = table2array(max_controls_t(:,1))';

for imus=1:model.nMus
    musname = model.muscles{imus}.osim_name;
    model.muscles{imus}.maxact = max_controls{1,musname};
end

save("model_struct_A_injured", "model");

% Now we can decide which muscles to include, which are paralysed, and which to stimulate via
% FES. 
new_model_file = "model_struct_A_injured";
paralysed_model_file = "model_struct_A_paralysed";
FES_model_file = "model_struct_A_FES";
include_subset_muscles("model_struct_A_injured",new_model_file,paralysed_model_file,FES_model_file);



