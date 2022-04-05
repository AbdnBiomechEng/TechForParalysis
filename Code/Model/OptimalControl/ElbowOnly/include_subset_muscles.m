function include_subset_muscles(full_model_file,new_model_file)
% A function that changes which muscles are included in the model, by
% creating a modified model structure to include a subset of the muscles
%
% The script get_muscle_names.m could be useful here, to find out the
% indices of the muscles in model_struct.mat

% Load the original model structure, that includes all the muscles
full_model = load(full_model_file);
model = full_model.model;

% Get muscle names
musnames = cell(model.nMus,1);
elemnames = cell(model.nMus,1);
for imus=1:model.nMus
    elemnames{imus} = model.muscles{imus}.osim_name;
    musnames{imus} = elemnames{imus}(1:end-2);
end

unique_musnames = unique(musnames,'stable');

% Select the muscles we want to include
[muscle_indices,~] = listdlg('PromptString','Select muscles to include:','ListString',unique_musnames);

selected_mus = unique_musnames(muscle_indices);
elem_indices = find(contains(elemnames,selected_mus));

model.nMus = length(elem_indices);

model.muscles = model.muscles(elem_indices);

% Get names of muscle subset
selected_elemnames = elemnames(elem_indices);

% Select the paralysed muscles
[par_muscle_indices,~] = listdlg('PromptString','Select paralysed muscles:','ListString',selected_elemnames);
par_muscles = model.muscles(par_muscle_indices);
par_muscles = cellfun(@(x) setfield(x,'maxact',0), par_muscles,'UniformOutput',false);
model.muscles(par_muscle_indices) = par_muscles;

% Select the FES muscles
[FES_muscle_indices,~] = listdlg('PromptString','Select muscles for FES:','ListString',selected_elemnames);
FES_muscles = model.muscles(FES_muscle_indices);
FES_muscles = cellfun(@(x) setfield(x,'maxact',0.5), FES_muscles,'UniformOutput',false);
model.muscles(FES_muscle_indices) = FES_muscles;

% Save the modified model structure
save(new_model_file,'model'); 




