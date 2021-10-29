% A script that changes which muscles are included in the model, by
% creating a modified model structure to include a subset of the muscles
%
% The script get_muscle_names.m could be useful here, to find out the
% indices of the muscles in model_struct.mat

% Load the original model structure, that includes all the muscles
load('model_struct.mat');

% Print out all muscle element names
das3('Initialize',model);
for i_mus = 1:model.nMus
    name = das3('Musclename', i_mus);
    fprintf('Element %d: %s\n', i_mus, name);
end

% Set the indices of muscles we want to include
muscle_indices = [83:85, 86:89, 109:115]; % biceps, triceps, brachialis
model.nMus = length(muscle_indices);

all_muscles = model.muscles;
model.muscles = all_muscles(muscle_indices);

% Save the modified model structure
save('model_struct_elbow.mat','model'); 




