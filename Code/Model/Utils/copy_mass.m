function copy_mass(model_struct)

% full_model = load('model_struct.mat');
% simplified_model = load('simplified_arm_model_struct.mat');
% 
% % Copy muscle mass for first 36 elements
% for imus=1:36
%     simplified_model.model.muscles{imus}.mass = full_model.model.muscles{imus}.mass;
% end
% 
% model = simplified_model.model;
% save simplified_arm_model_struct model

full_model = load('model_struct.mat');
simplified_model = load(model_struct);

% Copy muscle mass for first 83 elements
for imus=1:83
    simplified_model.model.muscles{imus}.mass = full_model.model.muscles{imus}.mass;
end

% For short head of biceps, sum the mass of the two short heads of biceps
simplified_model.model.muscles{84}.mass = full_model.model.muscles{84}.mass + full_model.model.muscles{85}.mass;

% Now the rest
for imus=86:138
    simplified_model.model.muscles{imus-1}.mass = full_model.model.muscles{imus}.mass;
end

model = simplified_model.model;
save(model_struct,'model');