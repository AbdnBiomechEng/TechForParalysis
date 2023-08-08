full_model = load('model_struct.mat');
simplified_model = load('simplified_arm_model_struct.mat');

% Copy muscle mass for first 36 elements
for imus=1:36
    simplified_model.model.muscles{imus}.mass = full_model.model.muscles{imus}.mass;
end

model = simplified_model.model;
save simplified_arm_model_struct model