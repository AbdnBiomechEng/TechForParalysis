function copy_parameters(model_struct)


full_model = load('full_model.mat');
simplified_model = load(model_struct);

% Copy PEE slack 
for imus=1:137
    simplified_model.model.muscles{imus}.PEEslack = full_model.model.muscles{imus}.PEEslack;
end

model = simplified_model.model;
save(model_struct,'model');