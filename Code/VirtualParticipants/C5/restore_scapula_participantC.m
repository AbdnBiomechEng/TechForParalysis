load('simplified_model_struct_partC.mat')

for imus=1:9
    model.muscles{imus}.fmax = 5000;
end

save simplified_model_struct_partCscap model