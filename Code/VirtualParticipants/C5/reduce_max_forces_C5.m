load('simplified_model_struct.mat')
max_force_reduction = readtable('maxforce_reduction_C5.csv');

for imus=1:model.nMus
    musname = model.muscles{imus}.osim_name;
    fmax_healthy = model.muscles{imus}.fmax;
    [~,ii] = ismember(musname,max_force_reduction.Var1);
    model.muscles{imus}.fmax = fmax_healthy - max_force_reduction.mc_new(ii)*fmax_healthy;
end

save simplified_model_struct_C5 model