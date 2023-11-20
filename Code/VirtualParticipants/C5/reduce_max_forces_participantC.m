load('simplified_model_struct.mat')
max_force_participantC = readtable('maxforce_participantC.csv');

for imus=1:model.nMus
    musname = model.muscles{imus}.osim_name;
    [~,ii] = ismember(musname,max_force_participantC.Var1);
    model.muscles{imus}.fmax = max_force_participantC.fmax_new(ii);
end

save simplified_model_struct_partC model