function print_maxforces(model_struct)
% Print out maximum forces in plaintext file from model structure .mat file

model = load(model_struct);

% Get muscle names and fmax
elemnames = cell(model.model.nMus,1);
maxforces = zeros(model.model.nMus,1);
for imus=1:model.model.nMus
    elemnames{imus} = model.model.muscles{imus}.osim_name;
    maxforces(imus) = model.model.muscles{imus}.fmax;
end

% Create a table of element names and max forces
mus_forces = table(elemnames,maxforces);

writetable(mus_forces, [model_struct '_MaxForces.csv']);