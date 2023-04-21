% A script that reads the muscle names from the model structure model.struct.mat
% (Can be useful if we need to know the muscle order in the structure.)

load('model_struct.mat')
musnames = cell(model.nMus,1);
for imus=1:model.nMus
    musnames{imus} = model.muscles{imus}.osim_name;
end