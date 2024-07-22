% Compare previous and new model

% Read in new model structure
modelstr = load('model_struct_clav_scap_orig.mat');
model_new = modelstr.model;

% Read in previous model structure
modelstr = load('model_struct.mat');
model_prev = modelstr.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% get muscle names, SEEslack, PEEslack and lceopt
musnames = cell(nmus,1);
SEEslack_prev = zeros(nmus,1);
PEEslack_prev = zeros(nmus,1);
LCEopt_prev = zeros(nmus,1);
SEEslack_new = zeros(nmus,1);
PEEslack_new = zeros(nmus,1);
LCEopt_new = zeros(nmus,1);
for imus=1:nmus
    musnames{imus} = model_prev.muscles{imus}.osim_name;
    SEEslack_prev(imus) = model_prev.muscles{imus}.lslack;
    PEEslack_prev(imus) = model_prev.muscles{imus}.PEEslack;
    LCEopt_prev(imus) = model_prev.muscles{imus}.lceopt;
    SEEslack_new(imus) = model_new.muscles{imus}.lslack;
    PEEslack_new(imus) = model_new.muscles{imus}.PEEslack;
    LCEopt_new(imus) = model_new.muscles{imus}.lceopt;
end



fprintf('\n\nMuscle           Lceopt_prev   LCeopt_new     PEEslack_prev    PEEslack_new     SEEslack_prev    SEEslack_new  \n');
fprintf('--------------- ------------ ----------- -------------- ---------------------  ----------   ----------\n');
for i=1:nmus
    fprintf('%-15s %9.3f    %9.3f    %9.3f      %9.3f           %9.3f     %9.3f\n',musnames{i}, LCEopt_prev(i), LCEopt_new(i), PEEslack_prev(i), PEEslack_new(i), SEEslack_prev(i), SEEslack_new(i));
end
