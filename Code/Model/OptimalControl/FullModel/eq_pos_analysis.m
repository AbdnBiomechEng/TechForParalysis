% Script to plot states and muscle forces of the initial equilibrium
% position for the simplified model

res_eq = load('simplified_initial_eq.mat');
modelstruct = load('simplified_model_struct.mat');

nmus = modelstruct.model.nMus;
ndof = modelstruct.model.nDofs;

iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

musclenames = cell(nmus,1);
PEEslack = zeros(nmus,1);
for imus=1:nmus
    musclenames{imus} = modelstruct.model.muscles{imus}.name;
    PEEslack(imus) = modelstruct.model.muscles{imus}.PEEslack;
end

dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = modelstruct.model.dofs{idof}.name;
end

lce_eq = res_eq.Result.x(2*ndof+1:2*ndof+nmus);
act_eq = res_eq.Result.x(2*ndof+nmus+1:end);
fmus_eq = res_eq.Result.mus_forces;
angles_eq = res_eq.Result.x(1:ndof)*180/pi;

figure; subplot(4,1,1); plot(act_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle activations');
subplot(4,1,2); plot(lce_eq); hold on; plot(PEEslack,'o');
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Normalised fibre lengths');
subplot(4,1,3); plot(fmus_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle forces (N)');
subplot(4,1,4); plot(angles_eq,'o'); 
xticks(1:ndof); xticklabels(dofnames); xtickangle(45); title('Angles (degrees)');
