% Script to plot states and muscle forces of the initial equilibrium
% position for the simplified model

res_eq = load('equilibrium.mat');
model = res_eq.Result.model;

nmus = model.nMus;
ndof = model.nDofs;

iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

musclenames = cell(nmus,1);
PEEslack = zeros(nmus,1);
SEEslack = zeros(nmus,1);
LCEopt = zeros(nmus,1);
for imus=1:nmus
    musclenames{imus} = model.muscles{imus}.name;
    PEEslack(imus) = model.muscles{imus}.PEEslack;
    SEEslack(imus) = model.muscles{imus}.lslack;
    LCEopt(imus) = model.muscles{imus}.lceopt;
end

dofnames = cell(ndof,1);
range = zeros(ndof,2);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.name;
    range(idof,:) = model.dofs{idof}.range;
end

lce_eq = res_eq.Result.x(2*ndof+1:2*ndof+nmus);
act_eq = res_eq.Result.x(2*ndof+nmus+1:end);
fmus_eq = res_eq.Result.mus_forces;
angles_eq = res_eq.Result.x(1:ndof)*180/pi;

SEE_elong = res_eq.Result.mus_lengths - res_eq.Result.x(2*ndof+(1:nmus)).*LCEopt - SEEslack;

figure; subplot(4,1,1); plot(act_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle activations');
subplot(4,1,2); plot(lce_eq); hold on; plot(PEEslack,'o');
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Normalised fibre lengths');
subplot(4,1,3); plot(fmus_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle forces (N)');
subplot(4,1,4); plot(angles_eq,'o'); 
xticks(1:ndof); xticklabels(dofnames); xtickangle(45); title('Angles (degrees)');


fprintf('\n\nDOF               angle(deg)    limits (deg)            ang.vel(deg/s)\n');
fprintf('--------------- --------------  ---------------------   --------------\n');
for i=1:ndof
    fprintf('%-15s %9.3f      %9.3f   %9.3f      %9.3f\n',dofnames{i}, angles_eq(i), 180/pi*range(i,1), 180/pi*range(i,2), 180/pi*res_eq.Result.x(ndof+i));
end


fprintf('\n\nMuscle           Lce/Lceopt   PEEslack    SEE elong     Muscletendon length    Activation    Force(N)  \n');
fprintf('--------------- ------------ ----------- -------------- ---------------------  ----------   ----------\n');
for i=1:nmus
    fprintf('%-15s %9.3f    %9.3f    %9.3f      %9.3f           %9.3f     %9.3f\n',musclenames{i}, res_eq.Result.x(2*ndof+i), PEEslack(i), SEE_elong(i), res_eq.Result.mus_lengths(i), res_eq.Result.x(2*ndof+nmus+i), res_eq.Result.mus_forces(i));
end

