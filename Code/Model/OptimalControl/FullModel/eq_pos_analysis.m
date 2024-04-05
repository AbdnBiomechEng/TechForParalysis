% Script to plot states and muscle forces of the initial equilibrium
% position for the simplified model

%res_eq = load('C:\Users\s04db9\Code\TechForParalysis\Code\Model\OptimalControl\FullModel\small_range');
%load(res_eq.Result.model);
%model = res_eq.Result.model;

% res_eq = load('dsem_out_new_forward.mat');
% x = res_eq.xout(end,:)';
% u = res_eq.u;

res_eq = load('equilibrium');
x = res_eq.x;
u = zeros(137,1);

modelstr = load('extended_workspace_model.mat');
model = modelstr.model;

% for imus=[37:40 48:51 95:102]  % posterior deltoid, anterior deltoid, pec major
%     model.muscles{imus}.PEEslack = 1.3;
% end

das3('Initialize',model);

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

lce_eq = x(2*ndof+1:2*ndof+nmus);
act_eq = x(2*ndof+nmus+1:end);
fmus_eq = das3('Muscleforces', x);
mus_lengths = das3('Musclelengths', x);
angles_eq = x(1:ndof)*180/pi;

SEE_elong = (mus_lengths - x(2*ndof+(1:nmus)).*LCEopt - SEEslack)./SEEslack;

moments = das3('Jointmoments', x);
momentarms = das3('Momentarms', x);

figure; 
subplot(4,1,1); plot(act_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle activations');
subplot(4,1,2); plot(lce_eq); hold on; plot(PEEslack,'o');
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Normalised fibre lengths');
subplot(4,1,3); plot(fmus_eq); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle forces (N)');
subplot(4,1,4); plot(angles_eq,'o'); 
xticks(1:ndof); xticklabels(dofnames); xtickangle(45); title('Angles (degrees)');

figure; % plot moments for clavicular and scapular dofs
for idof=4:9
    subplot(6,1,idof-3);
    momarms_dof = full(momentarms(:,idof));
    moments_dof = momarms_dof.*fmus_eq;
    bar(moments_dof);
    xticks(1:nmus); xticklabels(musclenames); xtickangle(45);
    title(dofnames{idof});
    ylabel('Moment (Nm)');
end

figure; % plot moments for humeral and forearm dofs
for idof=10:14
    subplot(5,1,idof-9);
    momarms_dof = full(momentarms(:,idof));
    moments_dof = momarms_dof.*fmus_eq;
    bar(moments_dof);
    xticks(1:nmus); xticklabels(musclenames); xtickangle(45);
    title(dofnames{idof});
    ylabel('Moment (Nm)');
end

fprintf('\n\nDOF               angle(deg)    limits (deg)            ang.vel(deg/s)   moment(Nm)  \n');
fprintf('--------------- --------------  ---------------------   -------------- --------------\n');
for i=1:ndof
    fprintf('%-15s %9.3f      %9.3f   %9.3f      %9.3f    %9.3f\n',dofnames{i}, angles_eq(i), 180/pi*range(i,1), 180/pi*range(i,2), 180/pi*x(ndof+i), moments(i));
end


fprintf('\n\nMuscle           Lce/Lceopt   PEEslack    SEE elong     Muscletendon length    Activation    Force(N)  \n');
fprintf('--------------- ------------ ----------- -------------- ---------------------  ----------   ----------\n');
for i=1:nmus
    fprintf('%-15s %9.3f    %9.3f    %9.3f      %9.3f           %9.3f     %9.3f\n',musclenames{i}, x(2*ndof+i), PEEslack(i), SEE_elong(i), mus_lengths(i), x(2*ndof+nmus+i), fmus_eq(i));
end

[f, ~, ~, ~, FGH] = das3('Dynamics', x, 0*x,u);
Fscap = das3('Scapulacontact', x);

aphi=tand(38.55);
atheta=tand(44.37);
Rgt = glenoid_scap;

Fgh0 = Rgt*FGH;  % take glenoid orientation into account
if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
% decompose into polar angles
thetar = asin(-Fgh0(2));
if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
else, phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
end
FGHcontact = (thetar/atheta)^2 + (phir/aphi)^2;