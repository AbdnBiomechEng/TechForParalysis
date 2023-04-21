% res_sh = load('simplemus_dsem\out_kin_halfact_80.mat');
% res_eq = load('simplemus_dsem\out_kin_equalact_5.mat');
% %res_eq = load('simplemus_dsem\out_kin_halfact_40.mat');
% res_nm = load('simplemus_dsem\out_kinact_softscap_softGH_5.mat');
% %res_nm = load('simplemus_dsem\out_kin_halfact_10.mat');
% 
% res_sh = load('simplemus_zero\out_ztwenty_5.mat');
% res_eq = load('simplemus_zero\out_ztwenty_5.mat');
% res_nm = load('simplemus_zero\out_ztwenty_nokin_5.mat');

res_nm = load('hanging.mat');
% res_nm = load('simulation_hanging.mat');
% res_nm.Result.x = res_nm.xout';

modelstruct = load('simplemus_model_struct.mat');

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

% lce_shrug = res_sh.Result.x(2*ndof+1:2*ndof+nmus,:);
% lce_equal = res_eq.Result.x(2*ndof+1:2*ndof+nmus,:);
lce_norm = res_nm.Result.x(2*ndof+1:2*ndof+nmus,end);

% act_shrug = res_sh.Result.x(2*ndof+nmus+1:end,:);
% act_equal = res_eq.Result.x(2*ndof+nmus+1:end,:);
act_norm = res_nm.Result.x(2*ndof+nmus+1:end,end);

figure; subplot(3,1,1); plot(act_norm); %hold on; plot(mean(act_shrug,2)); plot(mean(act_equal,2));
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle activations');
%legend('no kin','kin','kin');

subplot(3,1,2); plot(lce_norm); hold on; plot(PEEslack,'o'); % plot(mean(lce_shrug,2)); plot(mean(lce_equal,2)); 
xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Normalised fibre lengths');

% fmus_shrug = res_sh.Result.mus_forces;
% fmus_equal = res_eq.Result.mus_forces;
% fmus_norm = res_nm.Result.mus_forces;
% 
% subplot(4,1,3); plot(mean(fmus_norm,2)); %hold on; plot(mean(fmus_shrug,2)); plot(mean(fmus_equal,2)); 
% xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle forces (N)');

angles_norm = res_nm.Result.x(1:ndof,:)*180/pi;

res_nm = load('elevation/eq_pos_shoulder6_2.mat');
angles_input = res_nm.Result.input_data{1,:}(2:end)*180/pi;
% angles_shrug = res_sh.Result.x(1:ndof,:)*180/pi;
% angles_equal = res_eq.Result.x(1:ndof,:)*180/pi;

% subplot(4,1,4); plot(mean(angles_norm,2)-angles_input'); %hold on; plot(mean(angles_shrug,2)-angles_input'); plot(mean(angles_equal,2)-angles_input');
% xticks(1:ndof); xticklabels(dofnames); xtickangle(45); title('Output-input angles (degrees)');

subplot(3,1,3); plot(angles_input','*-'); hold on; plot(mean(angles_norm,2),'o-'); %hold on; plot(mean(angles_shrug,2)-angles_input'); plot(mean(angles_equal,2)-angles_input');
xticks(1:ndof); xticklabels(dofnames); xtickangle(45); title('Angles (degrees)'); legend('input','output');

% mus_forces = das3('Muscleforces', xout(end,:)');
% 
% xout = Result.x';
% mus_forces = Result.mus_forces(:,end);
% 
% figure; 
% subplot(3,1,1);
% plot(xout(end,iLce)); hold on; plot(PEEslack,'o'); xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Normalised fibre lengths');
% subplot(3,1,2);
% plot(xout(end,iAct)); xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle activations');
% subplot(3,1,3);
% plot(mus_forces); xticks(1:nmus); xticklabels(musclenames); xtickangle(45); title('Muscle forces (N)');

% figure; subplot(3,1,1); plot(res_nm.Result.x(2*ndof+nmus+1:end,:)'); title('activations');
% subplot(3,1,2); plot(res_nm.Result.u'); title('excitations');
% subplot(3,1,3); plot(res_nm.Result.x(2*ndof+1:2*ndof+nmus,:)','o-'); title('normalised fibre lengths');
% 