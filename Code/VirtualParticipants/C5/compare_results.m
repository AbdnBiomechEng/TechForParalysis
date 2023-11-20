healthy = load('high_reach/healthy_5.mat');
partC = load('high_reach/partC_20.mat');
partC_stim1 = load('high_reach/partC_stim1_20.mat');

nmus = healthy.Result.nmus;
ndof = healthy.Result.ndof;


tric_force_healthy = healthy.Result.mus_forces(21,:);
tric_stim_healthy = healthy.Result.u(21,:);
act_healthy = healthy.Result.x(2*ndof+nmus+1:end,:);


tric_force_partC = partC.Result.mus_forces(21,:);
tric_stim_partC = partC.Result.u(21,:);
act_partC = partC.Result.x(2*ndof+nmus+1:end,:);

tric_force_partC_stim1 = partC_stim1.Result.mus_forces(21,:);
tric_stim_partC_stim1 = partC_stim1.Result.u(21,:);
act_partC_stim1 = partC_stim1.Result.x(2*ndof+nmus+1:end,:);

figure; 
plot(healthy.Result.times,tric_force_healthy); hold on; plot(partC.Result.times, tric_force_partC); plot(partC_stim1.Result.times, tric_force_partC_stim1);
figure;
plot(healthy.Result.times,act_healthy(21,:)); hold on; plot(partC.Result.times, act_partC(21,:)); plot(partC_stim1.Result.times, act_partC_stim1(21,:));