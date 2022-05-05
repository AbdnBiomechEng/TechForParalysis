% script to make the graph for the BioMedEng22 conference abstract

% We want angles and muscle activations for (a) the paralysed model, and
% (b) the FES model

% The angles are elbow flexion-extension and pronation-supination
% The activations are the mean among the elements of: biceps, triceps, and
% pronator teres

% colours for angles
colours_angles = [         
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    ];

% colours for muscle activations    
colours_act = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    ];


res_FES = load('reach/A_FES_40.mat');
res_par = load('reach/A_paralysed_40.mat');

time = res_FES.Result.times;
angles = [res_FES.Result.x(13:14,:); res_par.Result.x(13:14,:)]'*180/pi;
activations = [
    mean(res_FES.Result.x(2*res_FES.Result.ndof+res_FES.Result.nmus+(1:3),:)); ... % biceps_FES
    mean(res_FES.Result.x(2*res_FES.Result.ndof+res_FES.Result.nmus+([4:12, 33:37]),:)); ... % triceps_FES
    mean(res_FES.Result.x(2*res_FES.Result.ndof+res_FES.Result.nmus+(23:24),:)); ... % pronator teres_FES
    mean(res_par.Result.x(2*res_par.Result.ndof+res_par.Result.nmus+(1:3),:)); ... % biceps_paralysed
    mean(res_par.Result.x(2*res_par.Result.ndof+res_par.Result.nmus+([4:12, 33:37]),:)); ... % triceps_paralysed
    mean(res_par.Result.x(2*res_par.Result.ndof+res_par.Result.nmus+(23:24),:)); ... % pronator teres_paralysed
   ]';

figure;
subplot(2,1,1); 
set(gca, 'ColorOrder',colours_angles, 'LineStyleOrder', {'-','--'});
hold all
for i = 1:4
plot(time,angles(:,i));
end
ylabel('angle (degrees)');
legend('flexion/extension','pronation/supination');
legend boxoff 
title('A')

subplot(2,1,2); 
set(gca, 'ColorOrder',colours_act, 'LineStyleOrder', {'-','--'});
hold all
for i = 1:6
    plot(time,activations(:,i)); 
end
xlabel('Time (s)');
ylabel('Activation');
legend('Biceps','Triceps','Pronator Teres'); 
legend boxoff 
title('B')

