% script to plot the optimisation results
% remember to update the muscle elements selected for plotting, if they
% were changed in include_subset_muscles.m

filename = 'flexion5_90/output_IPOPT_12.mat';
load(filename);

figure;
subplot(2,1,1); plot(Result.times,Result.x(13:14,:)'*180/pi); 
ylabel('angle (degrees)');
legend('flexion/extension','pronation/supination');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,...
   [mean(Result.x(2*Result.ndof+Result.nmus+(1:3),:)); ... % biceps
    mean(Result.x(2*Result.ndof+Result.nmus+(4:7),:)); ... % triceps
    mean(Result.x(2*Result.ndof+Result.nmus+(8:14),:)); ... % brachialis
    mean(Result.x(2*Result.ndof+Result.nmus+(15:17),:)); ... % brachioradialis
    mean(Result.x(2*Result.ndof+Result.nmus+(18:19),:)); ... % pronator teres
    ]); 
xlabel('time (s)');
ylabel('activation');
legend('biceps','triceps','brachialis','brachioradialis','pronator teres');

