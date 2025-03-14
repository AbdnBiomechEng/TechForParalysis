% script to plot the optimisation results
% remember to update the muscle elements selected for plotting, if they
% were changed in include_subset_muscles.m

%load('reach\A_40.mat');
%load('reach\A_paralysed_40.mat');
load('reach\A_FES_40.mat');

% figure;
% subplot(2,1,1); plot(Result.times,Result.x(13:14,:)'*180/pi); 
% ylabel('angle (degrees)');
% legend('flexion/extension','pronation/supination');
% title(Result.obj_str)
% subplot(2,1,2); plot(Result.times,...
%    [mean(Result.x(2*Result.ndof+Result.nmus+(1:3),:)); ... % biceps
%     mean(Result.x(2*Result.ndof+Result.nmus+(4:7),:)); ... % triceps
%     mean(Result.x(2*Result.ndof+Result.nmus+(8:14),:)); ... % brachialis
%     mean(Result.x(2*Result.ndof+Result.nmus+(15:17),:)); ... % brachioradialis
%     mean(Result.x(2*Result.ndof+Result.nmus+(18:19),:)); ... % pronator teres
%     ]); 
% xlabel('time (s)');
% ylabel('activation');
% legend('biceps','triceps','brachialis','brachioradialis','pronator teres');

figure;
subplot(2,1,1); plot(Result.times,Result.x(13:14,:)'*180/pi); 
ylabel('angle (degrees)');
legend('flexion/extension','pronation/supination');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,...
   [mean(Result.x(2*Result.ndof+Result.nmus+(1:3),:)); ... % biceps
    mean(Result.x(2*Result.ndof+Result.nmus+([4:12, 33:37]),:)); ... % triceps
    mean(Result.x(2*Result.ndof+Result.nmus+(13:19),:)); ... % brachialis
    mean(Result.x(2*Result.ndof+Result.nmus+(20:22),:)); ... % brachioradialis
    mean(Result.x(2*Result.ndof+Result.nmus+([23:24, 30:32]),:)); ... % pronators
    mean(Result.x(2*Result.ndof+Result.nmus+(25:29),:)); ... % supinator
    ]);
xlabel('time (s)');
ylabel('activation');
legend('biceps','triceps','brachialis','brachioradialis','pronators','supinator');
