% script to plot the optimisation results
% remember to update the muscle elements selected for plotting, if they
% were changed in include_subset_muscles.m

filename = 'flexion90/output_5.mat';
load(filename);

figure;
subplot(2,1,1); plot(Result.times,Result.x(13,:)'*180/pi); ylabel('elbow angle (degrees)');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,...
    [mean(Result.x(28+size(Result.u,1)+(4:7),:)); ...
    mean(Result.x(28+size(Result.u,1)+(8:14),:)); ...
    mean(Result.x(28+size(Result.u,1)+(1:3),:))]); 
xlabel('time (s)'); legend('triceps','brachialis','biceps'); ylabel('activation');
