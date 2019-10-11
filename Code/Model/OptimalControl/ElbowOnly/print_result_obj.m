load('flexion100-min-effort\output_17.mat')
figure;
subplot(2,1,1); plot(Result.times,Result.x(13,:)'*180/pi); ylim([20 105]); ylabel('elbow angle (degrees)');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,[mean(Result.x(40:43,:)); mean(Result.x(44:50,:))]); 
ylim([0 0.5]); xlabel('time (s)'); legend('triceps','brachialis'); ylabel('activation');
load('flexion-min-effort\output_17.mat')
figure;
subplot(2,1,1); plot(Result.times,Result.x(13,:)'*180/pi); ylim([20 105]); ylabel('elbow angle (degrees)');
title(Result.obj_str) 
subplot(2,1,2); plot(Result.times,[mean(Result.x(40:43,:)); mean(Result.x(44:50,:))]); 
ylim([0 0.5]); xlabel('time (s)'); legend('triceps','brachialis'); ylabel('activation');
load('flexion-no-effort\output_17.mat')
figure;
subplot(2,1,1); plot(Result.times,Result.x(13,:)'*180/pi); ylim([20 105]); ylabel('elbow angle (degrees)');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,[mean(Result.x(40:43,:)); mean(Result.x(44:50,:))]); 
ylim([0 0.5]); xlabel('time (s)'); legend('triceps','brachialis'); ylabel('activation');
load('no-flexion-min-effort\output_17.mat')
figure;
subplot(2,1,1); plot(Result.times,Result.x(13,:)'*180/pi); ylim([20 105]); ylabel('elbow angle (degrees)');
title(Result.obj_str)
subplot(2,1,2); plot(Result.times,[mean(Result.x(40:43,:)); mean(Result.x(44:50,:))]); 
ylim([0 0.5]); xlabel('time (s)'); legend('triceps','brachialis'); ylabel('activation');
 