N_files = 10;
iact = 2*14+35+1:98;
iLce = 2*14+1:63;
obj_f = zeros(N_files,3);
initial_guess = zeros(N_files,133);
x_out = zeros(N_files,133);
for i=1:N_files
    one_sim = load(['min_act/random_',num2str(i)]);
    if one_sim.Result.status ~=0, continue; end
    obj_f(i,:) = one_sim.Result.objective_f;
    initial_guess(i,:) = one_sim.Result.X0';
    x_out(i,:) = one_sim.Result.X';
end

remove_ind = initial_guess(:,4)==0;
initial_guess(remove_ind,:) = [];
x_out(remove_ind,:) = [];
obj_f(remove_ind,:) = [];

for i=1:35
figure; subplot(2,2,1); histogram(initial_guess(:,iact(i))); ylabel('activation in');
subplot(2,2,3); histogram(x_out(:,iact(i))); ylabel('activation out');
subplot(2,2,2); histogram(initial_guess(:,iLce(i))); ylabel('Lce in');
subplot(2,2,4); histogram(x_out(:,iLce(i))); ylabel('Lce out');
end

for i=1:11
    figure; subplot(2,1,1); histogram(initial_guess(:,i+3)*180/pi); ylabel('angle in');
    subplot(2,1,2); histogram(x_out(:,i+3)*180/pi); ylabel('angle out');
end

% sum_obj = sum(obj_f,2);
% [~,indeces] = sort(sum_obj);
% figure; plot(obj_f(indeces,:));
% 
[~,indeces2] = sort(obj_f(:,1));
figure; plot(obj_f(indeces2,:));
legend('activations','scapula','GH');
