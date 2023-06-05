init_pop = [];
kin_cost = zeros(200,1);
act_cost = zeros(200,1);

for i=1:200
    one_file = load(['feasible_solutions\random_feas_',num2str(i),'.mat']);
    kin_cost(i) = one_file.Result.objective_f(1);
    act_cost(i) = one_file.Result.objective_f(2);
    
    X = one_file.Result.X;
    
    if one_file.Result.status~=0, disp('Solution not found');
    else
        init_pop = [init_pop; X'];
    end
end

% Save matrix of feasible solutions
save init_pop_elbow_pareto init_pop

% Plot histograms of kinematic tracking and muscle effort cost functions
figure; subplot(2,1,1); histogram(kin_cost); title('Kinematic tracking cost function');
subplot(2,1,2); histogram(act_cost); title('Muscle effort cost function');

% Scatterplot of kinematic tracking vs muscle effort cost functions
figure; plot(kin_cost,act_cost,'o'); xlabel('Kinematic tracking'); ylabel('Muscle effort');


