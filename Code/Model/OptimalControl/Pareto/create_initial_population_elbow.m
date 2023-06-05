init_pop = [];

for i=1:100
    one_file = load(['feasible_solutions\random_feas_',num2str(i),'.mat']);
    X = one_file.Result.X;
    
    if one_file.Result.status~=0, disp('Solution not found');
    else
        init_pop = [init_pop; X'];
    end
    if size(init_pop,1)==200, break; end
end

save init_pop_elbow_pareto init_pop