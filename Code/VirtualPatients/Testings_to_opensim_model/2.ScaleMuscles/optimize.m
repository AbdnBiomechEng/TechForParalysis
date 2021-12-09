function [mc, J, exitflag, output] = optimize(M_vol, M_stim, r_mus, Fmax, x0, lb, ub)

    % Functions 
    fun = @(mc) (1/prod(mc.^2)); % maximises the values of mc and keeps them similar?
    nonlcon = @(mc) nonlinfun(mc,r_mus,Fmax,M_vol,M_stim); % anonymous function for non-linear constraints with additional arguments

    [mc, J, exitflag, output] = fmincon(fun, ... % 'fun' is the name of the objective function
        x0, ...  % Initial guess, or starting point
        [], [], ... % there are no linear inequality contraints
        [], [], ... % these are A and B in A[x] = B, the linear equality contraint
        lb, ub, ...  % lower and upper bounds on the values
        nonlcon,... % non-linear constraints
        []); % options

    %fprintf('%i\n', mc)

end
