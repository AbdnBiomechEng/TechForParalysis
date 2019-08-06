function Result = simple_opt_ipopt(initialguess,N,mass,stiffness)
% This program optimizes the motion of a mass-spring system to fit motion capture data

% optimization settings
MaxIterations = 2000;
OptimalityTolerance = 1e-12;
FeasibilityTolerance = 1e-6;


% some constants
nstates = 2;                        % mass displacement and linear velocity
nvarpernode = nstates;              % number of unknowns per time node
nvar = nvarpernode * N;				% total number of unknowns
ncon = nstates * (N-1);				% number of constraints due to discretized dynamics
    
% these are the range of motion limits
xlim =  [-100 100];

% Lower and upper bounds for the optimization variables X
% X will contain, for each node, the states and controls.

L = zeros(nvar,1);
U = zeros(nvar,1);
% Bounds are:
%   mass displacement 	-100 to 100 m
%   linear velocity -50 to 50 m/s

for i_node = 0:N-1
    L(i_node*nvarpernode + (1:nvarpernode) ) = [-100;-50];                                 
    
    U(i_node*nvarpernode + (1:nvarpernode) ) = [100; 50];
end

% precompute the indices for the displacement, so we can compute cost function quickly
ix = [];
for i_node=1:N
   ix = [ix nvarpernode*(i_node-1)+1];
end

% example data
t = 0:0.05:0.05*(1000-1);
data = 0.5*sin(t)+(0.2-0.4*rand(size(t)));

Result.mass = mass;
Result.stiffness = stiffness;
Result.input = data;
Result.input_t = t;
tdata = [t; data]';

% make the vector of times for direct collocation
t1 = min(tdata(:,1));
t2 = max(tdata(:,1));
times = t1 + (0:N-1)'*(t2-t1)/(N-1);

% resample the motion capture data onto the direct collocation times
% and put into a long vector to save time when computing cost function
data = interp1(tdata(:,1),tdata(:,2:2),times,'linear','extrap');
datavec = reshape(data', N, 1);

objeval = 0;
coneval = 0;
objlog = [];
conlog = [];
print = 0;

% make an initial guess
if numel(strfind(initialguess, 'mid')) > 0
    X0 = (L + U)/2;						% halfway between upper and lower bound
elseif numel(strfind(initialguess, 'random')) > 0
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
else
    % load a previous solution, initialguess contains file name
    initg = load(initialguess);
    t0 = initg.Result.times;
    x0 = initg.Result.x';
    % interpolate states from initial guess to the current time grid
    if length(t0)>1
        x0 = interp1(t0,x0,times,'linear','extrap');
    else
        x0 = repmat(x0,length(times),1);
    end
    X0 = reshape(x0',nvar,1);
end

% run optimization
if (MaxIterations > 0)
        
    % determine sparsity structure of Jacobian (midpoint discretization)
    % I have verified that we always get same structure by testing with random X
    Jnnz = 1;
    X = L + (U-L).*rand(size(L));		% a random vector of unknowns
    J = conjac(X);
    Jnnz = nnz(J);
    fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n',Jnnz, ncon*nvar, 100*Jnnz/(ncon*nvar));
    Jpattern = double(J~=0);
    fprintf('Nonzero elements in Jpattern: %d\n', nnz(Jpattern));
%     figure; spy(Jpattern);
        
    % check the derivatives
    checkderivatives = 0;

    if (checkderivatives)
        hh = 1e-7;
        X = X0;
        objf = objfun(X);
        grad = objgrad(X);
        cf = confun(X);
        cjac = conjac(X);
        cjac_num = zeros(ncon,nvar);
        grad_num = zeros(nvar,1);
        for ivar=1:nvar
            fprintf('checking derivatives for unknown %4d of %4d\n',ivar,nvar);
            Xisave = X(ivar);
            X(ivar) = X(ivar) + hh;
            cfi = confun(X);
            cjac_num(:,ivar) = (cfi - cf)/hh;
            grad_num(ivar) =   (objfun(X) - objf)/hh;
            X(ivar) = Xisave;
        end

        % find the max difference in constraint jacobian and objective gradient
        [maxerr,irow] = max(abs(cjac-cjac_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d\n', maxerr, irow(icol), icol);
        [maxerr,irow] = max(abs(grad_num-grad));
        fprintf('Max.error in objective gradient:  %8.5f at %d\n', maxerr, irow);
        keyboard
    end
    
    % evaluate initial guess
    print=1;
    objfun(X0);
    print=0;
    fprintf('Initial guess evaluation: f = %f   norm(c) = %f\n',objfun(X0),norm(confun(X0)));
        
    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    options.lb = L;
    options.ub = U;
    % constraints    
    options.cl = zeros(nstates*(N-1),1);
    options.cu = zeros(nstates*(N-1),1);

    % IPOPT options
    options.ipopt.print_level = 5;
    options.ipopt.max_iter = MaxIterations;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.tol = OptimalityTolerance;
	options.ipopt.acceptable_constr_viol_tol = FeasibilityTolerance;
    options.ipopt.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large

    [X, info] = ipopt(X0,funcs,options);
    Result.status = info.status;
    fprintf('IPOPT info status: %d\n', info.status);
    
else		% skip optimization
    X = X0;
end

print=1;
Result.objfun = objfun(X);

% save this result on a file
x = reshape(X,nvarpernode,N);
Result.times = times;
Result.x = x;

save opt_result_ipopt.mat Result

    % objective function
    function f = objfun(X)
                
        % cost function is mean of squared differences to measured data
        f =  mean((X(ix)-datavec).^2);
        
        if (print)
            fprintf('Objfun: %9.5f\n', f);
        end
        
        objeval = objeval + 1;
        objlog = [objlog ; objeval f];

	end


    % gradient of objective function
    function g = objgrad(X)
        g = zeros(nvar,1);                
        g(ix) = 2*(X(ix)-datavec)/numel(ix);
                
    end


    % constraints
    function c = confun(X)
        c = zeros(ncon,1);
        
        % indices of states of node 1
        inode = 1:nstates;
        % indices of constraints at node 1
        ic = 1:(ncon/(N-1));
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            x1 = X(inode);
            x2 = X(inode+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            
            c(ic) = implicit_state_equation((x1+x2)/2,(x2-x1)/h);            
            
            % advance the indices
            inode = inode + nvarpernode;
            ic = ic + (ncon/(N-1));
       end
        
        coneval = coneval + 1;
        conlog = [conlog ; coneval norm(c)];

        if (print)
            fprintf('Norm(c): %9.5f  \n', norm(c));
        end
    end

	
    function J = conjacstructure()
    % returns structure of constraint Jacobian matrix
        J = Jpattern;
    end


    % Jacobian of constraints
    function J = conjac(X)
         
        J = spalloc(ncon, nvar, Jnnz);

        % indices of states of node 1
        ix1 = 1:nstates;
        % indices of constraints at node 1
        ic = 1:nstates;
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            ix2 = ix1 + nvarpernode;
            x1 = X(ix1);
            x2 = X(ix2);
            h = times(i+1) - times(i);		% time interval
            
            % evaluate dynamics violation, and derivatives
            [~, dfdx, dfdxdot] = implicit_state_equation((x1+x2)/2,(x2-x1)/h);
            
            % which generates two blocks in the Jacobian:
             J(ic,ix1) = dfdx/2 - dfdxdot/h;
             J(ic,ix2) = dfdx/2 + dfdxdot/h;
                                   
            % advance the indices
            ix1 = ix1 + nvarpernode;
            ic = ic + nstates;          % there are nstates constraints for each node
        end  
        
    end


    function [f,dfdx, dfdxdot] = implicit_state_equation(s,sdot)        
        % equations describing the movement of a mass-spring system
        f(1) = mass*sdot(2)+stiffness*s(1);
        f(2) = sdot(1) - s(2);
        dfdx = [stiffness 0; 0 -1];
        dfdxdot = [0 mass;1 0];        
    end

end		% end of function simple_opt_ipopt
