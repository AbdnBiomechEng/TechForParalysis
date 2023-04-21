function find_das3_feas_solutions
    for i=1:1000
        disp(['Iteration ' num2str(i) '...']);
        find_one_solution(['pareto_con/random_feas_',num2str(i)]);
    end
end

function find_one_solution(out_filename)
% This program finds a feasible arm position from a random initial guess, ignoring cost function

modelparams = load('simplemus_model_struct.mat'); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% some other constants
nvar = nstates + nmus;              % number of unknowns 
ncon_eq = 2*ndof + 2*nmus;          % number of constraints due to discretized dynamics
ncon = ncon_eq;

% define DOF names (ensure their order remains consistent with q[] array defined in das3.al)
dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% Initialize the model
das3('Initialize',model);

% lock thorax dofs
lockeddofs = [1,2,3];
lockeddofvalues = [0;0;0];

% these are the range of motion limits 
xlim = das3('Limits')';

% Lower and upper bounds for the optimization variables X
% Bounds are:
%   joint angles 	xlim
%   angular velocities 0
%   CE lengths		0.3 to 1.7
%	active states	0 to 1
%	excitations     0 to 1

L = [xlim(1:ndof,1);            % q
    zeros(ndof,1);              % qdot
    zeros(nmus,1) + 0.3;        % Lce
    zeros(nmus,1);              % active states
    zeros(nmus,1)];             % neural excitations

U = [xlim(1:ndof,2);            % q
    zeros(ndof,1);              % qdot
    zeros(nmus,1) + 1.7;        % Lce
    ones(nmus,1);               % active states
    ones(nmus,1)];              % neural excitations

% Precompute the indices for activations and excitations
iact = 2*ndof+nmus+1:nstates;
iexc = nstates+(1:nmus);

% Further constrain the feasible solution space:
% Low shoulder elevation
L(11) = 0;
U(11) = 30*pi/180;
% Low elbow flexion
L(13) = 0;
U(13) = 50*pi/180;
% Mid pro/supination
L(14) = 80*pi/180;
U(14) = 100*pi/180;
% low activations
U(iact) = 0.2;
U(iexc) = 0.2;

L(lockeddofs) = lockeddofvalues;
U(lockeddofs) = lockeddofvalues;
L(ndof+lockeddofs) = zeros(length(lockeddofs),1);
U(ndof+lockeddofs) = zeros(length(lockeddofs),1);

print_flag = 0;

% Precompute the indices for muscles that cross GH
iGH = [];
for imus=1:nmus
    if(das3('crossGH', imus))
        iGH=[iGH imus];
    end
end

% Set initial guess
X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound

% Parameters required for glenohumeral constraint
aphi=tand(38.55);
atheta=tand(44.37);
Rgt = glenoid_scap;

% Optimisation parameters
MaxIterations = 10000;	% max number of iterations for each optimization
OptimalityTolerance = 1e-3;
FeasibilityTolerance = 1e-3;

% Run optimization
if (MaxIterations > 0)
    print_flag=0;    
    % determine sparsity structure of Jacobian (midpoint discretization) I
    % have verified that we always get same structure by testing with
    % random X
    Jnnz = 1;
    X = L + (U-L).*rand(size(L));		% a random vector of unknowns
    J = conjac(X);
    Jnnz = nnz(J);
    fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n',Jnnz, ncon*nvar, 100*Jnnz/(ncon*nvar));
    Jpattern = double(J~=0);
    fprintf('Nonzero elements in Jpattern: %d\n', nnz(Jpattern));
    %figure; spy(Jpattern);
    
    % Check the derivatives - only set to 1 during code development
    checkderivatives = 0;
    
    if (checkderivatives)
        hh = 1e-7;
        X = L + (U-L).*rand(size(L));
        cf = confun(X);
        cjac = conjac(X);
        cjac_num = zeros(ncon,nvar);
        for i_var=1:nvar
            fprintf('checking derivatives for unknown %4d of %4d\n',i_var,nvar);
            Xisave = X(i_var);
            X(i_var) = X(i_var) + hh;
            cfi = confun(X);
            cjac_num(:,i_var) = (cfi - cf)/hh;
            X(i_var) = Xisave;
        end
        cjac_full = full(cjac);
        % Find the max difference in constraint jacobian
        [maxerr,irow] = max(abs(cjac-cjac_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d, %8.8f percent error\n', maxerr, irow(icol), icol, maxerr/cjac_full(irow(icol),icol));
        keyboard
    end
    
    % Evaluate initial guess
    print_flag=1;
    confun(X0);
    print_flag=0;
    
    % Set up constraints
    cl = zeros(ncon_eq,1);
    cu = zeros(ncon_eq,1);

    funcs.objective = @objfun_ipopt;
    funcs.gradient  = @objgrad_ipopt;
    
    options.cl = cl;
    options.cu = cu;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    
    options.lb = L;
    options.ub = U;

    % IPOPT options
    options.ipopt.print_level = 0;
    options.ipopt.max_iter = MaxIterations;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.tol = OptimalityTolerance;
    options.ipopt.acceptable_constr_viol_tol = FeasibilityTolerance;
    options.ipopt.acceptable_tol = FeasibilityTolerance;
    options.ipopt.constr_viol_tol = FeasibilityTolerance;
    options.ipopt.print_info_string = 'yes';
    options.ipopt.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
    options.ipopt.dual_inf_tol = 1;

    [X, info] = ipopt(X0,funcs,options);
    Result.status = info.status;
    fprintf('IPOPT info status: %d\n', info.status);
    Result.status = info.status;
    Result.info = info;
    Result.X = X;
    Result.X0 = X0;

else		% Skip optimization
    Result.X = X0;
end

obj_str = '';
print_flag=1;
confun(Result.X);
print_flag=0;

% Save this result on a file
Result.ndof = ndof;
Result.nmus = nmus;
Result.states = X(1:nstates);
Result.muscles = model.muscles;
Result.obj_str = obj_str;
Result.MaxIterations = MaxIterations;
Result.OptimalityTolerance = OptimalityTolerance;
Result.FeasibilityTolerance = FeasibilityTolerance;
Result.times = 0.1;
Result.u = X(iexc);
Result.objective_f = objfun(X);
angles = X(1:ndof);

% Calculate muscle lengths and forces, GH and scapula stability
Result.mus_lengths = das3('Musclelengths', Result.states);
Result.mus_forces = das3('Muscleforces', Result.states);
[~, ~, ~, ~, FGH, ~, Result.thor_hum] = das3('Dynamics', Result.states, 0*Result.states,Result.u);
Result.Fscap = das3('Scapulacontact', Result.states)';
Result.FGHcontact = calculate_FGH(FGH);

save([out_filename '.mat'],'Result');
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};
make_osimm(out_filename,dofnames,angles);

% Objective function (not used in optimisation, but calculated afterwards)
    function f = objfun(X)
        
        % First term is mean squared muscle activation
        wf1 = mean(X(iact).^2);

        % Second term is thorax-scapula constraint
        Fscap = das3('Scapulacontact', X(1:nstates));
        wf2 = Fscap(1).^2+Fscap(2).^2;
                
        % Third term is glenohumeral stability constraint
        [~, ~, ~, ~, FGH] = das3('Dynamics', X(1:nstates), zeros(nstates,1),X(iexc));
        FGHcontact = calculate_FGH(FGH);
        wf3 = mean(FGHcontact.^2);
        
        f = [wf1,wf2,wf3];            
    end

% No objective function for IPOPT
    function f = objfun_ipopt(X)
        f = 0;
    end

    function g = objgrad_ipopt(X)
        g = zeros(nvar,1);
    end

% Constraints
    function ceq = confun(X)
        
        % Evaluate dynamics
        f = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iexc));
        f(lockeddofs) = zeros(length(lockeddofs),1);
        f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
        ceq = f';        
                        
        if (print_flag)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq));
        end
    end

    function J = conjacstructure()
        % Returns structure of constraint Jacobian matrix
        J = Jpattern;
    end


% Jacobian of constraints
    function J = conjac(X)
        
        allrows = zeros(Jnnz,1);
        allcols = zeros(Jnnz,1);
        allvals = zeros(Jnnz,1);
        index=1;
        
        % evaluate dynamics         
        [~, dfdx, ~, dfdu] = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iexc));
        dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
        dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
        dfdu(lockeddofs,:) = zeros(length(lockeddofs),nmus);
        dfdu(ndof+lockeddofs,:) = zeros(length(lockeddofs),nmus);
                      
        % which generates two blocks in the Jacobian:
        [r,c,v] = find(dfdx);
        datal = length(v);
        allrows(index:index+datal-1) = r;
        allcols(index:index+datal-1) = c;
        allvals(index:index+datal-1) = v;
        index = index+datal;
        
        [r,c,v] = find(dfdu);
        datal = length(v);
        allrows(index:index+datal-1) = r;
        allcols(index:index+datal-1) = nstates+c;
        allvals(index:index+datal-1) = v;
        index = index+datal;
                    
        if exist('Jpattern','var')
            [r,c,~] = find(Jpattern);
            nz_rows = ismember([allrows(1:index-1),allcols(1:index-1)],[r,c], 'rows');
            J = sparse(allrows(nz_rows),allcols(nz_rows),allvals(nz_rows),ncon, nvar);
        else
            J = sparse(allrows(1:index-1),allcols(1:index-1),allvals(1:index-1),ncon, nvar);
        end
               
    end

% Function to calculate orientation of glenohumeral stability vector
    function FGHcontact = calculate_FGH(FGH)
        Fgh0 = Rgt*FGH;  % take glenoid orientation into account
        if norm(Fgh0), Fgh0 = Fgh0/norm(Fgh0); end
        % decompose into polar angles
        thetar = asin(-Fgh0(2));
        if ~(sqrt(Fgh0(1)^2+Fgh0(3)^2)), phir = 0.0;
        else, phir=asin(Fgh0(3)/sqrt(Fgh0(1)^2+Fgh0(3)^2));
        end
        FGHcontact = (thetar/atheta)^2 + (phir/aphi)^2;
    end


end
