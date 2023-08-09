function find_das3elbow_feas_solutions

if ~exist('feasible_solutions', 'dir')
    mkdir('feasible_solutions');
end

for i=1:200
    disp(['Iteration ' num2str(i) '...']);
    find_one_solution(['feasible_solutions/random_feas_',num2str(i)]);
end
end

function find_one_solution(out_filename)
% This program finds a feasible arm position from a random initial guess, ignoring cost function

% Which model to use
OptSetup.model_file = 'simplified_model_struct.mat';

% Input data
OptSetup.t = [0;0.5];  % the time vector
OptSetup.data_init = zeros(2,8);

% Thoracohumeral angles from 5 to 60 degrees of flexion
OptSetup.data_init(1,4:6) = [90,5,0]*pi/180;
OptSetup.data_init(2,4:6) = [90,60,0]*pi/180;

% Elbow flexion-extension from 60 to 0 degrees
OptSetup.data_init(1,7) = 60*pi/180;
OptSetup.data_init(2,7) = 0*pi/180;

% Pronation-supination from 120 to 90 degrees
OptSetup.data_init(1,8) = 120*pi/180;
OptSetup.data_init(2,8) = 90*pi/180;

% Should angular velocities be zero at the start and end of movement?
% (true/false)
OptSetup.start_at_rest = true;
OptSetup.end_at_rest = true;

% Number of nodes
OptSetup.N = 10;

% optimization settings
Result.OptSetup = OptSetup;

% lock thorax and shoulder dofs
lockeddofs = 1:12;
nlockeddofs = length(lockeddofs);

nmeas_dof = 14;

N = OptSetup.N;

% Load structure with model related variables
modelparams = load(OptSetup.model_file);
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
ncontrols = nmus;
nstates = 2*ndof + 2*nmus;

% some other constants
nvarpernode = nstates + ncontrols;          % number of unknowns per time node
nvar = nvarpernode * N;                     % total number of unknowns
ncon_eq = nstates*(N-1);                    % number of constraints due to discretized dynamics
ncon = ncon_eq;

% define DOF names (ensure their order remains consistent with q[] array defined in das3.al)
dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% Initialize the model
das3('Initialize',model);

% these are the range of motion limits
xlims = das3('Limits')';

% Lower and upper bounds for the optimization variables X
% X will contain, for each node, the states and controls.
L = zeros(nvar,1);
U = zeros(nvar,1);

% Bounds are:
%   joint angles 	xlimdeg
%   angular velocities -40 to 40 degrees/s
%   CE lengths		0 to 2
%	active states	0 to 1
%	neural controls 0 to 1
for i_node = 0:N-1
    L(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,1)-0.1;         % q
        (zeros(ndof,1) - 40);                                           % qdot
        zeros(nmus,1) + 0.3;                                            % Lce
        zeros(nmus,1);                                                  % active states
        zeros(ncontrols,1) ];                                           % neural excitations
    
    U(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,2)+0.1;         % q
        (zeros(ndof,1) + 40);                                           % qdot
        zeros(nmus,1) + 1.7;                                            % Lce
        ones(nmus,1);                                                   % active states
        ones(ncontrols,1)];                                             % neural excitations
end

% Should angular velocities be zero at the start and end of movement?
if OptSetup.start_at_rest
    % First node
    L(ndof+1:2*ndof) = zeros(ndof,1);
    U(ndof+1:2*ndof) = zeros(ndof,1);
end

if OptSetup.end_at_rest
    % Final node
    final_node_start_index = (N-1)*nvarpernode;
    L(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);
    U(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);
end

% Reaching data: 3 thorax angles (zero), 3 thoracohumeral angles, elbow flexion/extension and pronation/supination
t = OptSetup.t;
data_init = OptSetup.data_init;

Result.input = data_init;
Result.input_t = t;

% make the vector of times for direct collocation
t1 = min(t);
t2 = max(t);
times = t1 + (0:N-1)'*(t2-t1)/(N-1);

% resample the movement data onto the direct collocation times
data_N = interp1(t,data_init,times);

% estimate clavicular and scapular angles from thoraco-humeral angles
shoulder_angles = zeros(length(times),9);
for i_t = 1:length(times)
    shoulder_angles(i_t,:) = estimate_shoulder_angles(data_N(i_t,4),data_N(i_t,5),data_N(i_t,6));
end

data = [data_N(:,1:3),shoulder_angles,data_N(:,7:8)];

data_deriv = diff(data)./diff(times);
data_deriv = [data_deriv;data_deriv(end,:)];

% convert into a single column to save time when computing cost function
datavec = reshape(data', [], 1);
data_deriv_vec = reshape(data_deriv', [], 1);

% define the estimated uncertainty in each measured variable (in radians)
datasd = ones(nmeas_dof,1);         % for one node
datasd = repmat(datasd,N,1);		% replicate for all nodes

% precompute the indices for activations and some angles, so we can compute cost function quickly
iact = [];
iLce = [];
iexc = [];

idof = [];
iddof = [];

ilocked_dof = [];
ilocked_ddof = [];
ilocked_meas = [];

for i_node=0:N-1
    iLce = [iLce nvarpernode*i_node+2*ndof+(1:nmus)];
    iact = [iact nvarpernode*i_node+2*ndof+nmus+(1:ncontrols)];
    iexc = [iexc nvarpernode*i_node+2*ndof+2*nmus+(1:ncontrols)];
    
    idof = [idof nvarpernode*i_node+(1:ndof)];
    iddof = [iddof nvarpernode*i_node+ndof+(1:ndof)];
    
    ilocked_dof = [ilocked_dof nvarpernode*i_node+lockeddofs];
    ilocked_ddof = [ilocked_ddof nvarpernode*i_node+ndof+lockeddofs];
    ilocked_meas = [ilocked_meas nmeas_dof*i_node+lockeddofs];
end

% Locked dofs and derivatives of locked dofs stay at measured values
U(ilocked_dof) = datavec(ilocked_meas);
L(ilocked_dof) = datavec(ilocked_meas);
U(ilocked_ddof) = data_deriv_vec(ilocked_meas);
L(ilocked_ddof) = data_deriv_vec(ilocked_meas);

% For the kinematics we are tracking (elbow), we only care about the
% time points where the movement is measured (or set), not all nodes
iElbow = [];
imeas_elbow = [];

for i_node=ceil(linspace(0,N-1,length(t)))
    iElbow = [iElbow nvarpernode*i_node+nlockeddofs+(1:2)];
    imeas_elbow = [imeas_elbow nmeas_dof*i_node+nlockeddofs+(1:2)];
end


% Set initial guess
X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound


% Optimisation parameters
MaxIterations = 10000;	% max number of iterations for each optimization
OptimalityTolerance = 1e-3;
FeasibilityTolerance = 1e-3;

% Run optimization
if (MaxIterations > 0)
    fprintf('Running optimisation...\n');
    print_flag=0;
    % determine sparsity structure of Jacobian (midpoint discretization) I
    % have verified that we always get same structure by testing with
    % random X
    Jnnz = 1;
    X = L + (U-L).*rand(size(L));		% a random vector of unknowns
    J = conjac(X);
    Jnnz = nnz(J);
    %fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n',Jnnz, ncon*nvar, 100*Jnnz/(ncon*nvar));
    Jpattern = double(J~=0);
    %fprintf('Nonzero elements in Jpattern: %d\n', nnz(Jpattern));
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
objfun(Result.X);
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

% Calculate muscle lengths and forces, GH and scapula stability
Result.mus_lengths = das3('Musclelengths', Result.states);
Result.mus_forces = das3('Muscleforces', Result.states);

save([out_filename '.mat'],'Result');

x = reshape(Result.X,[],N);
angles = x(1:ndof,:)';
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};
make_osimm(out_filename,dofnames,angles);

% Objective function (not used in optimisation, but calculated afterwards)
    function f = objfun(X)
        
        % First term of cost function is mean of squared differences to measured data
        f1 = mean(((X(iElbow)-datavec(imeas_elbow))./datasd(imeas_elbow)).^2);
        
        % Second term is mean squared muscle activation
        f2 = mean(X(iact).^2);
        
        f = [f1,f2];
        
        if print_flag
            obj_str = sprintf('Cost functions (not optimised here): %9.5f (fit) , %9.5f (effort)', f1,f2);
            fprintf(obj_str);
            fprintf('\n');
        end
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
        
        % Linear constraints: ceq(x) = 0
        ceq = zeros(ncon_eq,1);
        
        % Calculate constraints due to discretized dynamics
        
        % indices of states and controls of node 1
        ix = 1:nstates;
        ius = nstates+(1:ncontrols);
        % indices of equality constraints at node 1
        iceq = 1:(ncon_eq/(N-1));
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            x1 = X(ix);
            x2 = X(ix+nvarpernode);
            
            u1 = X(ius);
            u2 = X(ius+nvarpernode);
            
            h = times(i+1) - times(i);		% time interval
            
            f = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2);
            
            f(lockeddofs) = zeros(length(lockeddofs),1);
            f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
            ceq(iceq) = f;
            
            % advance the indices
            ix = ix + nvarpernode;
            ius = ius + nvarpernode;
            iceq = iceq + (ncon_eq/(N-1));
        end
        
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
        
        % indices of states and controls of node 1
        ix1 = 1:nstates;
        iu1 = nstates+(1:ncontrols);
        % indices of equality constraints at node 1
        iceq = 1:(ncon_eq/(N-1));
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            ix2 = ix1 + nvarpernode;
            iu2 = iu1 + nvarpernode;
            
            x1 = X(ix1);
            x2 = X(ix2);
            
            u1 = X(iu1);
            u2 = X(iu2);
            
            h = times(i+1) - times(i);		% time interval
            
            [~, dfdx, dfdxdot, dfdu] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2);
            
            dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
            dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
            dfdxdot(lockeddofs,:) = zeros(length(lockeddofs),nstates);
            dfdxdot(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
            dfdu(lockeddofs,:) = zeros(length(lockeddofs),ncontrols);
            dfdu(ndof+lockeddofs,:) = zeros(length(lockeddofs),ncontrols);
            
            
            % which generates four blocks in the Jacobian:
            %J(iceq,ix1) = dfdx/2 - dfdxdot/h;
            [r,c,v] = find(dfdx/2 - dfdxdot/h);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = ix1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            %J(iceq,ix2) = dfdx/2 + dfdxdot/h;
            [r,c,v] = find(dfdx/2 + dfdxdot/h);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = ix2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            %J(iceq,iu1) = dfdu/2;
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = iu1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            %J(iceq,iu2) = dfdu/2;
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = iu2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            
            % advance the indices
            ix1 = ix1 + nvarpernode;
            iu1 = iu1 + nvarpernode;
            iceq = iceq + (ncon_eq/(N-1));
        end
        
        if exist('Jpattern','var')
            [r,c,~] = find(Jpattern);
            nz_rows = ismember([allrows(1:index-1),allcols(1:index-1)],[r,c], 'rows');
            J = sparse(allrows(nz_rows),allcols(nz_rows),allvals(nz_rows),ncon, nvar);
        else
            J = sparse(allrows(1:index-1),allcols(1:index-1),allvals(1:index-1),ncon, nvar);
        end
        
    end

end
