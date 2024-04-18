function Result = das3_optimize(data,filename,OptSetup)
% This program optimizes das3 movements to fit motion capture data
% It uses das3_simplified.osim, which includes DoFs at the thorax, but locks those to the input data
%
%
% Inputs:
%		data: 	The input data (in a table), containing time, angles and forces at the hand
%		filename:	the output file name
%		OptSetup: settings for the optimization:
%					weights for the terms in the objective function:
%						Wdata (kinematic tracking), Weffort (minimization of muscle effort)
%						Wscap (scapulo-thoracic gliding plane) and Whum (glenohumeral stability).
%						Wscap and Whum are optional: if they are not included,
%						scapular and humeral stability are hard constraints
%					N: the number of nodes for Direct Collocation.
%					initialguess, with options: "random", "mid",
% 						and if none of the above, assumed to be name of file with previous solution
%					MaxIter: maximum number of iterations
%					OptimTol and FeasTol: Optimality and feasibility tolerances
%                   Indices in input data and model state for locked and tracking angles,
%                   thoraco-humeral angles, and hand forces, if applicable
% Output:
%		Result:		The output of the optimization, including the settings in OptSetup
%						and the input data. The output is also saved in <filename>.mat
%						(which can be used as a previous solution in initialguess)
%						and <filename>.mot, for visualization in OpenSim

% optimization settings
Result.OptSetup = OptSetup;

% weights for terms in objective function
% kinematic tracking
Wdata = OptSetup.Wdata;
% minimization of muscle effort
Weffort = OptSetup.Weffort;
% scapular stability, if it is a soft constraint
if isfield(OptSetup, 'Wscap')
    Wscap = OptSetup.Wscap;
    Scapcon_flag = 0;
else
    Scapcon_flag = 1;
end
% humeral stability, if it is a soft constraint
if isfield(OptSetup, 'Whum')
    Whum = OptSetup.Whum;
    Humcon_flag = 0;
else
    Humcon_flag = 1;
end

N = OptSetup.N;
initialguess = OptSetup.initialguess;

% Additional model inputs
% Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
ex_mom = zeros(5,1);
% Vertical force of amplitude arm_support(2) applied to the ulna at a distance of arm_support(1) (meters) from the elbow
% to simulate arm support
arm_support = zeros(2,1);

% Load structure with model related variables
modelparams = load(OptSetup.model); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;
ncontrols = 2*nmus;

% get muscle names, SEEslack and lceopt
musnames = cell(nmus,1);
SEEslack = zeros(nmus,1);
LCEopt = zeros(nmus,1);
for imus=1:nmus
    musnames{imus} = model.muscles{imus}.osim_name;
    SEEslack(imus) = model.muscles{imus}.lslack;
    LCEopt(imus) = model.muscles{imus}.lceopt;
end
muscles_table = cell2table(musnames);

% get dof names
dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% some constants
nvarpernode = nstates + ncontrols;                  % number of unknowns per time node
nvar = nvarpernode * N;                             % total number of unknowns
ncon_eq = nstates*(N-1);                            % number of constraints due to discretized dynamics
ncon_eq_elong = nmus*N;                             % number of constraints due to SEE elongation
ncon_neq = Humcon_flag*(N-1) + Scapcon_flag*2*N;    % number of constraints due to glenohumeral and scapular stability
ncon = ncon_eq + ncon_eq_elong + ncon_neq;

% Initialize the model
das3('Initialize',model);

% Range of motion limits
xlims = das3('Limits')';
% allow AC dofs to go below the limit
xlims(8,1) = -0.45;
xlims(9,1) = -0.32;

% Approximate hand position
local_hand_coord = model.joints{13}.location';

% maximum activations/excitations (to simulate injury)
if isfield(OptSetup, 'max_act_table')
    % make sure the muscles are in the correct order
    [~,ii] = ismember(muscles_table.musnames,OptSetup.max_act_table.Var1);
    new_max_act_table = OptSetup.max_act_table(ii,:);
    max_act = new_max_act_table.mc_new;
else
    max_act = ones(nmus,1);
end

% Lower and upper bounds for the optimization variables X
% X will contain, for each node, the states and controls.
L = zeros(nvar,1);
U = zeros(nvar,1);

% Bounds are:
%   joint angles 	xlimdeg
%   angular velocities -40 to 40 rad/s
%   CE lengths		0 to 2
%	active states	0 to 1
%	neural controls 0 to 1
%   normalised SEE elongation: from -0.1 to 0.06 -> max force is
%   ~2.25 times the maximum isometric force

max_elong = 0.06;

for i_node = 0:N-1
    L(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,1)    ;         % q
        (zeros(ndof,1) - 40);                                           % qdot
        zeros(nmus,1) + 0.3;                                            % Lce
        zeros(nmus,1);                                                  % active states
        zeros(nmus,1);                                                  % neural excitations
        -0.01*ones(nmus,1)];                                             % normalised SEE elongation
    
    U(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,2)    ;         % q
        (zeros(ndof,1) + 40);                                           % qdot
        zeros(nmus,1) + 1.7;                                            % Lce
        max_act;                                                        % active states
        max_act;                                                        % neural excitations
        max_elong*ones(nmus,1)];                                        % normalised SEE elongation

end

% Should angular velocities be zero at the start and end of movement?
if OptSetup.start_at_rest
    % First node
    L(ndof+1:2*ndof) = zeros(ndof,1);
    U(ndof+1:2*ndof) = zeros(ndof,1);
end

% Should movement start at initial state?
if OptSetup.initial_state
    init_state = load('initial_state.mat');
    mus_lengths = das3('Musclelengths', init_state.x);
    SEE_elong = (mus_lengths - init_state.x(2*ndof+(1:nmus)).*LCEopt - SEEslack)./SEEslack;
    L(1:nvarpernode) = [init_state.x; init_state.x(2*ndof+nmus+(1:nmus)); SEE_elong];
    U(1:nvarpernode) = L(1:nvarpernode);
end

if OptSetup.end_at_rest
    % Final node
    final_node_start_index = (N-1)*nvarpernode;
    L(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);
    U(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);
end

% Load the input data
Result.input = data(:,OptSetup.input_variables);
Result.input_t = data.time;

% Make the vector of times for direct collocation
t1 = min(data.time);
t2 = max(data.time);
times = t1 + (0:N-1)'*(t2-t1)/(N-1);

% Resample the input data onto the direct collocation times
% and put into a long vector to save time when computing cost function
ntrackdofs = length(OptSetup.tracking_indata);

if ntrackdofs == 0
    data_dofs = [];
    data_dofs_sd = [];
else
    datadofs_m = data(:,OptSetup.tracking_indata);
    datadofs = interp1(data.time,table2array(datadofs_m),times);
    
    Result.resampled_data = datadofs;
    
    datadofs_deriv = diff(datadofs)./diff(times);
    datadofs_deriv = [datadofs_deriv; datadofs_deriv(end,:)];
    
    % convert into a single column
    data_dofs = reshape(datadofs', [], 1);
    data_dofs_deriv = reshape(datadofs_deriv', [], 1);
    
    % define the estimated uncertainty in each measured variable (in radians)
    if ntrackdofs == 11
        data_dofs_sd =1000*[1 1 1000 1 1 1 1 1 1 1 1]';   % for one node 
        data_dofs_sd = repmat(data_dofs_sd,N,1);		% replicate for all nodes
        data_dofs_sd(1:ntrackdofs) = [1 1 1000 1 1 1 1 1 1 1 1]'; % set to lower number to track start of movement better
        data_dofs_sd(ntrackdofs*(N-1)+1:end) = [1 1 1000 1 1 1 1 1 1 1 1]'; % set to lower number to track end of movement better
    else
        data_dofs_sd =1000*ones(ntrackdofs,1);     % for one node 
        data_dofs_sd = repmat(data_dofs_sd,N,1); % replicate for all nodes
        data_dofs_sd(1:ntrackdofs) = ones(ntrackdofs,1); % set to lower number to track start of movement better
        data_dofs_sd(ntrackdofs*(N-1)+1:end) = ones(ntrackdofs,1); % set to lower number to track end of movement better
    end
end

% Do the same for thoraco-humeral angles, if they were provided
if isempty(OptSetup.thorhum_indata)
    thorhum_flag = 0;
    nthorhum = 0;
else
    thorhum_flag = 1;
    nthorhum = length(OptSetup.thorhum_indata);
    datathorhum_m = data(:,OptSetup.thorhum_indata);
    datathorhum = interp1(data.time,table2array(datathorhum_m),times);
    Result.resampled_thorhum = datathorhum;
    data_thorhum = reshape(datathorhum', [], 1);
    
    % define the estimated uncertainty in each measured variable (in radians)
%    data_thorhum_sd = 1*ones(length(OptSetup.thorhum_indata),1);    % for one node (set to higher number to only track start and end of movement)
    data_thorhum_sd = 1000*[1 1 10]';    % for one node (set to higher number to only track start and end of movement)
    data_thorhum_sd = repmat(data_thorhum_sd,N,1);		% replicate for all nodes
    data_thorhum_sd(1:length(OptSetup.thorhum_indata))=1;
    data_thorhum_sd(length(OptSetup.thorhum_indata)*(N-1)+1:end)=1;
    
    % also estimate clavicular and scapular angles, in case we want to use
    % them for an initial guess
    shoulder_angles = zeros(length(times),9);
    for i_t = 1:length(times)
        shoulder_angles(i_t,:) = estimate_shoulder_angles(datathorhum(i_t,1),datathorhum(i_t,2),datathorhum(i_t,3));
    end
    data_clavscaphum = reshape(shoulder_angles', [], 1);    
    
    shoulder_angles_deriv = diff(shoulder_angles)./diff(times);
    shoulder_angles_deriv = [shoulder_angles_deriv; shoulder_angles_deriv(end,:)];
    data_clavscaphum_deriv = reshape(shoulder_angles_deriv', [], 1);    

end

% Do the same for hand position, if provided
if isempty(OptSetup.handpos_indata)
    handpos_flag = 0;
else
    handpos_flag = 1;
    datahandpos_m = data(:,OptSetup.handpos_indata);
    datahandpos = interp1(data.time,table2array(datahandpos_m),times);
    data_handpos = reshape(datahandpos', [], 1);
    
    % define the estimated uncertainty in each measured variable (in radians)
    data_handpos_sd = 1*ones(length(OptSetup.handpos_indata),1);    % for one node (set to higher number to only track start and end of movement)
    data_handpos_sd = repmat(data_handpos_sd,N,1);		% replicate for all nodes
    data_handpos_sd(1:length(OptSetup.handpos_indata))=1;
    data_handpos_sd(length(OptSetup.handpos_indata)*(N-1)+1:end)=1;
end


% precompute the indices for activations and some angles, so we can compute cost function quickly
iact = zeros(1,N*nmus);
iu = zeros(1,N*nmus);
iSEEelong = zeros(1,N*nmus);
iLce = 2*ndof + (1:nmus);

imeas = zeros(1,N*ntrackdofs);
idmeas = zeros(1,N*ntrackdofs);

ishoulder = zeros(1,N*9);
idshoulder = zeros(1,N*9);

nlockeddofs = length(OptSetup.lockeddof_indata);
ilocked_dof = zeros(1,N*nlockeddofs);
ilocked_ddof = zeros(1,N*nlockeddofs);

for i_node=0:N-1
    iact(i_node*nmus+(1:nmus)) = nvarpernode*i_node+2*ndof+nmus+(1:nmus);
    iu(i_node*nmus+(1:nmus)) = nvarpernode*i_node+ nstates + (1:nmus);
    iSEEelong(i_node*nmus+(1:nmus)) = nvarpernode*i_node+ nstates + nmus +(1:nmus);
    
    imeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+OptSetup.tracking_inx;
    idmeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+ndof+OptSetup.tracking_inx;
    
    ishoulder(i_node*9+(1:9)) = nvarpernode*i_node+3+(1:9);
    idshoulder(i_node*9+(1:9)) = nvarpernode*i_node+ndof+3+(1:9);
    
    ilocked_dof(i_node*nlockeddofs+(1:nlockeddofs)) = nvarpernode*i_node+OptSetup.lockeddof_inx;
    ilocked_ddof(i_node*nlockeddofs+(1:nlockeddofs)) = nvarpernode*i_node+ndof+OptSetup.lockeddof_inx;
end

% Locked dofs stay at measured values, and have zero velocity
datalockeddofs_m = data(:,OptSetup.lockeddof_indata);
datalockeddofs = interp1(data.time,table2array(datalockeddofs_m),times);

Result.resampled_lockeddata = datalockeddofs;

datalockeddofs_deriv = diff(datalockeddofs)./diff(times);
datalockeddofs_deriv = [datalockeddofs_deriv; datalockeddofs_deriv(end,:)];

% convert into a single column
data_lockeddofs = reshape(datalockeddofs', [], 1);
data_lockeddofs_deriv = reshape(datalockeddofs_deriv', [], 1);

U(ilocked_dof) = data_lockeddofs;
L(ilocked_dof) = data_lockeddofs;
U(ilocked_ddof) = data_lockeddofs_deriv;
L(ilocked_ddof) = data_lockeddofs_deriv;

% Precompute the indices for muscles that cross GH
iGH = [];
for imus=1:nmus
    if(das3('crossGH', imus))
        iGH=[iGH imus];
    end
end

% For glenohumeral constraint:
aphi=tand(38.55);
atheta=tand(44.37);
Rgt = glenoid_scap;

% Put hand forces in separate matrix
hand_force = zeros(length(data.time),3);
hand_force(:,OptSetup.handF_inF) = table2array(data(:,OptSetup.handF_indata));
hand_force = interp1(data.time,hand_force,times);
hand_force = reshape(hand_force', [], 1);

%% Make an initial guess

% Initial guess for kinematics + velocities is the interpolated measured data
% For activations and fibre lengths there are different options: mid and random

% If we are tracking thoraco-humeral angles, use the "normal"
% scapulohumeral rhythm as the initial guess

% Alternatively, we load an earlier solution

if strcmp(initialguess, 'mid')
    X0 = (L + U)/2;						% halfway between upper and lower bound
    X0 = X0 + 0.001;					% to avoid exact zero velocity, for which Autolev equations do not work
    X0(imeas) = data_dofs;
    X0(idmeas) = data_dofs_deriv;
    if thorhum_flag
        X0(ishoulder(4:end)) = data_clavscaphum(4:end);
        X0(idshoulder(4:end)) = data_clavscaphum_deriv(4:end);        
    end
    ix = 1:nstates;
    for i_node=1:N
        i_node_states = X0(ix);
        i_node_lengths = das3('Musclelengths',i_node_states);
        X0(ix(end)+nmus+(1:nmus)) = 0; % SEE elongation
        X0(ix(2*ndof + (1:nmus))) = (i_node_lengths - SEEslack)./LCEopt; % normalised fibre length
        ix = ix + nvarpernode;
    end

 elseif numel(strfind(initialguess, 'init_dsem')) > 0
    x0 = zeros(nstates,1);	
    x0(4:14) = [-0.3802;0.11;0;0.808;0.0855;-0.0206;0;-0.0873;0;0.0873;1.5708]; % initial DSEM position
    u0 = zeros(nmus,1);    
    lengths0 = das3('Musclelengths',x0);
    elong0 = zeros(nmus,1); % SEE elongation
    x0(iLce) = (lengths0 - SEEslack)./LCEopt;
    
    x_0 = repmat(x0',length(times),1);
    u_0 = repmat(u0',length(times),1);
    elong_0 = repmat(elong0',length(times),1);
    X0 = reshape([x_0 u_0 elong_0]',nvar,1);

elseif numel(strfind(initialguess, 'random')) > 0
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
    X0(iact) = 0;
    X0(iu) = 0;
    X0(imeas) = data_dofs;
    X0(idmeas) = data_dofs_deriv;
    if thorhum_flag
        X0(ishoulder) = data_clavscaphum;
        X0(idshoulder) = data_clavscaphum_deriv;        
    end
    ix = 1:nstates;
    for i_node=1:N
        i_node_states = X0(ix);
        i_node_lengths = das3('Musclelengths',i_node_states);
        X0(ix(end)+nmus+(1:nmus)) = 0;% SEE elongation
        X0(ix(2*ndof + (1:nmus))) = (i_node_lengths - SEEslack)./LCEopt; % normalised fibre length
        ix = ix + nvarpernode;
    end

elseif numel(strfind(initialguess, 'equilibrium')) > 0
    init_state = load('initial_state.mat');
    x0 = init_state.x';
    u0 = init_state.x(2*ndof+nmus+(1:nmus))';    
    mus_lengths = das3('Musclelengths', init_state.x);
    elong0 = ((mus_lengths - init_state.x(2*ndof+(1:nmus)).*LCEopt - SEEslack)./SEEslack)';    
    x_0 = repmat(x0,length(times),1);
    u_0 = repmat(u0,length(times),1);
    elong_0 = repmat(elong0,length(times),1);
    X0 = reshape([x_0 u_0 elong_0]',nvar,1);

else
    % load a previous solution, initialguess contains file name
    initg = load(initialguess);
    t0 = initg.Result.times;
    x0 = initg.Result.x';
    u0 = initg.Result.u';
    if isfield(initg.Result,'elong')
        elong0 = initg.Result.elong';
    else
        elong0 = zeros(1,nmus);
    end
    
    % interpolate states and controls from initial guess to the current time grid
    if length(t0)>1
        x_0 = interp1(t0,x0,times,'linear','extrap');
        u_0 = interp1(t0,u0,times,'linear','extrap');
        elong_0 = interp1(t0,elong0,times,'linear','extrap');
   else
        x_0 = repmat(x0,length(times),1);
        u_0 = repmat(u0,length(times),1);
        elong_0 = repmat(elong0,length(times),1);
    end
    X0 = reshape([x_0 u_0 elong_0]',nvar,1);
end

tic

%% Run optimization
if (OptSetup.MaxIter > 0)
    
    print_flag = 0;
    
    % determine sparsity structure of Jacobian (midpoint discretization)
    
    Jnnz = 1;
    X = L + (U-L).*rand(size(L));		% a random vector of unknowns
    J = conjac(X);
    Jnnz = nnz(J);
    fprintf('Jacobian sparsity: %d nonzero elements out of %d (%5.3f%%).\n',Jnnz, ncon*nvar, 100*Jnnz/(ncon*nvar));
    Jpattern = double(J~=0);
    fprintf('Nonzero elements in Jpattern: %d\n', nnz(Jpattern));
    %figure; spy(Jpattern);
    
    % check the derivatives?
    checkderivatives = 0;
    
    if (checkderivatives)
        hh = 1e-7;
        X = L + (U - L).*rand(size(L));
        %X = X0;
        objf = objfun(X);
        grad = objgrad(X);
        cf = confun(X);
        cjac = conjac(X);
        cjac_num = zeros(ncon,nvar);
        grad_num = zeros(nvar,1);
        for i_var=1:nvar
            fprintf('checking derivatives for unknown %4d of %4d\n',i_var,nvar);
            Xisave = X(i_var);
            X(i_var) = X(i_var) + hh;
            cfi = confun(X);
            cjac_num(:,i_var) = (cfi - cf)/hh;
            grad_num(i_var) =   (objfun(X) - objf)/hh;
            X(i_var) = Xisave;
        end
        
        % find the max difference in constraint jacobian and objective gradient
        figure; spy(cjac_num);
        cjac_full = full(cjac);
        figure; spy(cjac_full);
        [maxerr,irow] = max(abs(cjac-cjac_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max error in constraint jacobian: %8.5f at %d %d\n', maxerr, irow(icol), icol);
        [maxerr,irow] = max(abs(grad_num-grad));
        fprintf('Max error in objective gradient:  %8.5f at %d\n', maxerr, irow);
        keyboard
    end
    
    % evaluate initial guess
    print_flag=1;
    objfun(X0);
    confun(X0);
    print_flag=0;
    
    % set up constraints
    cl = zeros(ncon_eq+ncon_eq_elong,1);
    if Humcon_flag, cl = [cl; zeros(N-1,1)]; end
    if Scapcon_flag, cl = [cl; -0.19*ones(2*N,1)]; end
    cu = zeros(ncon_eq+ncon_eq_elong,1);
    if Humcon_flag, cu = [cu; ones(N-1,1)]; end
    if Scapcon_flag, cu = [cu; 0.21*ones(2*N,1)]; end
        
    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;
    
    options.cl = cl;
    options.cu = cu;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    
    options.lb = L;
    options.ub = U;
    
    % IPOPT options
    options.ipopt.print_level = 0;
    options.ipopt.max_iter = OptSetup.MaxIter;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.tol = OptSetup.OptimTol;
    options.ipopt.acceptable_constr_viol_tol = OptSetup.FeasTol;
    options.ipopt.acceptable_tol = OptSetup.FeasTol;
    options.ipopt.constr_viol_tol = OptSetup.FeasTol;
    options.warm_start_init_point = 'yes';
    options.warm_start_bound_push = 1e-12;
    options.warm_start_bound_frac = 1e-12;
    options.warm_start_slack_bound_frac = 1e-12;
    options.warm_start_slack_bound_push = 1e-12;
    options.warm_start_mult_bound_push = 1e-12;
    options.warm_start_same_structure = 'yes';
    options.ipopt.print_info_string = 'yes';
    options.ipopt.limited_memory_max_history = 12;
    options.ipopt.dual_inf_tol = 1;
    
    [X, info] = ipopt(X0,funcs,options);
    fprintf('IPOPT info status: %d\n', info.status);
    Result.status = info.status;
    Result.info = info;
    Result.X = X;
    Result.X0 = X0;
    
else		% skip optimization
    Result.X = X0;
end

Result.duration = toc;
disp(['Time elapsted: ' num2str(Result.duration) ' seconds']);
beep;

obj_str = '';
print_flag=1;
objfun(Result.X);
confun(Result.X);
print_flag=0;

%% Save this result in a mat file and an Opensim motion file

x = reshape(Result.X,nvarpernode,N);
Result.ndof = ndof;
Result.nmus = nmus;
Result.muscles = model.muscles;
Result.times = times;
Result.input_data = data;
Result.x = x(1:nstates,:);
Result.obj_str = obj_str;
Result.u = x(nstates+(1:nmus),:);
Result.elong = x(nstates+nmus+(1:nmus),:);
angles = x(1:ndof,:);

% Calculate muscle lengths and forces
Result.mus_lengths = zeros(size(Result.u));
Result.mus_forces = zeros(size(Result.u));
Result.thor_hum = zeros(3,N);
Result.hand_pos = zeros(3,N);
Result.Fscap = zeros(2,N);
Result.FGHcontact = zeros(1,N);

for i_node=1:Result.OptSetup.N
    Result.mus_lengths(:,i_node) = das3('Musclelengths', Result.x(:,i_node));
    Result.mus_forces(:,i_node) = das3('Muscleforces', Result.x(:,i_node));
    [~, ~, ~, ~, ~, ~, Result.thor_hum(:,i_node)] = das3('Dynamics', Result.x(:,i_node), 0*Result.x(:,i_node),zeros(nmus,1));
    Result.Fscap(:,i_node) = das3('Scapulacontact', Result.x(:,i_node))';
    Result.hand_pos(:,i_node) = das3_handpos(Result.x(:,i_node),local_hand_coord);
end

ixs = 1:nstates;
ius = nstates+(1:nmus);
% evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
for i=1:N-1
    x1 = Result.X(ixs);
    x2 = Result.X(ixs+nvarpernode);
    u1 = Result.X(ius);
    u2 = Result.X(ius+nvarpernode);
    h = times(i+1) - times(i);		% time interval

    handF1 = hand_force(i*3-2:i*3);
    handF2 = hand_force((i+1)*3-2:(i+1)*3);
    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);

    Result.FGHcontact(i) = calculate_FGH(FGH);

    % advance the indices
    ixs = ixs + nvarpernode;
    ius = ius + nvarpernode;
end

Result.FGHcontact(N) = Result.FGHcontact(N-1);

save([filename '.mat'],'Result');
make_osimm(filename,dofnames,angles,times);

%% Plot elbow angle and velocity, and muscle excitations and forces
% figure;
% ndof = Result.ndof;
% nmus = Result.nmus;
% subplot(6,1,1); plot(Result.times,Result.x(13:14,:)'*180/pi,'o-'); ylabel('Angle (degrees)');
% title(Result.obj_str)
% subplot(6,1,2); plot(Result.times,Result.x(ndof+(13:14),:)'*180/pi,'o-'); ylabel('Velocity (degrees/second)'); legend('Flexion/extension','Pronation/supination');
% subplot(6,1,3); plot(Result.times,Result.u,'o-'); ylabel('Excitations');
% subplot(6,1,4); plot(Result.times,Result.x(2*ndof+nmus+1:2*ndof+2*nmus,:),'o-'); ylabel('Activations');
% subplot(6,1,5); plot(Result.times,Result.x(2*ndof+1:2*ndof+nmus,:),'o-'); ylabel('Norm fibre lengths');
% subplot(6,1,6); plot(Result.times,Result.mus_forces','o-'); ylabel('Forces (N)');
% xlabel('time (s)');


%% Functions to specify objective and constraints (and their derivatives)

% Objective function
    function f = objfun(X)
        
        % First term of cost function is mean of squared differences to
        % measured data, including thoraco-humeral data if provided
        if ntrackdofs
            wf1 = mean(((X(imeas)-data_dofs)./data_dofs_sd).^2);
        else
            wf1 = 0;
        end
        
        % If required, calculate thoracohumeral angles
        if thorhum_flag
            thorhum = zeros(nthorhum*N,1);
            
            ixs = 1:nstates;
            for inode=1:N
                [~, ~, ~, ~, ~, ~, thorhum_x] = das3('Dynamics',X(ixs),zeros(size(X(ixs))),zeros(nmus,1));
                thorhum((inode-1)*nthorhum+(1:nthorhum)) = thorhum_x(OptSetup.thorhum_inx);
                ixs = ixs + nvarpernode;
            end
            wf1 = wf1 + mean(((thorhum-data_thorhum)./data_thorhum_sd).^2);
        end
        
        % If required, calculate hand position
        if handpos_flag
            handpos = zeros(3*N,1);
            
            ixs = 1:nstates;
            for inode=1:N
                handpos((inode-1)*3+(1:3)) = das3_handpos(X(ixs),local_hand_coord);
                ixs = ixs + nvarpernode;
            end
            wf1 = wf1 + mean(((handpos-data_handpos)./data_handpos_sd).^2);
        end
        f1 = Wdata * wf1;

        % Second term is mean squared muscle activation        
        wf2 = mean(X(iact).^2);
        f2 = Weffort * wf2;
                
        f = f1 + f2;
        
        if print_flag
            obj_str = sprintf('Objfun (weighted): %9.5f = %9.5f (fit) + %9.5f (effort) ', f,f1,f2);
            wobj_str = sprintf('Objfun (unweighted): %9.5f (fit) + %9.5f (effort) ', wf1,wf2);
        end
        
        % Third term is thorax-scapula constraint
        if ~Scapcon_flag && Wscap
            Fscap = zeros(N,2);            
            ixs = 1:nstates;
            for inode=1:N
                Fscap(inode,:) = das3('Scapulacontact', X(ixs));
                ixs = ixs + nvarpernode;
            end            
            wf3 = mean(Fscap(:,1).^2)+mean(Fscap(:,2).^2);
            f3 = Wscap * wf3;
            
            f = f + f3;
            
            if print_flag
                obj_str3 = sprintf('+ %9.5f (scapula)', f3);
                wobj_str3 = sprintf('+ %9.5f (scapula)', wf3);
                obj_str = strcat(obj_str,obj_str3);
                wobj_str = strcat(wobj_str,wobj_str3);
            end
            Result.Fscap = Fscap;
        end
        
        % Fourth term is glenohumeral stability constraint
        if ~Humcon_flag && Whum
            
            FGHcontact = zeros(N-1,1);
            
            % indices of states and controls of node 1
            ixs = 1:nstates;
            ius = nstates+(1:nmus);
            % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
            for inode=1:N-1
                x1 = X(ixs);
                x2 = X(ixs+nvarpernode);
                u1 = X(ius);
                u2 = X(ius+nvarpernode);
                h = times(inode+1) - times(inode);		% time interval
                
                handF1 = hand_force(inode*3-2:inode*3);
                handF2 = hand_force((inode+1)*3-2:(inode+1)*3);
                [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                
                FGHcontact(inode) = calculate_FGH(FGH);
                
                % advance the indices
                ixs = ixs + nvarpernode;
                ius = ius + nvarpernode;
            end
            
            wf4 = mean(FGHcontact.^2);
            f4 = Whum * wf4;
            f = f + f4;
            if print_flag
                obj_str4 = sprintf('+ %9.5f (humerus)', f4);
                wobj_str4 = sprintf('+ %9.5f (humerus)', wf4);
                obj_str = strcat(obj_str,obj_str4);
                wobj_str = strcat(wobj_str,wobj_str4);
            end
            
            Result.FGHcontact = FGHcontact;
        end
        
        if print_flag
            fprintf(obj_str);
            fprintf('\n');
            fprintf(wobj_str);
            fprintf('\n');
        end
        
    end


% Gradient of objective function
    function g = objgrad(X)

        g = zeros(nvar,1);
        
        % Perturbation to estimate derivative
        dx = 1e-7;
        
        % First term of cost function is mean of squared differences to
        % measured data (thoraco-humeral angles are calculated separate,
        % below)
        if ntrackdofs
            g(imeas) = g(imeas) + 2*Wdata*(X(imeas)-data_dofs)./data_dofs_sd.^2/length(imeas);
        end
        
        % If required, calculate scapulothoracic term derivatives
        if thorhum_flag
            ixs = 1:nstates;
            ithorhum = 1:nthorhum;
            for inode=1:N
                [~, ~, ~, ~, ~, ~, thorhum] = das3('Dynamics',X(ixs),zeros(size(X(ixs))),zeros(nmus,1));
                thorhum = thorhum(OptSetup.thorhum_inx);
                
                % thoraco-humeral angles only depend on dofs 4-12
                % (SC_y,SC_z,SC_x,AC_y,AC_z,AC_x,GH_y,GH_z,GH_yy)
                for istate = 4:12
                    xisave = X(ixs(istate));
                    X(ixs(istate)) = X(ixs(istate)) + dx;
                    
                    [~, ~, ~, ~, ~, ~, thorhum_dx] = das3('Dynamics',X(ixs),zeros(size(X(ixs))),zeros(nmus,1));
                    thorhum_dx = thorhum_dx(OptSetup.thorhum_inx);
                    dthorhum = (thorhum_dx-thorhum)/dx;
                    g(ixs(istate)) = g(ixs(istate)) + ...
                        sum(((2*Wdata*(thorhum-data_thorhum(ithorhum)))./data_thorhum_sd(ithorhum).^2).*dthorhum/(N*nthorhum));
                    X(ixs(istate)) = xisave;
                end
                ixs = ixs + nvarpernode;
                ithorhum = ithorhum + nthorhum;
            end
        end
 
         % If required, calculate hand position term derivatives
        if handpos_flag
            ixs = 1:nstates;
            ihandpos = 1:3;
            for inode=1:N
                handpos = das3_handpos(X(ixs),local_hand_coord);
                
                % thoraco-humeral angles only depend on dofs 1-14
                for istate = 1:14
                    xisave = X(ixs(istate));
                    X(ixs(istate)) = X(ixs(istate)) + dx;
                    
                    handpos_dx = das3_handpos(X(ixs),local_hand_coord);
                    dhandpos = (handpos_dx-handpos)/dx;
                    g(ixs(istate)) = g(ixs(istate)) + ...
                        sum(((2*Wdata*(handpos-data_handpos(ihandpos)))./data_handpos_sd(ihandpos).^2).*dhandpos/(N*3));
                    X(ixs(istate)) = xisave;
                end
                ixs = ixs + nvarpernode;
                ihandpos = ihandpos + 3;
            end
        end
       
        
        % Second term is mean squared muscle activation
        g(iact) = g(iact) + 2*Weffort*(X(iact))./(N*nmus);
                
        % If required, calculate scapula term derivatives
        if ~Scapcon_flag && Wscap
            
            ixs = 1:nstates;
            for inode=1:N
                Fscap = das3('Scapulacontact', X(ixs));
                
                % Scapulo-thoracic glinding plane only depends on dofs 4-9
                % (SC_y,SC_z,SC_x,AC_y,AC_z,AC_x)
                for istate = 4:9
                    xisave = X(ixs(istate));
                    X(ixs(istate)) = X(ixs(istate)) + dx;
                    Fscap_dx = das3('Scapulacontact', X(ixs));
                    dFscap = (Fscap_dx-Fscap)/dx;
                    g(ixs(istate)) = g(ixs(istate)) + Wscap*(2*Fscap(1)*dFscap(1)/N + 2*Fscap(2)*dFscap(2)/N);
                    X(ixs(istate)) = xisave;
                end
                ixs = ixs + nvarpernode;
            end
        end
        
        % For humeral stability we need to evaluate dynamics
        if ~Humcon_flag && Whum
            ix1 = 1:nstates;
            iu1 = nstates+(1:nmus);
            for inode=1:N-1
                ix2 = ix1 + nvarpernode;
                iu2 = iu1 + nvarpernode;
                x1 = X(ix1);
                x2 = X(ix2);
                u1 = X(iu1);
                u2 = X(iu2);
                h = times(inode+1) - times(inode);		% time interval
                
                handF1 = hand_force(inode*3-2:inode*3);
                handF2 = hand_force((inode+1)*3-2:(inode+1)*3);
                
                [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                
                % calculate GH constraint
                FGHcontact = calculate_FGH(FGH);
                
                % Perturb x1 and x2 by dx to estimate derivative
                for istate = [1:2*ndof 2*ndof+iGH 2*ndof+nmus+iGH]
                    xisave = x1(istate);
                    x1(istate) = x1(istate) + dx;
                    
                   % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH);
                    dFGHcontact = (FGHcontact_dx-FGHcontact)/dx;
                    g(ix1(istate)) = g(ix1(istate)) + Whum * 2 * FGHcontact * dFGHcontact / (N-1);
                    x1(istate) = xisave;
                    
                    xisave = x2(istate);
                    x2(istate) = x2(istate) + dx;
                    
                    % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH);
                    dFGHcontact = (FGHcontact_dx-FGHcontact)/dx;
                    g(ix2(istate)) = g(ix2(istate)) + Whum * 2 * FGHcontact * dFGHcontact / (N-1);
                    x2(istate) = xisave;
                end
                
                ix1 = ix1 + nvarpernode;
                iu1 = iu1 + nvarpernode;
            end
        end
        
    end


% Constraints
    function allcon = confun(X)

        c = zeros(ncon_neq,1);
        c_flag = 0*c;
        ceq = zeros(ncon_eq+ncon_eq_elong,1);
        
        % Calculate constraints due to discretized dynamics
        % indices of states and controls of node 1
        ixs = 1:nstates;
        ius = nstates+(1:nmus);
        ielong = nstates+nmus+(1:nmus);
        % indices of equality constraints at node 1
        iceq_f = 1:(ncon_eq/(N-1));
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for inode=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            u1 = X(ius);
            u2 = X(ius+nvarpernode);
            h = times(inode+1) - times(inode);		% time interval

            handF1 = hand_force(inode*3-2:inode*3);
            handF2 = hand_force((inode+1)*3-2:(inode+1)*3);
            
            [f, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);

            f(OptSetup.lockeddof_inx) = zeros(nlockeddofs,1);
            f(ndof+OptSetup.lockeddof_inx) = zeros(nlockeddofs,1);
            ceq(iceq_f) = f;
            
            if Humcon_flag 
                c(inode) = calculate_FGH(FGH);
                c_flag(inode) = c(inode)>0 && c(inode)<1;
            end
            
            % advance the indices
            ixs = ixs + nvarpernode;
            ius = ius + nvarpernode;
            iceq_f = iceq_f + (ncon_eq/(N-1));
        end
        
        % calculate SEE elongation at each node
        ixs = 1:nstates;
        iceq_elong = ncon_eq+(1:(ncon_eq_elong/N));
        for inode = 1:N
            x1 = X(ixs);
            elong = X(ielong);
            lengths = das3('Musclelengths',x1);
            ceq(iceq_elong) = elong - (lengths - x1(iLce).*LCEopt - SEEslack)./SEEslack;
            
            % advance the indices
            ixs = ixs + nvarpernode;
            ielong = ielong + nvarpernode;
            iceq_elong = iceq_elong + (ncon_eq_elong/N);
        end
                    
        % Now calculate constraints due to scapular stability at each node
        if Scapcon_flag
            % indices of states and constraints of node 1
            ixs = 1:nstates;
            ic = Humcon_flag*(N-1) + (1:2);
            for inode=1:N
                % solve thorax ellipsoid surface equation for TS and AI, to find
                % out whether they are inside, on, or outside the thorax
                c(ic) = (das3('Scapulacontact', X(ixs)))';
                c_flag(ic) = c(ic)>-0.19 & c(ic)<0.21;
                
                % advance the indices
                ixs = ixs + nvarpernode;
                ic = ic + 2;
            end
        end
        
        allcon = [ceq; c];
        
        if (print_flag)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq(1:ncon_eq)));
            fprintf('Norm(ceq SEE elong): %9.5f  \n', norm(ceq(ncon_eq+(1:ncon_eq_elong))));
            if Scapcon_flag || Humcon_flag
                fprintf('Proportion of satisfied inequality constraints: %9.5f \n', sum(c_flag)/length(c_flag));
            end
        end
        
    end

% Structure of Jacobian matrix
    function J = conjacstructure()
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
        iu1 = nstates+(1:nmus);
        ielong = nstates+nmus+(1:nmus);
        % indices of equality constraints at node 1
        iceq_f = 1:(ncon_eq/(N-1));
        % index of non-equality glenohumeral constraint at node 1
        ic = (ncon_eq+ncon_eq_elong)+1;
        dx = 1e-7;
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for inode=1:N-1
            ix2 = ix1 + nvarpernode;
            iu2 = iu1 + nvarpernode;
            x1 = X(ix1);
            x2 = X(ix2);
            u1 = X(iu1);
            u2 = X(iu2);
            h = times(inode+1) - times(inode);		% time interval
            
            handF1 = hand_force(inode*3-2:inode*3);
            handF2 = hand_force((inode+1)*3-2:(inode+1)*3);
            
            [~, dfdx, dfdxdot, dfdu, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                       
            dfdx(OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nstates);
            dfdx(ndof+OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nstates);
            dfdxdot(OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nstates);
            dfdxdot(ndof+OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nstates);
            dfdu(OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nmus);
            dfdu(ndof+OptSetup.lockeddof_inx,:) = zeros(nlockeddofs,nmus);
            
            % which generates four blocks in the Jacobian:
            [r,c,v] = find(dfdx/2 - dfdxdot/h);
            datal = length(v);
            allrows(index:index+datal-1) = iceq_f(1)+r-1;
            allcols(index:index+datal-1) = ix1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdx/2 + dfdxdot/h);
            datal = length(v);
            allrows(index:index+datal-1) = iceq_f(1)+r-1;
            allcols(index:index+datal-1) = ix2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq_f(1)+r-1;
            allcols(index:index+datal-1) = iu1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq_f(1)+r-1;
            allcols(index:index+datal-1) = iu2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
                                   
            if Humcon_flag
                % Glenohumeral stability constraint
                FGHcontact = calculate_FGH(FGH);
                
                % Perturb x1 and x2 by dx to estimate derivative
                for istate = [1:2*ndof 2*ndof+iGH 2*ndof+nmus+iGH]
                    xisave = x1(istate);
                    x1(istate) = x1(istate) + dx;
                    
                    % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH);
                    
                    allrows(index) = ic;
                    allcols(index) = ix1(1)+istate-1;
                    allvals(index) = (FGHcontact_dx-FGHcontact)/dx;
                    index=index+1;
                    x1(istate) = xisave;
                    
                    xisave = x2(istate);
                    x2(istate) = x2(istate) + dx;
                    
                    % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH);
                    
                    allrows(index) = ic;
                    allcols(index) = ix2(1)+istate-1;
                    allvals(index) = (FGHcontact_dx-FGHcontact)/dx;
                    index=index+1;
                    x2(istate) = xisave;
                end
                ic = ic + 1;
            end
            
            % advance the indices
            ix1 = ix1 + nvarpernode;
            iu1 = iu1 + nvarpernode;
            iceq_f = iceq_f + (ncon_eq/(N-1));
        end
        
        % SEE elongation
        ix1 = 1:nstates;
        iceq_elong = ncon_eq+(1:(ncon_eq_elong/N));

        for inode=1:N
            x1 = X(ix1);
            elong = X(ielong);
            lengths = das3('Musclelengths',x1);
            elong_diff = elong - (lengths - x1(iLce).*LCEopt - SEEslack)./SEEslack;
 
            % Perturb x1 to estimate derivative
            for istate = 1:nstates
                xisave = x1(istate);
                x1(istate) = x1(istate) + dx;

                % estimate delta(elong_diff)/delta(x)
                lengths = das3('Musclelengths',x1);
                elong_diff_dx = elong - (lengths - x1(iLce).*LCEopt - SEEslack)./SEEslack;

                allrows(index:index+nmus-1) = iceq_elong;
                allcols(index:index+nmus-1) = ix1(1)+istate-1;
                allvals(index:index+nmus-1) = (elong_diff_dx-elong_diff)/dx;
                index=index+nmus;
                x1(istate) = xisave;
            end
            
            allrows(index:index+nmus-1) = iceq_elong;
            allcols(index:index+nmus-1) = ielong;
            allvals(index:index+nmus-1) = ones(nmus,1);
            index=index+nmus;
                       
            % advance the indices
            ix1 = ix1 + nvarpernode;
            ielong = ielong + nvarpernode;
            iceq_elong = iceq_elong + (ncon_eq_elong)/N;
       end
        
        if Scapcon_flag
            % indices of states of node 1
            ixs = 1:nstates;
            ic = ncon_eq + ncon_eq_elong + Humcon_flag*(N-1) + (1:2);
            dx = 1e-7;
            
            for inode=1:N
                % Scapular stability constraint
                Fscap = das3('Scapulacontact', X(ixs));
                
                % Fscap only depends on the clavicular and scapular dofs (4-9):
                for ivar = 4:9
                    xisave = X(ixs(ivar));
                    X(ixs(ivar)) = X(ixs(ivar)) + dx;
                    Fscap_dx = das3('Scapulacontact', X(ixs));
                    
                    allrows(index:index+1) = ic';
                    allcols(index:index+1) = [ixs(ivar);ixs(ivar)];
                    allvals(index:index+1) = (Fscap_dx-Fscap)/dx;
                    index=index+2;
                    
                    X(ixs(ivar)) = xisave;
                end
                
                ixs = ixs + nvarpernode;
                ic = ic + 2;
            end
        end
        
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

% Function to calculate location of hand in global coordinate frame
    function handpos = das3_handpos(x,local_hand_coord)

        % get the 3D bone position and orientation 
        d = das3('Visualization',x);
        hand_coord = d(13,1:3)';		% position vector of hand
        R = reshape(d(13,4:12),3,3)';	% orientation matrix

        % transform to global coordinates
        handpos = hand_coord + R*local_hand_coord;
    end


end % end of function das3_optimize
