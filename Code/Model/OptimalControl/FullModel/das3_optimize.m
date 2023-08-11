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
%					initialguess, with options: "init" (the passive equilibrium), "random","mid",
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
ncontrols = nmus;

% get muscle names
musnames = cell(nmus,1);
for imus=1:nmus
    musnames{imus} = model.muscles{imus}.osim_name;
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
ncon_neq = Humcon_flag*(N-1) + Scapcon_flag*2*N;    % number of constraints due to glenohumeral and scapular stability
ncon = ncon_eq + ncon_neq;

% Initialize the model
das3('Initialize',model);

% Range of motion limits
xlims = das3('Limits')';

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
for i_node = 0:N-1
    L(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,1)-0.1;         % q
        (zeros(ndof,1) - 40);                                           % qdot
        zeros(nmus,1) + 0.3;                                            % Lce
        zeros(nmus,1);                                                  % active states
        zeros(nmus,1) ];                                                % neural excitations
    
    U(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,2)+0.1;         % q
        (zeros(ndof,1) + 40);                                           % qdot
        zeros(nmus,1) + 1.7;                                            % Lce
        max_act;                                                        % active states
        max_act];                                                       % neural excitations
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

% Load the input data
Result.input = data(:,OptSetup.input_variables);
Result.input_t = data.time;

% Make the vector of times for direct collocation
t1 = min(data.time);
t2 = max(data.time);
times = t1 + (0:N-1)'*(t2-t1)/(N-1);

% Resample the input data onto the direct collocation times
% and put into a long vector to save time when computing cost function
datadofs_m = data(:,OptSetup.tracking_indata);
datadofs = interp1(data.time,table2array(datadofs_m),times);

Result.resampled_data = datadofs;

datadofs_deriv = diff(datadofs)./diff(times);
datadofs_deriv = [datadofs_deriv; datadofs_deriv(end,:)];

% convert into a single column
data_dofs = reshape(datadofs', [], 1);
data_dofs_deriv = reshape(datadofs_deriv', [], 1);

% define the estimated uncertainty in each measured variable (in radians)
ntrackdofs = length(OptSetup.tracking_indata);
data_dofs_sd = 1*ones(ntrackdofs,1);              % for one node (set to higher number to only track start and end of movement)
data_dofs_sd = repmat(data_dofs_sd,N,1);		% replicate for all nodes
data_dofs_sd(1:ntrackdofs) = 1;
data_dofs_sd(ntrackdofs*(N-1)+1:end) = 1;

% Do the same for thoraco-humeral angles, if they were provided
if isempty(OptSetup.thorhum_indata)
    thorhum_flag = 0;
    nthorhum = 0;
else
    thorhum_flag = 1;
    nthorhum = length(OptSetup.thorhum_indata);
    datathorhum_m = data(:,OptSetup.thorhum_indata);
    datathorhum = interp1(data.time,table2array(datathorhum_m),times);
    data_thorhum = reshape(datathorhum', [], 1);
    
    % define the estimated uncertainty in each measured variable (in radians)
    data_thorhum_sd = 1*ones(length(OptSetup.thorhum_indata),1);    % for one node (set to higher number to only track start and end of movement)
    data_thorhum_sd = repmat(data_thorhum_sd,N,1);		% replicate for all nodes
    data_thorhum_sd(1:length(OptSetup.thorhum_indata))=1;
    data_thorhum_sd(length(OptSetup.thorhum_indata)*(N-1)+1:end)=1;
end

% precompute the indices for activations and some angles, so we can compute cost function quickly
iact = zeros(1,N*nmus);
iLce = zeros(1,N*nmus);

imeas = zeros(1,N*ntrackdofs);
idmeas = zeros(1,N*ntrackdofs);

nlockeddofs = length(OptSetup.lockeddof_indata);
ilocked_dof = zeros(1,N*nlockeddofs);
ilocked_ddof = zeros(1,N*nlockeddofs);

for i_node=0:N-1
    iLce(i_node*nmus+(1:nmus)) = nvarpernode*i_node+2*ndof+(1:nmus);
    iact(i_node*nmus+(1:nmus)) = nvarpernode*i_node+2*ndof+nmus+(1:nmus);
    
    imeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+OptSetup.tracking_inx;
    idmeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+ndof+OptSetup.tracking_inx;
    
    ilocked_dof(i_node*nlockeddofs+(1:nlockeddofs)) = nvarpernode*i_node+OptSetup.lockeddof_inx;
    ilocked_ddof(i_node*nlockeddofs+(1:nlockeddofs)) = nvarpernode*i_node+ndof+OptSetup.lockeddof_inx;
end

% Locked dofs stay at measured values, and have zero velocity
U(ilocked_dof) = repmat(table2array(data(1,OptSetup.lockeddof_indata)),1,N);
L(ilocked_dof) = repmat(table2array(data(1,OptSetup.lockeddof_indata)),1,N);
U(ilocked_ddof) = zeros(N*nlockeddofs,1);
L(ilocked_ddof) = zeros(N*nlockeddofs,1);

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

% Alternatively, we load an earlier solution

if strcmp(initialguess, 'mid')
    X0 = (L + U)/2;						% halfway between upper and lower bound
    X0 = X0 + 0.001;					% to avoid exact zero velocity, for which Autolev equations do not work
    X0(imeas) = data_dofs;
    X0(idmeas) = data_dofs_deriv;
elseif numel(strfind(initialguess, 'random')) > 0
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
    X0(imeas) = data_dofs;
    X0(idmeas) = data_dofs_deriv;
elseif numel(strfind(initialguess, 'init')) > 0
    eq_pos = load('equilibrium'); % equilibrium position
    X0 = zeros(nvar,1);
    ix = 1:nstates;
    for i_node=1:N
        X0(ix) = eq_pos.Result.x;
        ix = ix + nvarpernode;
    end
else
    % load a previous solution, initialguess contains file name
    initg = load(initialguess);
    t0 = initg.Result.times;
    x0 = initg.Result.x';
    u0 = initg.Result.u';
    
    % interpolate states and controls from initial guess to the current time grid
    if length(t0)>1
        x_0 = interp1(t0,x0,times,'linear','extrap');
        u_0 = interp1(t0,u0,times,'linear','extrap');
    else
        x_0 = repmat(x0,length(times),1);
        u_0 = repmat(u0,length(times),1);
    end
    X0 = reshape([x_0 u_0]',nvar,1);
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
        %figure; spy(cjac_num);
        cjac_full = full(cjac);
        [maxerr,irow] = max(abs(cjac-cjac_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d, %8.8f percent error\n', maxerr, irow(icol), icol, maxerr/cjac_full(irow(icol),icol));
        [maxerr,irow] = max(abs(grad_num-grad));
        fprintf('Max.error in objective gradient:  %8.5f at %d\n', maxerr, irow);
        keyboard
    end
    
    % evaluate initial guess
    print_flag=1;
    objfun(X0);
    confun(X0);
    print_flag=0;
    
    % set up constraints
    cl = zeros(ncon_eq,1);
    if Humcon_flag, cl = [cl; -Inf*ones(N-1,1)]; end
    if Scapcon_flag, cl = [cl; -Inf*ones(2*N,1)]; end
    cu = zeros(ncon_eq,1);
    if Humcon_flag, cu = [cu; zeros(N-1,1)]; end
    if Scapcon_flag, cu = [cu; zeros(2*N,1)]; end
    
    
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
    options.ipopt.print_info_string = 'yes';
    options.ipopt.limited_memory_max_history = 6;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
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
Result.u = x(nstates+(1:ncontrols),:);
angles = x(1:ndof,:);

% Calculate muscle lengths and forces
Result.mus_lengths = zeros(size(Result.u));
Result.mus_forces = zeros(size(Result.u));
Result.thor_hum = zeros(3,N);
Result.Fscap = zeros(2,N);
Result.FGHcontact = zeros(1,N);

for i_node=1:Result.OptSetup.N
    Result.mus_lengths(:,i_node) = das3('Musclelengths', Result.x(:,i_node));
    Result.mus_forces(:,i_node) = das3('Muscleforces', Result.x(:,i_node));
    [~, ~, ~, ~, ~, ~, Result.thor_hum(:,i_node)] = das3('Dynamics', Result.x(:,i_node), 0*Result.x(:,i_node),zeros(nmus,1));
    Result.Fscap(:,i_node) = das3('Scapulacontact', Result.x(:,i_node))';
end

ixs = 1:nstates;
ius = nstates+(1:ncontrols);
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

    Result.FGHcontact(i_node) = calculate_FGH(FGH);

    % advance the indices
    ixs = ixs + nvarpernode;
    ius = ius + nvarpernode;
end

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
        wf1 = mean(((X(imeas)-data_dofs)./data_dofs_sd).^2);
        
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
        f1 = Wdata * wf1;
        
        % Second term is mean squared muscle activation
        wf2 = mean(X(iact).^2);
        f2 = Weffort * wf2;
        
        % Third term is mean squared fibre velocity
        ixs = 2*ndof+1:2*ndof+nmus; % indices of fibre lengths of node 1
        wf3 = 0;
        for i=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            h = times(i+1) - times(i);
            wf3 = wf3 + mean(((x2-x1)/h).^2);
            % advance the indices
            ixs = ixs + nvarpernode;
        end
        f3 = 0.001 * wf3;
        
        % Fourth term is mean squared angular acceleration
        ixs = ndof+1:2*ndof; % indices of angular velocities of node 1
        wf4 = 0;
        for i=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            h = times(i+1) - times(i);
            wf4 = wf4 + mean(((x2-x1)/h).^2);
            % advance the indices
            ixs = ixs + nvarpernode;
        end
        f4 = 0.001 * wf4;
        
        f = f1 + f2 + f3 + f4;
        
        if print_flag
            obj_str = sprintf('Objfun (weighted): %9.5f = %9.5f (fit) + %9.5f (effort) + %9.5f (fibre vel) + %9.5f (angular acc)', f,f1,f2,f3,f4);
            wobj_str = sprintf('Objfun (unweighted): %9.5f (fit) + %9.5f (effort) + %9.5f (fibre vel) + %9.5f (angular acc)', wf1,wf2,wf3,wf4);
        end
        
        % Fifth term is thorax-scapula constraint
        if ~Scapcon_flag && Wscap
            Fscap = zeros(N,2);
            
            ixs = 1:nstates;
            for inode=1:N
                Fscap(inode,:) = das3('Scapulacontact', X(ixs));
                ixs = ixs + nvarpernode;
            end
            
            wf5 = mean(Fscap(:,1).^2)+mean(Fscap(:,2).^2);
            f5 = Wscap * wf5;
            f = f + f5;
            if print_flag
                obj_str5 = sprintf('+ %9.5f (scapula)', f5);
                wobj_str5 = sprintf('+ %9.5f (scapula)', wf5);
                obj_str = strcat(obj_str,obj_str5);
                wobj_str = strcat(wobj_str,wobj_str5);
            end
            Result.Fscap = Fscap;
        end
        
        % Sixth term is glenohumeral stability constraint
        if ~Humcon_flag && Whum
            
            FGHcontact = zeros(N-1,1);
            
            % indices of states and controls of node 1
            ixs = 1:nstates;
            ius = nstates+(1:ncontrols);
            % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
            for i=1:N-1
                x1 = X(ixs);
                x2 = X(ixs+nvarpernode);
                u1 = X(ius);
                u2 = X(ius+nvarpernode);
                h = times(i+1) - times(i);		% time interval
                
                handF1 = hand_force(i*3-2:i*3);
                handF2 = hand_force((i+1)*3-2:(i+1)*3);
                [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                
                FGHcontact(i) = calculate_FGH(FGH);
                
                % advance the indices
                ixs = ixs + nvarpernode;
                ius = ius + nvarpernode;
            end
            
            wf6 = mean(FGHcontact.^2);
            f6 = Whum * wf6;
            f = f + f6;
            if print_flag
                obj_str6 = sprintf('+ %9.5f (humerus)', f6);
                wobj_str6 = sprintf('+ %9.5f (humerus)', wf6);
                obj_str = strcat(obj_str,obj_str6);
                wobj_str = strcat(wobj_str,wobj_str6);
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
        g(imeas) = g(imeas) + 2*Wdata*(X(imeas)-data_dofs)./data_dofs_sd.^2/length(imeas);
        
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
        
        % Second term is mean squared muscle activation
        g(iact) = g(iact) + 2*Weffort*(X(iact))./(N*nmus);
        
        % Third term is mean squared fibre velocity
        ixs = 2*ndof+1:2*ndof+nmus;
        for i=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            h = times(i+1) - times(i);		% time interval
             g(ixs) = g(ixs) - 2*0.001*((x2-x1)/h^2)/nmus;
             g(ixs+nvarpernode) = g(ixs+nvarpernode) + 2*0.001*((x2-x1)/h^2)/nmus;
            % advance the indices
            ixs = ixs + nvarpernode;
        end
        
        % Fourth term is mean squared angular acceleration
        ixs = ndof+1:2*ndof;
        for i=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            g(ixs) = g(ixs) - 2*0.001*((x2-x1)/h^2)/ndof;
            g(ixs+nvarpernode) = g(ixs+nvarpernode) + 2*0.001*((x2-x1)/h^2)/ndof;
            % advance the indices
            ixs = ixs + nvarpernode;
        end
        
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
            iu1 = nstates+(1:ncontrols);
            for i=1:N-1
                ix2 = ix1 + nvarpernode;
                iu2 = iu1 + nvarpernode;
                x1 = X(ix1);
                x2 = X(ix2);
                u1 = X(iu1);
                u2 = X(iu2);
                h = times(i+1) - times(i);		% time interval
                handF1 = hand_force(i*3-2:i*3);
                handF2 = hand_force((i+1)*3-2:(i+1)*3);
                
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
        ceq = zeros(ncon_eq,1);
        
        % Calculate constraints due to discretized dynamics
        % indices of states and controls of node 1
        ixs = 1:nstates;
        ius = nstates+(1:ncontrols);
        % indices of equality constraints at node 1
        iceq = 1:(ncon_eq/(N-1));
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            x1 = X(ixs);
            x2 = X(ixs+nvarpernode);
            u1 = X(ius);
            u2 = X(ius+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            
            handF1 = hand_force(i*3-2:i*3);
            handF2 = hand_force((i+1)*3-2:(i+1)*3);
            [f, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);            

            f(OptSetup.lockeddof_inx) = zeros(nlockeddofs,1);
            f(ndof+OptSetup.lockeddof_inx) = zeros(nlockeddofs,1);
            ceq(iceq) = f;
            
            if Humcon_flag, c(i) = calculate_FGH(FGH)-1; end
            
            % advance the indices
            ixs = ixs + nvarpernode;
            ius = ius + nvarpernode;
            iceq = iceq + (ncon_eq/(N-1));
        end
        
        % Now calculate constraints due to scapular stability at each node
        if Scapcon_flag
            % indices of states and constraints of node 1
            ixs = 1:nstates;
            ic = Humcon_flag*(N-1) + (1:2);
            for inode=0:N-1
                % solve thorax ellipsoid surface equation for TS and AI, to find
                % out whether they are inside, on, or outside the thorax
                c(ic) = (das3('Scapulacontact', X(ixs)))';
                
                % advance the indices
                ixs = ixs + nvarpernode;
                ic = ic + 2;
            end
        end
        
        allcon = [ceq; c];
        
        if (print_flag)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq));
            fprintf('Norm(cnoneq): %9.5f  \n', norm(c));
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
        iu1 = nstates+(1:ncontrols);
        % indices of equality constraints at node 1
        iceq = 1:(ncon_eq/(N-1));
        % index of non-equality glenohumeral constraint at node 1
        ic = ncon_eq+1;
        dx = 1e-7;
        
        % evaluate dynamics at each pair of successive nodes, using trapezoidal integration formula
        for i=1:N-1
            ix2 = ix1 + nvarpernode;
            iu2 = iu1 + nvarpernode;
            x1 = X(ix1);
            x2 = X(ix2);
            u1 = X(iu1);
            u2 = X(iu2);
            h = times(i+1) - times(i);		% time interval
            handF1 = hand_force(i*3-2:i*3);
            handF2 = hand_force((i+1)*3-2:(i+1)*3);
            
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
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = ix1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdx/2 + dfdxdot/h);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = ix2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = iu1(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            [r,c,v] = find(dfdu/2);
            datal = length(v);
            allrows(index:index+datal-1) = iceq(1)+r-1;
            allcols(index:index+datal-1) = iu2(1)+c-1;
            allvals(index:index+datal-1) = v;
            index = index+datal;
            
            
            if Humcon_flag
                % Glenohumeral stability constraint
                FGHcontact = calculate_FGH(FGH)-1;
                
                % Perturb x1 and x2 by dx to estimate derivative
                for istate = [1:2*ndof 2*ndof+iGH 2*ndof+nmus+iGH]
                    xisave = x1(istate);
                    x1(istate) = x1(istate) + dx;
                    
                    % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH)-1;
                    
                    allrows(index) = ic;
                    allcols(index) = ix1(1)+istate-1;
                    allvals(index) = (FGHcontact_dx-FGHcontact)/dx;
                    index=index+1;
                    x1(istate) = xisave;
                    
                    xisave = x2(istate);
                    x2(istate) = x2(istate) + dx;
                    
                    % estimate delta(GHconstraint)/delta(x)
                    [~, ~, ~, ~, FGH] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,(handF1+handF2)/2);
                    FGHcontact_dx = calculate_FGH(FGH)-1;
                    
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
            iceq = iceq + (ncon_eq/(N-1));
        end
        
        if Scapcon_flag
            % indices of states of node 1
            ixs = 1:nstates;
            ic = ncon_eq + Humcon_flag*(N-1) + (1:2);
            dx = 1e-7;
            
            for inode=0:N-1
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

end % end of function das3_optimize
