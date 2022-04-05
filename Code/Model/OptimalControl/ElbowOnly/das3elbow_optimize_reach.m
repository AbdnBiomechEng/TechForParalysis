function Result = das3elbow_optimize_reach(filename,OptSetup)
% This program optimizes das3 movements to fit motion capture data
% It uses das3_thor.osim, which includes DoFs at the thorax, but locks those to the input data
% It also locks shoulder angles, and just optimises elbow flexion/extension
% and pronation/supination
%
%
% Inputs:
%		filename:	the output file name
%
%		OptSetup: settings for the optimization:
%
%					weights for the terms in the objective function: 
%						Wdata (kinematic tracking), Weffort (minimization of muscle effort)
%%
%					initialguess, with options: "init" (the equilibrium), "random"
%						"mid", "eqLce" (matches the kinematic input, with tendons at slack length)
% 						and if none of the above, assumed to be name of file with previous solution
%									
%					MaxIter: maximum number of iterations
%
%					"equality_constraints" flag to turn
% 						the constraints on and off (for testing). It is
% 						optional: if it is missing it is assumed to
% 						be 1 (constraints on)
%
%					OptimTol and FeasTol: Optimality and feasibility tolerances
%
% Output: 
%		Result:		The output of the optimization, including the settings in OptSetup
%						and the input data. The output is also saved in <filename>.mat 
%						(which can be used as a previous solution in initialguess)
%						and <filename>.sto, for visualization in OpenSim 				
%
% Dimitra Blana, March 2022, based on das3_optimize.m
	
% optimization settings
Result.OptSetup = OptSetup;

% weights for terms in objective function
% kinematic tracking
Wdata = OptSetup.Wdata;
% minimization of muscle effort
Weffort = OptSetup.Weffort;

% check whether constraint flag is present
if ~isfield(OptSetup, 'equality_constraints')
    OptSetup.equality_constraints = 1;
end

% lock thorax and shoulder dofs
lockeddofs = 1:12;
nlockeddofs = length(lockeddofs);

nmeas_dof = 14;   

N = OptSetup.N;
initialguess = OptSetup.initialguess;

% additional model inputs

% Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
ex_mom = zeros(5,1); 
% Vertical force of amplitude arm_support(2) applied to the ulna at a distance of arm_support(1) (meters) from the elbow
% to simulate arm support
arm_support = zeros(2,1);
% Force at the hand is read in 
hand_force = OptSetup.hand_force;

% Load structure with model related variables
modelparams = load(OptSetup.model_file); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;
ncontrols = nmus;

% some constants

nvarpernode = nstates + ncontrols;   % number of unknowns per time node
nvar = nvarpernode * N;              % total number of unknowns
ncon_eq = nstates*(N-1);             % number of constraints due to discretized dynamics
ncon_neq = 0;
ncon = ncon_eq + ncon_neq;

dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% Get maximum muscle activation (normally 1, but 0 for paralysed muscles,
% and 0.5 for FES muscles)
maxact = zeros(nmus,1);
for imus=1:mus
    maxact(imus) = model.muscles{imus}.maxact;
end


% Initialize the model
das3('Initialize',model);

% these are the range of motion limits
xlims = das3('Limits')';

LCEopt = das3('LCEopt');
LCEopt_N = repmat(LCEopt,N,1);

SEEslack = das3('SEEslack');

muscle_mass = das3('MuscleMass');
mass_N = repmat(muscle_mass,N,1);

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
        zeros(nmus,1) ];                                                % neural excitations

    U(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,2)+0.1;         % q
        (zeros(ndof,1) + 40);                                           % qdot
        zeros(nmus,1) + 1.7;                                            % Lce
        maxact;                                                         % active states
        maxact ];                                                       % neural excitations
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
t = [0;1];  % the time vector
data_init = zeros(2,8);

% Thoracohumeral angles from 5 to 60 degrees of flexion
data_init(1,4:6) = [90,5,0]*pi/180;
data_init(2,4:6) = [90,60,0]*pi/180;

% Elbow flexion-extension from 60 to 0 degrees 
data_init(1,7) = 60*pi/180;  
data_init(2,7) = 0*pi/180; 

% Pronation-supination set at 90 degrees
data_init(1:2,8) = 90*pi/180;  
    
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
print_flag = 0;

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
    iact = [iact nvarpernode*i_node+2*ndof+nmus+(1:nmus)];
    iexc = [iexc nvarpernode*i_node+2*ndof+2*nmus+(1:nmus)];
    
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

%% Make an initial guess

% Initial guess for kinematics + velocities is the interpolated measured data
% For activations and fibre lengths there are different options: mid,
% random and eqLce

% Alternatively, we load an earlier solution

if strcmp(initialguess, 'mid')
    X0 = (L + U)/2;						% halfway between upper and lower bound
    X0 = X0 + 0.001;					% to avoid exact zero velocity, for which Autolev equations do not work
    X0(idof) = datavec;
    X0(iddof) = data_deriv_vec;
elseif numel(strfind(initialguess, 'random')) > 0
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
    X0(idof) = datavec;
    X0(iddof) = data_deriv_vec;
elseif numel(strfind(initialguess, 'eqLce')) > 0  % muscle lenghts set so that tendons are slack
    X0 = zeros(nvar,1);
    X0(idof) = datavec;
    X0(iddof) = data_deriv_vec;
    ix = 1:nstates;
    equil_Lce = zeros(nmus*N,1);
    for i_node=1:N
        mustendon_lengths = das3('Musclelengths', X0(ix));
        equil_Lce(i_node*nmus-nmus+1:i_node*nmus) = (mustendon_lengths-SEEslack)./LCEopt;
        ix = ix + nvarpernode;
    end
    X0(iLce) = equil_Lce;
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
    if N>1
        X0 = reshape([x_0 u_0]',nvar,1);
    else
        X0 = reshape(x_0',nvar,1);
    end
end

%% Run optimization
if (OptSetup.MaxIter > 0)
    
    % determine sparsity structure of Jacobian (midpoint discretization)
    % Ton: I have verified that we always get same structure by testing with random X
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
        X = X0;
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
    cu = zeros(ncon_eq,1);
				
    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;

    if OptSetup.equality_constraints || OptSetup.inequality_constraints
        options.cl = cl;
        options.cu = cu;
        funcs.constraints = @confun;
        funcs.jacobian    = @conjac;
        funcs.jacobianstructure = @conjacstructure;
    end

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
    X = X0;
end

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

% Calculate muscle lengths, moment arms and forces
Result.elflex_mom_arms = zeros(size(Result.u));
Result.mus_lengths = zeros(size(Result.u));
Result.mus_forces = zeros(size(Result.u));

for i_node=1:Result.OptSetup.N
    momentarms = full(das3('Momentarms', Result.x(:,i_node)));
    Result.elflex_mom_arms(:,i_node) = momentarms(:,13);
    Result.mus_lengths(:,i_node) = das3('Musclelengths', Result.x(:,i_node));
    Result.mus_forces(:,i_node) = das3('Muscleforces', Result.x(:,i_node));
end

save([filename '.mat'],'Result');
make_osimm(filename,dofnames,angles,times);

%% Plot elbow angle and velocity, and muscle excitations and forces
figure;
ndof = Result.ndof;
nmus = Result.nmus;
subplot(6,1,1); plot(Result.times,Result.x(13:14,:)'*180/pi,'o-'); ylabel('Angle (degrees)');
title(Result.obj_str)
subplot(6,1,2); plot(Result.times,Result.x(ndof+(13:14),:)'*180/pi,'o-'); ylabel('Velocity (degrees/second)'); legend('Flexion/extension','Pronation/supination');
subplot(6,1,3); plot(Result.times,Result.u,'o-'); ylabel('Excitations');
subplot(6,1,4); plot(Result.times,Result.x(2*ndof+nmus+1:2*ndof+2*nmus,:),'o-'); ylabel('Activations');
subplot(6,1,5); plot(Result.times,Result.x(2*ndof+1:2*ndof+nmus,:),'o-'); ylabel('Norm fibre lengths');
subplot(6,1,6); plot(Result.times,Result.mus_forces','o-'); ylabel('Forces (N)');
xlabel('time (s)');


%% Functions to specify objective and constraints (and their derivatives)

% Objective function
    function f = objfun(X)
        
%         indices of states of node 1
%         ix = 1:nstates;
%         forces = zeros(nmus*N,1);
%         
%         for inode=1:N
%             forces(inode*nmus-nmus+1:inode*nmus) = das3('Muscleforces', X(ix));
%             ix = ix + nvarpernode;
%         end
        
        % First term of cost function is mean of squared differences to measured data
		wf1 = mean(((X(iElbow)-datavec(imeas_elbow))./datasd(imeas_elbow)).^2);
        f1 = Wdata * wf1;
        
        % Second term is mean squared muscle activation
        wf2 = mean(X(iact).^2);
%         
%         % Energy-consumption term 1
%         ect1 = mean(forces.*LCEopt_N);
%         
%         % Energy-consumption term 2
%         ect2 = mean(100*mass_N.*X(iact));
%         
%         % Energy-consumption term 3
%         ect3 = mean(400*mass_N.*X(iact).^2);
%         
%         wf2 = ect1 + ect2 + ect3;
        
        f2 = Weffort * wf2;

        % Third term is mean squared fibre velocity
        % indices of fibre lengths of node 1
        ix = 2*ndof+1:2*ndof+nmus;
        wf3 = 0;
        for i=1:N-1
            x1 = X(ix);
            x2 = X(ix+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            wf3 = wf3 + mean(((x2-x1)/h).^2);
            % advance the indices
            ix = ix + nvarpernode;
        end        
        f3 = 0.001 * wf3;

        % Fourth term is mean squared angular acceleration
        % indices of angular velocities of node 1
        ix = ndof+1:2*ndof;
        wf4 = 0;
        for i=1:N-1
            x1 = X(ix);
            x2 = X(ix+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            wf4 = wf4 + mean(((x2-x1)/h).^2);
            % advance the indices
            ix = ix + nvarpernode;
        end        
        f4 = 0.001 * wf4;
        
        
        f = f1 + f2 + f3 + f4;
        
		if print_flag
			obj_str = sprintf('Objfun (weighted): %9.5f = %9.5f (fit) + %9.5f (effort) + %9.5f (fibre vel) + %9.5f (angular acc)', f,f1,f2,f3,f4);
            wobj_str = sprintf('Objfun (unweighted): %9.5f (fit) + %9.5f (effort) + %9.5f (fibre vel) + %9.5f (angular acc)', wf1,wf2,wf3,wf4);
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
        
%        % indices of states of node 1
%         ix = 1:nstates;
%         ixus = 1:nvarpernode;
%         dx = 1e-7;
        nmeas = numel(iElbow);
        nact = numel(iact);
        
        % First term of cost function is mean of squared differences to
        % measured data
        g(iElbow) = g(iElbow) + 2*Wdata*(X(iElbow)-datavec(imeas_elbow))./datasd(imeas_elbow).^2/nmeas;
        
%         % Second term is energy consumption, in three parts:
%         for inode=1:N					
%             forces = das3('Muscleforces', X(ix));
% 			            			
%             for ivar=1:nvarpernode
%                 xisave = X(ixus(ivar));
%                 X(ixus(ivar)) = X(ixus(ivar)) + dx;
%                 
%                 forces_dx = das3('Muscleforces', X(ix));
%                 dforces = (forces_dx-forces)/dx;
%                 
%                 % Energy-consumption part 1: mean(forces.*LCEopt_N);
%                 g(ixus(ivar)) = g(ixus(ivar)) + Weffort * sum(LCEopt.*dforces/nact);
%                 				
%                 X(ixus(ivar)) = xisave;
%             end
%             ix = ix + nvarpernode;
%             ixus = ixus + nvarpernode;
%         end
        
        % Energy-consumption part 2: mean(100*mass_N*X(iact));
        %g(iact) = g(iact) + 100*Weffort*mass_N/nact;
        
        % Energy-consumption part 3: mean(400*mass_N*X(iact).^2);
        %g(iact) = g(iact) + 2*400*Weffort*mass_N.*X(iact)/nact;
        		
        
        % Second term is mean squared muscle activation
        g(iact) = g(iact) + 2*Weffort*(X(iact))./nact;
        
        % Third term is mean squared fibre velocity
        ix = 2*ndof+1:2*ndof+nmus;
        for i=1:N-1
            x1 = X(ix);
            x2 = X(ix+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            g(ix) = g(ix) - 2*0.001*((x2-x1)/h^2)/nmus;
            g(ix+nvarpernode) = g(ix+nvarpernode) + 2*0.001*((x2-x1)/h^2)/nmus;
            % advance the indices
            ix = ix + nvarpernode;
        end        

        % Fourth term is mean squared angular acceleration
        ix = ndof+1:2*ndof;
        for i=1:N-1
            x1 = X(ix);
            x2 = X(ix+nvarpernode);
            h = times(i+1) - times(i);		% time interval
            g(ix) = g(ix) - 2*0.001*((x2-x1)/h^2)/ndof;
            g(ix+nvarpernode) = g(ix+nvarpernode) + 2*0.001*((x2-x1)/h^2)/ndof;
            % advance the indices
            ix = ix + nvarpernode;
        end        
        
		        
    end

% Constraints
    function [c,ceq] = confun_fmincon(X)
        
        % Linear constraints: ceq(x) = 0
        c = [];
        ceq = zeros(ncon_eq,1);
        
        if ~OptSetup.equality_constraints
            return;
        end
                       
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

            f = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,hand_force);

            if OptSetup.equality_constraints
                f(lockeddofs) = zeros(length(lockeddofs),1);
                f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
                ceq(iceq) = f;
            end

            % advance the indices
            ix = ix + nvarpernode;
            ius = ius + nvarpernode;
            iceq = iceq + (ncon_eq/(N-1));
        end

        
        if (print_flag)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq));
        end
    end

% Linear and nonlinear constraints in one output, for IPOPT
    function allcon = confun(X)
        [c,ceq] = confun_fmincon(X);
        allcon = [ceq; c];
    end
        
% Structure of Jacobian matrix        
    function J = conjacstructure()
        J = Jpattern;
    end


% Jacobian of constraints
    function [nonJ,J] = conjac_fmincon(X)
        
        nonJ=[];
        %J = spalloc(ncon, nvar, Jnnz);
        
        if ~OptSetup.equality_constraints 
            J = [];
            return;
        end
        
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

            [~, dfdx, dfdxdot, dfdu] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2,ex_mom,arm_support,hand_force);

            if OptSetup.equality_constraints

                dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
                dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
                dfdxdot(lockeddofs,:) = zeros(length(lockeddofs),nstates);
                dfdxdot(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
                dfdu(lockeddofs,:) = zeros(length(lockeddofs),nmus);
                dfdu(ndof+lockeddofs,:) = zeros(length(lockeddofs),nmus);

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

            end

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

% Returns jacobians of linear and nonlinear constraints separately, for IPOPT
    function J = conjac(X)
        [~,J] = conjac_fmincon(X);
    end

end		% end of function das3elbow_optimize_reach
