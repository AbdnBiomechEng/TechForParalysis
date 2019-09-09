function Result = das3_optimize_elbow(data,t,filename,OptSetup)
% This program optimizes das3 movements to fit motion capture data
% It uses das3_thor.osim, which includes DoFs at the thorax, but locks those to the input data
%
%
% Inputs:
%		data: 	The kinematic data with one angle per column (14 in total)
%
%		t:		the time vector (same length as the number of rows in data)
%
%		filename:	the output file name
%
%		OptSetup: settings for the optimization:
%
%					weights for the terms in the objective function: 
%						Wdata (kinematic tracking), Weffort (minimization of muscle effort)
%
%					N: the number of nodes for Direct Collocation. 
%						If N==1, it finds a static equilibrium.
%
%					initialguess, with options: "init" (the equilibrium), "random"
%						"mid", "eqLce" (matches the kinematic input, with tendons at slack length)
% 						and if none of the above, assumed to be name of file with previous solution
%									
%					MaxIter: maximum number of iterations
%
%					solver: IPOPT or SNOPT			
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
% Dimitra Blana, September 2019, based on das3_optimize.m
	

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

% Load structure with model related variables
% modelparams = load('model_struct.mat'); this has all the muscles
modelparams = load('model_struct_elbow.mat'); % this has only brachialis
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;
ncontrols = nmus;

% some constants
if N>1
    nvarpernode = nstates + ncontrols;   % number of unknowns per time node
    nvar = nvarpernode * N;              % total number of unknowns
    ncon_eq = nstates*(N-1);             % number of constraints due to discretized dynamics
    ncon_neq = 0;
else
    nvarpernode = nstates;              % number of unknowns per time node
    nvar = nvarpernode * N;             % total number of unknowns
    ncon_eq = 2*ndof + nmus;            % number of constraints due to discretized dynamics
    ncon_neq = 0;
end
ncon = ncon_eq + ncon_neq;

dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
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
        (zeros(ndof,1) - 40)*pi/180;                                    % qdot
        zeros(nmus,1) + 0.3;                                            % Lce
        zeros(nmus,1);                                                  % active states
        zeros(nmus,1) ];                                                % neural excitations

    U(i_node*nvarpernode + (1:nvarpernode) ) = [xlims(:,2)+0.1;         % q
        (zeros(ndof,1) + 40)*pi/180;                                    % qdot
        zeros(nmus,1) + 1.7;                                            % Lce
        ones(nmus,1);                                                   % active states
        ones(nmus,1) ];                                                 % neural excitations
end

% Bounds for angular velocities are different in first and final node: they
% should be zero

% First node
L(ndof+1:2*ndof) = zeros(ndof,1);
U(ndof+1:2*ndof) = zeros(ndof,1);

% Final node
final_node_start_index = (N-1)*nvarpernode;
L(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);
U(final_node_start_index + (ndof+1:2*ndof)) = zeros(ndof,1);

% load the motion capture data, columns are time and 14 angles
Result.input = data;
Result.input_t = t;

% make the vector of times for direct collocation
t1 = min(t);
t2 = max(t);
if N>1
    times = t1 + (0:N-1)'*(t2-t1)/(N-1);
else
    times = t1 + (t2-t1)/2;
end

% resample the motion capture data onto the direct collocation times
% and put into a long vector to save time when computing cost function
if N>1
    data = interp1(t,data,times,'pchip','extrap');
else
    if size(data,1)>1, data = mean(data); end
end

datavec = reshape(data', N*(nmeas_dof), 1);

% define the estimated uncertainty in each measured variable (in radians)
datasd = ones(nmeas_dof,1);         % for one node
datasd = repmat(datasd,N,1);		% replicate for all nodes
print = 0;

% precompute the indices for activations and some angles, so we can compute cost function quickly
iact = [];
iLce = [];

iThorSh = [];
imeas_ThorSh = [];
iElbow = [];
imeas_elbow = [];

for i_node=0:N-1
    iact = [iact nvarpernode*i_node+2*ndof+nmus+(1:nmus)];
    iLce = [iLce nvarpernode*i_node+2*ndof+(1:nmus)];
    
    iThorSh = [iThorSh nvarpernode*i_node+lockeddofs];
    iElbow = [iElbow nvarpernode*i_node+nlockeddofs+(1:2)];
    imeas_ThorSh = [imeas_ThorSh nmeas_dof*i_node+lockeddofs];
	imeas_elbow = [imeas_elbow nmeas_dof*i_node+nlockeddofs+(1:2)];
end

% the thorax and shoulder are set at input data
U(iThorSh) = datavec(imeas_ThorSh);
L(iThorSh) = datavec(imeas_ThorSh);

% make an initial guess
if strcmp(initialguess, 'mid')
    X0 = (L + U)/2;						% halfway between upper and lower bound
    X0 = X0 + 0.001;					% to avoid exact zero velocity, for which Autolev equations do not work
elseif numel(strfind(initialguess, 'random')) > 0
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
elseif numel(strfind(initialguess, 'init')) > 0
    % xeq = load('equilibrium.mat'); this has all the muscles
    xeq = load('equilibrium_brac.mat'); % this has only brachialis
    if N>1
        X0 = reshape(repmat([xeq.x; xeq.x(2*ndof+nmus+1:end)],1,N),nvar,1);
    else
        X0 = xeq.x;
    end
elseif numel(strfind(initialguess, 'eqLce')) > 0
    % xeq = load('equilibrium.mat'); this has all the muscles
    xeq = load('equilibrium_brac.mat'); % this has only brachialis
    if N>1
        X0 = reshape(repmat([xeq.x; xeq.x(2*ndof+nmus+1:end)],1,N),nvar,1);
    else
        X0 = xeq.x;
    end
    X0(iElbow) = datavec(imeas_elbow);
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

X0(iThorSh) = datavec(imeas_ThorSh);

% run optimization
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
    print=1;
    objfun(X0);
    confun(X0);
    print=0;
    
	% set up constraints
   if N>1
		cl = zeros(ncon_eq,1);
		cu = zeros(ncon_eq,1);
		
	else
		cl = zeros(ncon_eq,1);
		cu = zeros(ncon_eq,1);
   end
		
    if strcmp(OptSetup.solver,'IPOPT')
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
        options.ipopt.print_level = 5;
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
        Result.status = info.status;
        Result.info = info;
        fprintf('IPOPT info status: %d\n', info.status);
    elseif strcmp(OptSetup.solver,'SNOPT')
        if OptSetup.equality_constraints || OptSetup.inequality_constraints
			Prob = conAssign(@objfun, @objgrad, [], [], L, U, 'das3', X0, ...
				[], 0, ...
				[], [], [], @confun, @conjac, [], Jpattern, ...
				cl, cu, ...
				[], [], [],[]);
        else
            Prob = conAssign(@objfun, @objgrad, [], [], L, U, 'das3', X0, ...
                [], 0, ...
                [], [], [], [], [], [], [], ...
                [], [], ...
                [], [], [],[]);
        end
        % Prob.SOL.optPar(1)= 1;
        Prob.PriLevOpt = 1;
        Prob.SOL.optPar(9) = OptSetup.FeasTol;
        Prob.SOL.optPar(10) = OptSetup.OptimTol;
        Prob.SOL.optPar(11) = 1e-6; % Minor feasibility tolerance (1e-6)
        Prob.SOL.optPar(30) = 1000000; % maximal sum of minor iterations (max(10000,20*m))
        Prob.SOL.optPar(35) = OptSetup.MaxIter;
        Prob.SOL.optPar(36) = 40000; % maximal number of minor iterations in the solution of the QP problem (500)
        Prob.SOL.moremem = 10000000; % increase internal memory
        %print = 1;		% we want something on the screen
        ResultS = tomRun('snopt',Prob);
        X = ResultS.x_k;
        fprintf('--------------------------------------\n');
        fprintf('SNOPT finished:\n')
        disp(ResultS.ExitText);
        fprintf('Number of iterations: %d\n',ResultS.Iter);
        Result.ExitText = ResultS.ExitText;
        Result.ExitFlag = ResultS.ExitFlag;
        Result.Inform = ResultS.Inform;
        Result.Iter = ResultS.Iter;
    end
    
else		% skip optimization
    X = X0;
end

beep;

print=1;
objfun(X);
confun(X);

% save this result on a file
x = reshape(X,nvarpernode,N);
Result.times = times;
Result.x = x(1:nstates,:);
if N>1
    Result.u = x(nstates+(1:ncontrols),:);
else
    Result.u = x(iact);
end
angles = x(1:ndof,:);

save([filename '.mat'],'Result');
make_osimm(filename,dofnames,angles,times);
beep;

% objective function
    function f = objfun(X)
        
        % indices of states of node 1
        ix = 1:nstates;
        forces = zeros(nmus*N,1);
		if N>1, ius = nstates+(1:ncontrols); else ius = iact; end
        
        for inode=1:N
            if N==1
                xdot = zeros(nstates,1);
            elseif inode==1
                xdot = (X(ix + nvarpernode) - X(ix))/(times(inode+1) - times(inode));
            elseif inode==N
                xdot = (X(ix) - X(ix - nvarpernode))/(times(inode) - times(inode-1));
            else
                xdot = (X(ix + nvarpernode) - X(ix - nvarpernode))/(times(inode+1) - times(inode-1));
            end			
            forces(inode*nmus-nmus+1:inode*nmus) = das3('Muscleforces', X(ix));

            ix = ix + nvarpernode;
            ius = ius + nvarpernode;
        end
        
        % first term of cost function is mean of squared differences to measured data
		wf1 = mean(((X(iElbow)-datavec(imeas_elbow))./datasd(imeas_elbow)).^2);
        f1 = Wdata * wf1;
        
        % second term is energy consumption, with Weffort weight
        
        % Energy-consumption term 1
        ect1 = mean(forces.*LCEopt_N);
        
        % Energy-consumption term 2
        ect2 = mean(100*mass_N.*X(iact));
        
        % Energy-consumption term 3
        ect3 = mean(400*mass_N.*X(iact).^2);
        
        wf2 = ect1 + ect2 + ect3;
        f2 = Weffort * wf2;
        
        f = f1 + f2;
        
		if print
			obj_str = sprintf('Objfun: %9.5f = %9.5f (fit) + %9.5f (effort)', f,f1,f2);
            wobj_str = sprintf('%9.5f (fit) + %9.5f (effort)', wf1,wf2);
		end
						
        if (print)
            fprintf(obj_str);
            fprintf('\n');
            fprintf(wobj_str);
            fprintf('\n');
        end
        
    end


% gradient of objective function
    function g = objgrad(X)
        g = zeros(nvar,1);
        
        % indices of states of node 1
        ix = 1:nstates;
		if N>1, ius = nstates+(1:ncontrols); else ius = iact; end
        ixus = 1:nvarpernode;
        dx = 1e-7;
        nmeas = numel(iElbow);
        nact = numel(iact);
        
        % first term of cost function is mean of squared differences to
        % measured data
        for inode=1:N
		
            if N==1
                xdot = zeros(nstates,1);
            elseif inode==1
                xdot = (X(ix + nvarpernode) - X(ix))/(times(inode+1) - times(inode));
            elseif inode==N
                xdot = (X(ix) - X(ix - nvarpernode))/(times(inode) - times(inode-1));
            else
                xdot = (X(ix + nvarpernode) - X(ix - nvarpernode))/(times(inode+1) - times(inode-1));
            end
			
            forces = das3('Muscleforces', X(ix));
			            			
            for ivar=1:nvarpernode
                xisave = X(ixus(ivar));
                X(ixus(ivar)) = X(ixus(ivar)) + dx;
                
                forces_dx = das3('Muscleforces', X(ix));
                dforces = (forces_dx-forces)/dx;
                % Energy-consumption term 1: mean(forces.*LCEopt_N);
                g(ixus(ivar)) = g(ixus(ivar)) + Weffort * sum(LCEopt.*dforces/nact);
                				
                X(ixus(ivar)) = xisave;
            end
            ix = ix + nvarpernode;
            ixus = ixus + nvarpernode;
        end
        		
        g(iElbow) = g(iElbow) + 2*Wdata*(X(iElbow)-datavec(imeas_elbow))./datasd(imeas_elbow).^2/nmeas;
        
        % Energy-consumption term 2: mean(100*mass_N*X(iact));
        g(iact) = g(iact) + 100*Weffort*mass_N/nact;
        
        % Energy-consumption term 3: mean(400*mass_N*X(iact).^2);
        g(iact) = g(iact) + 2*400*Weffort*mass_N.*X(iact)/nact;
		        
    end

% constraints
    function allcon = confun(X)
        
        % Linear constraints: ceq(x) = 0
        ceq = zeros(ncon_eq,1);
        
        if ~OptSetup.equality_constraints
            allcon = ceq;
            return;
        end
        
        if N>1
                
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

        else
            % Calculate constraints due to discretized dynamics

            % evaluate dynamics
            f = das3('Dynamics',X,zeros(nstates,1),X(iact));

            if OptSetup.equality_constraints
                f(lockeddofs) = zeros(length(lockeddofs),1);
                f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
                ceq = f(1:2*ndof+nmus);
            end

        end
        allcon = ceq;
        
        if (print)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq));
        end
    end


    function J = conjacstructure()
        % returns structure of constraint Jacobian matrix
        J = Jpattern;
    end


% Jacobian of constraints
    function J = conjac(X)
        
        J = spalloc(ncon, nvar, Jnnz);
        
        if ~OptSetup.equality_constraints 
            return;
        end
        
        if N>1
            %         allrows = zeros(Jnnz+con_lenJ,1);
            %         allcols = zeros(Jnnz+con_lenJ,1);
            %         allvals = zeros(Jnnz+con_lenJ,1);
            %         index=1;

            % indices of states and controls of node 1
            ix1 = 1:nstates;
            iu1 = nstates+(1:ncontrols);
            % indices of equality constraints at node 1
            iceq = 1:(ncon_eq/(N-1));
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

                [~, dfdx, dfdxdot, dfdu] = das3('Dynamics',(x1+x2)/2,(x2-x1)/h,(u1+u2)/2);

                if OptSetup.equality_constraints

                    dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
                    dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
                    dfdxdot(lockeddofs,:) = zeros(length(lockeddofs),nstates);
                    dfdxdot(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
                    dfdu(lockeddofs,:) = zeros(length(lockeddofs),nmus);
                    dfdu(ndof+lockeddofs,:) = zeros(length(lockeddofs),nmus);

                    % which generates four blocks in the Jacobian:
                    J(iceq,ix1) = dfdx/2 - dfdxdot/h;
                    %                 [r,c,v] = find(dfdx/2 - dfdxdot/h);
                    %                 datal = length(v);
                    %                 allrows(index:index+datal-1) = iceq(1)+r-1;
                    %                 allcols(index:index+datal-1) = ix1(1)+c-1;
                    %                 allvals(index:index+datal-1) = v;
                    %                 index = index+datal;

                    J(iceq,ix2) = dfdx/2 + dfdxdot/h;
                    %                 [r,c,v] = find(dfdx/2 + dfdxdot/h);
                    %                 datal = length(v);
                    %                 allrows(index:index+datal-1) = iceq(1)+r-1;
                    %                 allcols(index:index+datal-1) = ix2(1)+c-1;
                    %                 allvals(index:index+datal-1) = v;
                    %                 index = index+datal;

                    J(iceq,iu1) = dfdu/2;
                    %                 [r,c,v] = find(dfdu/2);
                    %                 datal = length(v);
                    %                 allrows(index:index+datal-1) = iceq(1)+r-1;
                    %                 allcols(index:index+datal-1) = iu1(1)+c-1;
                    %                 allvals(index:index+datal-1) = v;
                    %                 index = index+datal;

                    J(iceq,iu2) = dfdu/2;
                    %                 [r,c,v] = find(dfdu/2);
                    %                 datal = length(v);
                    %                 allrows(index:index+datal-1) = iceq(1)+r-1;
                    %                 allcols(index:index+datal-1) = iu2(1)+c-1;
                    %                 allvals(index:index+datal-1) = v;
                    %                 index = index+datal;

                end

                % advance the indices
                ix1 = ix1 + nvarpernode;
                iu1 = iu1 + nvarpernode;
                iceq = iceq + (ncon_eq/(N-1));
            end
        
        else
            % evaluate dynamics
            [~, dfdx, ~ ,~] = das3('Dynamics',X,zeros(nstates,1),X(iact));

            if OptSetup.equality_constraints

                dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
                dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);

                dfdxstar = dfdx(1:2*ndof+nmus,:);
                % Gradient matrix: transpose of the form of Jacobians
                % (one column per constraint)
                J(1:2*ndof+nmus,:) = dfdxstar;

            end
        end
        %        J = sparse(allrows(1:index-1),allcols(1:index-1),allvals(1:index-1),ncon, nvar);
    end

end		% end of function das3_optimize_elbow
