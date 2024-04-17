function Result = das3_optimize_DS(data,filename,OptSetup)
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
%					N: the number of nodes for Direct Shooting.
%					initialguess, which is either "initial_state" for the
%					first optimisation, or the name of file with previous solution
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
% scapular stability
Wscap = OptSetup.Wscap;
% humeral stability
Whum = OptSetup.Whum;

N = OptSetup.N;

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

% Indeces of unlocked dofs
ix = setdiff(1:nstates, [OptSetup.lockeddof_inx  ndof+OptSetup.lockeddof_inx]);
nstates_ix = length(ix);  % number of states in the reduced-dof system

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
nvarpernode = nstates + ncontrols;                % number of unknowns per time node
nvar = ncontrols * N;                             % total number of unknowns

% Initialize the model
das3('Initialize',model);

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

% Lower and upper bounds for the optimization variables X (the muscle
% excitations)
L = zeros(1,nvar);
U = repmat(max_act',1,N);

% Load the input data
Result.input = data(:,OptSetup.input_variables);
Result.input_t = data.time;

% Make the vector of times for direct shooting
t_start = min(data.time);
t_end = max(data.time);
tstep_nodes = (t_end-t_start)/(N-1);
times = t_start + (0:N-1)'*tstep_nodes;

% Resample the input data onto the direct shooting times
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
    data_thorhum_sd = 1000*[1 1 10]';    % for one node (set to higher number to only track start and end of movement)
    data_thorhum_sd = repmat(data_thorhum_sd,N,1);		% replicate for all nodes
    data_thorhum_sd(1:length(OptSetup.thorhum_indata))=1;
    data_thorhum_sd(length(OptSetup.thorhum_indata)*(N-1)+1:end)=1;
    
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

imeas = zeros(1,N*ntrackdofs);
idmeas = zeros(1,N*ntrackdofs);

ishoulder = zeros(1,N*9);
idshoulder = zeros(1,N*9);

for i_node=0:N-1
    iact(i_node*nmus+(1:nmus)) = nvarpernode*i_node+2*ndof+nmus+(1:nmus);
    iu(i_node*nmus+(1:nmus)) = nvarpernode*i_node+ nstates + (1:nmus);
    
    imeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+OptSetup.tracking_inx;
    idmeas(i_node*ntrackdofs+(1:ntrackdofs)) = nvarpernode*i_node+ndof+OptSetup.tracking_inx;
    
    ishoulder(i_node*9+(1:9)) = nvarpernode*i_node+3+(1:9);
    idshoulder(i_node*9+(1:9)) = nvarpernode*i_node+ndof+3+(1:9);
    
end

% Locked dofs stay at measured values, and have zero velocity
ilocked_dof = OptSetup.lockeddof_inx;
ilocked_ddof = ndof+OptSetup.lockeddof_inx;
nlockeddofs = length(ilocked_dof);

datalockeddofs_m = data(:,OptSetup.lockeddof_indata);
datalockeddofs = interp1(data.time,table2array(datalockeddofs_m),times);

datalockeddofs_deriv = diff(datalockeddofs)./diff(times);
datalockeddofs_deriv = [datalockeddofs_deriv; datalockeddofs_deriv(end,:)];

% convert into a single column
data_lockeddofs = reshape(datalockeddofs', [], 1);
data_lockeddofs_deriv = reshape(datalockeddofs_deriv', [], 1);

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

%% Load initial position
init_state = load('initial_state.mat');
x_init = init_state.x;
x0 = repmat(x_init',length(times),1);
u0 = init_state.x(2*ndof+nmus+(1:nmus))';    
U0 = repmat(u0,length(times),1);
U0long = reshape(U0',ncontrols*N,1);
X = reshape([x0 U0]',nvarpernode*N,1);

%% Set up optimization
tstep = tstep_nodes/(ceil(tstep_nodes/0.001));
opt_times = t_start:tstep:t_end;
nsteps = length(opt_times);
mult_N_steps = tstep_nodes / tstep;

%% Run optimization
print_flag=0;
cmaes;
Result.X = X;
Result.U0 = U0;
    
print_flag=1;
objfun(Result.X);
print_flag=0;

%% Save this result in a mat file and an Opensim motion file

x = reshape(Result.X,nstates+ncontrols,N);
Result.ndof = ndof;
Result.nmus = nmus;
Result.muscles = model.muscles;
Result.times = times;
Result.input_data = data;
Result.x = x(1:nstates,:);
Result.u = x(nstates+(1:nmus),:);
angles = x(1:ndof,:);

% Calculate muscle lengths and forces
Result.mus_lengths = zeros(size(Result.u));
Result.mus_forces = zeros(size(Result.u));

for i_node=1:Result.OptSetup.N
    Result.mus_lengths(:,i_node) = das3('Musclelengths', Result.x(:,i_node));
    Result.mus_forces(:,i_node) = das3('Muscleforces', Result.x(:,i_node));
end

save([filename '.mat'],'Result');
make_osimm(filename,dofnames,angles,times);


% Objective function
    function f = objfun(U)

        % Run forward simulation with current guess, and save all states
        [X,thorhum,handpos,Fscap,FGHcontact] = das3_forward_DS(U');
        
        % First term of cost function is mean of squared differences to
        % measured data, including thoraco-humeral data if provided
        if ntrackdofs
            wf1 = mean(((X(imeas)-data_dofs)./data_dofs_sd).^2);
        else
            wf1 = 0;
        end
        
        % If required, calculate thoracohumeral angles
        if thorhum_flag
            wf1 = wf1 + mean(((thorhum-data_thorhum)./data_thorhum_sd).^2);
        end
        
        % If required, calculate hand position
        if handpos_flag
            wf1 = wf1 + mean(((handpos-data_handpos)./data_handpos_sd).^2);
        end
        f1 = Wdata * wf1;

        % Second term is mean squared muscle activation        
        wf2 = mean(X(iact).^2);
        f2 = Weffort * wf2;
                                
        % Third term is thorax-scapula constraint
        wf3 = mean(Fscap(:,1).^2)+mean(Fscap(:,2).^2);
        f3 = Wscap * wf3;
                        
        % Fourth term is glenohumeral stability constraint            
        wf4 = mean(FGHcontact.^2);
        f4 = Whum * wf4;

        f = f1 +f2 + f3 + f4;

        % If required, minimize dof velocity at final node
        if OptSetup.end_at_rest
            % Final node
            final_node_start_index = (N-1)*nmus;
            dof_vel = X(final_node_start_index + (ndof+1:2*ndof));

            f5 = mean(dof_vel.^2);

            f = f + f5;
        end

        if print_flag
            obj_str = sprintf('Objfun (weighted): %9.5f = %9.5f (fit) + %9.5f (effort) + %9.5f (scapula) + %9.5f (GH)', f, f1, f2, f3, f4);
            wobj_str = sprintf('Objfun (unweighted): %9.5f (fit) + %9.5f (effort) + %9.5f (scapula) + %9.5f (GH)', wf1, wf2, wf3, wf4);
        end
                
        if print_flag
            fprintf(obj_str);
            fprintf('\n');
            fprintf(wobj_str);
            fprintf('\n');
        end
        
    end


    function cmaes
        %% Problem Settings

        CostFunction=@objfun;   % Cost Function
        %nVar=nvar;              % Number of Unknown (Decision) Variables
        nVar=6*N;              % Number of Unknown (Decision) Variables
        VarSize=[1 nVar];       % Decision Variables Matrix Size
        % VarMin= L;              % Lower Bound of Decision Variables
        % VarMax= U;              % Upper Bound of Decision Variables
        VarMin= 0;              % Lower Bound of Decision Variables
        VarMax= 1;              % Upper Bound of Decision Variables

        %% CMA-ES Settings

        % Maximum Number of Iterations
        MaxIt=OptSetup.MaxIter;

        % Population Size (and Number of Offsprings)
        lambda=(4+round(3*log(nVar)))*10;

        % Number of Parents
        mu=round(lambda/2);

        % Parent Weights
        w=log(mu+0.5)-log(1:mu);
        w=w/sum(w);

        % Number of Effective Solutions
        mu_eff=1/sum(w.^2);

        % Step Size Control Parameters (c_sigma and d_sigma);
        sigma0=0.3*(VarMax-VarMin);
        cs=(mu_eff+2)/(nVar+mu_eff+5);
        ds=1+cs+2*max(sqrt((mu_eff-1)/(nVar+1))-1,0);
        ENN=sqrt(nVar)*(1-1/(4*nVar)+1/(21*nVar^2));

        % Covariance Update Parameters
        cc=(4+mu_eff/nVar)/(4+nVar+2*mu_eff/nVar);
        c1=2/((nVar+1.3)^2+mu_eff);
        alpha_mu=2;
        cmu=min(1-c1,alpha_mu*(mu_eff-2+1/mu_eff)/((nVar+2)^2+alpha_mu*mu_eff/2));
        hth=(1.4+2/(nVar+1))*ENN;

        %% Initialization

        ps=cell(MaxIt,1);
        pc=cell(MaxIt,1);
        C=cell(MaxIt,1);
        sigma=cell(MaxIt,1);

        ps{1}=zeros(VarSize);
        pc{1}=zeros(VarSize);
        C{1}=eye(nVar);
        sigma{1}=sigma0;

        empty_individual.Position=[];
        empty_individual.Step=[];
        empty_individual.Cost=[];

        M=repmat(empty_individual,MaxIt,1);
        M(1).Position=unifrnd(VarMin,VarMax,VarSize);
        M(1).Step=zeros(VarSize);
        M(1).Cost=CostFunction(M(1).Position);

        BestSol=M(1);

        BestCost=zeros(MaxIt,1);

        %% CMA-ES Main Loop

        for g=1:MaxIt

            % Generate Samples
            pop=repmat(empty_individual,lambda,1);
            for i=1:lambda
                pop(i).Step=mvnrnd(zeros(VarSize),C{g});
                pop(i).Position=M(g).Position+sigma{g}*pop(i).Step;
                pop(i).Position(pop(i).Position>1) = 1;
                pop(i).Position(pop(i).Position<0) = 0;
                pop(i).Cost=CostFunction(pop(i).Position);

                % Update Best Solution Ever Found
                if pop(i).Cost<BestSol.Cost
                    BestSol=pop(i);
                end
                disp(['Sample generation ' num2str(i) ' out of ' num2str(lambda) ' finished ']);
            end

            % Sort Population
            Costs=[pop.Cost];
            [Costs, SortOrder]=sort(Costs);
            pop=pop(SortOrder);

            % Save Results
            BestCost(g)=BestSol.Cost;

            % Display Results
            disp(['Iteration ' num2str(g) ': Best Cost = ' num2str(BestCost(g))]);

            % Exit At Last Iteration
            if g==MaxIt
                break;
            end

            % Update Mean
            M(g+1).Step=0;
            for j=1:mu
                M(g+1).Step=M(g+1).Step+w(j)*pop(j).Step;
            end
            M(g+1).Position=M(g).Position+sigma{g}*M(g+1).Step;
            M(g+1).Cost=CostFunction(M(g+1).Position);
            if M(g+1).Cost<BestSol.Cost
                BestSol=M(g+1);
            end

            % Update Step Size
            ps{g+1}=(1-cs)*ps{g}+sqrt(cs*(2-cs)*mu_eff)*M(g+1).Step/chol(C{g})';
            sigma{g+1}=sigma{g}*exp(cs/ds*(norm(ps{g+1})/ENN-1))^0.3;

            % Update Covariance Matrix
            if norm(ps{g+1})/sqrt(1-(1-cs)^(2*(g+1)))<hth
                hs=1;
            else
                hs=0;
            end
            delta=(1-hs)*cc*(2-cc);
            pc{g+1}=(1-cc)*pc{g}+hs*sqrt(cc*(2-cc)*mu_eff)*M(g+1).Step;
            C{g+1}=(1-c1-cmu)*C{g}+c1*(pc{g+1}'*pc{g+1}+delta*C{g});
            for j=1:mu
                C{g+1}=C{g+1}+cmu*w(j)*pop(j).Step'*pop(j).Step;
            end

            % If Covariance Matrix is not Positive Defenite or Near Singular
            [V, E]=eig(C{g+1});
            if any(diag(E)<0)
                E=max(E,0);
                C{g+1}=V*E/V;
            end

        end
    end

    function [X,thorhum,handpos,Fscap,FGHcontact] = das3_forward_DS(Unew)
    % Runs a forward dynamic simulation from the initial state using the muscle
    % excitations in U
        
    % reserve space to store results
    xout = zeros(N, nstates);
    uall = zeros(N,nmus);
    thorhum = zeros(nthorhum*N,1);
    handpos = zeros(3*N,1);
    Fscap = zeros(N,2);
    FGHcontact = zeros(N,1);

    % run simulation

    x = x_init(ix);    
    % Initialize variables
    step_xdot = zeros(nstates_ix,1);
    step_u = u0';
    step_u(83:88) = Unew(1:6);
    %step_u = Unew(1:ncontrols);
    y = zeros(nstates,1);               % state vector for full system
    ydot = zeros(nstates,1);            % state vector derivative for full system   
    xfull = zeros(nstates,1);           % for storage
    inode = 1;

    for istep=1:nsteps
    
        % Evaluate dynamics in current x and xdot (states of the reduced-dof system)
        y(ilocked_dof) = data_lockeddofs((inode-1)*nlockeddofs+1:inode*nlockeddofs);
        y(ilocked_ddof) = data_lockeddofs_deriv((inode-1)*nlockeddofs+1:inode*nlockeddofs);
        xfull(ilocked_dof) = data_lockeddofs((inode-1)*nlockeddofs+1:inode*nlockeddofs); 

        y(ix) = x;              % put the unlocked state variables in their proper place inside y
        ydot(ix) = step_xdot;   % same for ydot
        % compute dynamics residuals for the full model, and ignore those that correspond to locked joints
        u1 = u0';
        u1(83:88) = Unew((inode-1)*6+1:inode*6);
        %u1 = Unew((inode-1)*ncontrols+1:inode*ncontrols);
        [g, dgdy, dgdydot, dgdu, F_GH, ~, thorhum_x] = das3('Dynamics',y, ydot, u1, ex_mom,arm_support,hand_force(inode*3-2:inode*3));  

        f = g(ix);
        dfdx = dgdy(ix,ix);
        dfdxdot = dgdydot(ix,ix);
        dfdu = dgdu(ix,:);
    
        % Solve the change in x from the 1st order Rosenbrock formula
        du = u1 - step_u;
        dx = (dfdx + dfdxdot/tstep)\(dfdxdot*step_xdot - f - dfdu*du);
    
        x = x + dx;
        
        % Update variables for the next simulation step
        step_xdot = dx/tstep;
        step_u = u1;

        if ~mod(istep-1,mult_N_steps)
            thorhum((inode-1)*nthorhum+(1:nthorhum)) = thorhum_x(OptSetup.thorhum_inx);
            handpos((inode-1)*3+(1:3)) = das3_handpos(y,local_hand_coord);
            Fscap(inode,:) = das3('Scapulacontact', y);
            FGHcontact(inode) = calculate_FGH(F_GH);      
            xfull(ix) = x;            % put the unlocked state variables in their proper place inside x
            xout(inode,:) = xfull';   % store result
            uall(inode,:) = u1';
            inode = inode+1;
        end
        
    end

    X = reshape([xout uall]',nvarpernode*N,1);
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
