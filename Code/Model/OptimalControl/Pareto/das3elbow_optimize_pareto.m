function Result = das3elbow_optimize_pareto(out_filename)
% This program finds pareto front of two cost functions: kinematic tracking and muscle effort
% 
% It uses das3_thor.osim, which includes DoFs at the thorax, but locks those to the input data
% It also locks shoulder angles, and just optimises elbow flexion/extension
% and pronation/supination	

% Which model to use
OptSetup.model_file = 'model_struct_elbow.mat'; 

% Input data
OptSetup.t = [0;1];  % the time vector
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

% some constants

nvarpernode = nstates + ncontrols;          % number of unknowns per time node
nvar = nvarpernode * N;                     % total number of unknowns
ncon_eq = nstates*(N-1);                    % number of constraints due to discretized dynamics

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

% Read partial initial population from previous optimisation runs
init_pop = load('init_pop_elbow_pareto');

options = optimoptions('gamultiobj','ConstraintTolerance',0.1,'PopulationSize',200,...
    'InitialPopulationMatrix',init_pop.init_pop);

fprintf('Running optimisation...\n');
[Result.x,Result.fval,Result.exitflag,Result.output,Result.population,Result.scores] = gamultiobj(@objfun,nvar,[],[],[],[],L,U,@confun,options);

save([out_filename '.mat'],'Result');


%% Functions to specify objective and constraints 

% Objective function
    function f = objfun(X)
                
        % First term of cost function is mean of squared differences to measured data
		f1 = mean(((X(iElbow)'-datavec(imeas_elbow))./datasd(imeas_elbow)).^2);
         
        % Second term is mean squared muscle activation
        f2 = mean(X(iact).^2);
           
        f = [f1,f2];
                
    end

% Constraints
    function [c,ceq] = confun(X)
        
        % Linear constraints: ceq(x) = 0
        c = [];
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

            f = das3('Dynamics',(x1+x2)'/2,(x2-x1)'/h,(u1+u2)'/2);
            
            f(lockeddofs) = zeros(length(lockeddofs),1);
            f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
            ceq(iceq) = f;

            % advance the indices
            ix = ix + nvarpernode;
            ius = ius + nvarpernode;
            iceq = iceq + (ncon_eq/(N-1));
        end
    end

end		% end of function das3elbow_optimize_pareto
