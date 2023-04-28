function Result = das3_optimize_pareto(out_filename)
% This program finds pareto front of three cost functions: muscle effort,
% GH, and scapular stability

modelparams = load('simplemus_model_struct.mat'); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;
nvar = nstates + nmus;              % number of unknowns 

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

% % Further constrain the feasible solution space:
% % Low shoulder elevation
% L(11) = 0;
% U(11) = 30*pi/180;
% % Low elbow flexion
% L(13) = 0;
% U(13) = 50*pi/180;
% % Mid pro/supination
% L(14) = 80*pi/180;
% U(14) = 100*pi/180;
% % low activations
% U(iact) = 0.2;
% U(iexc) = 0.2;

L(lockeddofs) = lockeddofvalues;
U(lockeddofs) = lockeddofvalues;
L(ndof+lockeddofs) = zeros(length(lockeddofs),1);
U(ndof+lockeddofs) = zeros(length(lockeddofs),1);

% Parameters required for glenohumeral constraint
aphi=tand(38.55);
atheta=tand(44.37);
Rgt = glenoid_scap;

% Read partial initial population from previous optimisation runs
init_pop = load('init_pop_pareto');

options = optimoptions('gamultiobj','ConstraintTolerance',0.1,'PopulationSize',200,...
    'InitialPopulationMatrix',init_pop.init_pop);
[Result.x,Result.fval,Result.exitflag,Result.output,Result.population,Result.scores] = gamultiobj(@objfun,nvar,[],[],[],[],L,U,@confun,options);

save([out_filename '.mat'],'Result');

% Objective function
    function f = objfun(X)
        % First term is mean squared muscle activation
        cost1 = mean(X(iact).^2);

        % Second term is thorax-scapula separation
        Fscap = das3('Scapulacontact', X(1:nstates)');
        cost2 = Fscap(1).^2+Fscap(2).^2;
                
        % Third term is glenohumeral instability
        [~, ~, ~, ~, FGH] = das3('Dynamics', X(1:nstates)', zeros(nstates,1),X(iexc)');
        FGHcontact = calculate_FGH(FGH);
        cost3 = mean(FGHcontact.^2);
        
        f = [cost1,cost3];            
            
    end

% Constraints
    function [c,ceq] = confun(X)
        
        c = [];
        
        % Evaluate dynamics
        f = das3('Dynamics',X(1:nstates)',zeros(nstates,1),X(iexc)');
        f(lockeddofs) = zeros(length(lockeddofs),1);
        f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
        ceq = f';  
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
