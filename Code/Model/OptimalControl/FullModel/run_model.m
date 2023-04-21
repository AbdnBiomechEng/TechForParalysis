% For dynamic simulation.  Advances the model by one time step.
    %
	% Inputs
	%	x		(298 x 1) System state
	%	u		(138 x 1) Muscle excitations, assumed constant during the time step
	%	h		(scalar)  Time step size
    %   xdot    (298 x 1) System state derivatives (initially zero)
	%	step_u	(138 x 1) Muscle excitations from previous step (initially zero)
    % Optional inputs
    % 	M		(5 x 1)	  Moments applied to the thorax-humerus YZY axes and the elbow flexion and supination axes
	%   exF     (2 x 1)   Vertical force of amplitude exF(2) applied to the ulna at a distance of exF(1) from the elbow 
	%						(to simulate a mobile arm support)
    %   handF   (3 x 1)   Force at the CoM of the hand
    %
	% Outputs
	%	xnew	(298 x 1) The system state at the end of the time step
	%	FGH		(3 x 1)   The glenohumeral reaction force, acting on the scapula


modelparams = load('simplemus_model_struct.mat'); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% Initialize the model
das3('Initialize',model);

x_eq = load('equilibrium');
g = das3('Dynamics',x_eq.x,zeros(nstates,1),zeros(nmus,1));

u = zeros(nmus,1); 
step_u = zeros(nmus,1);
h = 0.0001;
xdot = zeros(nstates,1);
    
[xnew, xdot, step_u, FGH] = das3step(x, u, h, xdot, step_u);