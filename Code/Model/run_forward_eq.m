% Runs a forward dynamic simulation from the equilibrium posture for 3
% seconds

modelparams = load('simplified_model_struct.mat'); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% Initialize the model
das3('Initialize',model);

% lock thorax dofs
lockeddofs = 1:3;
lockeddofvalues = zeros(length(lockeddofs),1);
ix = setdiff(1:nstates, [lockeddofs  ndof+lockeddofs]);
nstates_ix = length(ix);  % number of states in the reduced-dof system

% set simulation parameters
t = 0;
tend = 2;
tstep = .001;
nsteps = round((tend-t)/tstep);

% reserve space to store results
tout = tstep*(0:nsteps)';
xout = zeros(nsteps+1, nstates);
tout(1) = t;

% load equilbrium state x
x_eq = load('equilibrium');
x = x_eq.Result.x;
u = x_eq.Result.u;

% Selectively change some muscle activations

%u(1:3) = 0.2; % traps
%u(7:9) = 0.2; % Serr ant
%u(10:12) = 0.4; % Delts
%u(29) = 0.2; % brd

% run simulation
xout(1,:) = x';
x = x(ix);    
% Initialize variables
step_xdot = zeros(nstates_ix,1);
step_u = u;
y = zeros(nstates,1);               % state vector for full system
y(lockeddofs) = lockeddofvalues;    % insert the joint angle values for the locked joints
ydot = zeros(nstates,1);            % state vector derivative for full system   
xfull = zeros(nstates,1);           % for storage
xfull(lockeddofs) = lockeddofvalues;        
tic
for i=1:nsteps
%    u = stimfun(t);

    % Evaluate dynamics in current x and xdot (states of the reduced-dof system)
    y(ix) = x;              % put the unlocked state variables in their proper place inside y
    ydot(ix) = step_xdot;   % same for ydot
    % compute dynamics residuals for the full model, and ignore those that correspond to locked joints
    [g, dgdy, dgdydot, dgdu] = das3('Dynamics',y, ydot, u, zeros(5,1),[0;0],[0;0;0]);       
    f = g(ix);
    dfdx = dgdy(ix,ix);
    dfdxdot = dgdydot(ix,ix);
    dfdu = dgdu(ix,:);

    % Solve the change in x from the 1st order Rosenbrock formula
    du = u - step_u;
    dx = (dfdx + dfdxdot/tstep)\(dfdxdot*step_xdot - f - dfdu*du);
    x = x + dx;

    % Update variables for the next simulation step
    step_xdot = dx/tstep;
    step_u = u;

    xfull(ix) = x;          % put the unlocked state variables in their proper place inside x
    xout(i+1,:) = xfull';   % store result

    t = t + tstep;
end
simtime = toc;

% report computation time on the screen
fprintf('CPU time per time step: %8.3f ms\n', 1000*simtime/nsteps);
fprintf('Simulation speed is %8.3f times real time\n',tend/simtime);
    
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};
make_osimm('equilibrium_forward',dofnames,xout(:,1:ndof),tout);
save equilibrium_forward.mat xout

