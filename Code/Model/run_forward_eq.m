function run_forward_eq(filename)
% Runs a forward dynamic simulation from a previous solution for 3
% seconds

% Load results and model
res_sim = load(filename);
model = res_sim.Result.model;

% eq = load(filename);
% x = eq.x;

% modelstr = load('extended_workspace_model.mat');
% model = modelstr.model;

for imus=[37:40 48:51 95:102]  % posterior deltoid, anterior deltoid, pec major
    model.muscles{imus}.PEEslack = 1.2;
end

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);

% Initialize the model
das3('Initialize',model);

LCEopt = das3('LCEopt');
SEEslack = das3('SEEslack');

% lock thorax dofs
lockeddofs = 1:3;
lockeddofvalues = zeros(length(lockeddofs),1);
ix = setdiff(1:nstates, [lockeddofs  ndof+lockeddofs]);
nstates_ix = length(ix);  % number of states in the reduced-dof system

% set simulation parameters
t = 0;
tend = 10;
tstep = .001;
nsteps = round((tend-t)/tstep);

% reserve space to store results
tout = tstep*(0:nsteps)';
xout = zeros(nsteps+1, nstates);
SEEelong = zeros(nmus,nsteps+1);
tout(1) = t;

% load state x
x = res_sim.Result.x(:,end);
u = res_sim.Result.u(:,end);

u([83:88 103:137]) = 0.0001;
% u([83:137]) = 0.0001;
% u(3:15) = 0.0001; % trapezius + levator scapulae
% u(29:36) = 0.0001; % serratus
% u(37:40) = 0.0001; % posterior deltoid
% u(48:51) = 0.0001; % anterior deltoid

% run simulation
xout(1,:) = x';
lengths = das3('Musclelengths', x);
SEEelong(:,1) = (lengths - SEEslack - x(iLce).* LCEopt)./SEEslack;

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

    lengths = das3('Musclelengths', xfull);
    SEEelong(:,i+1) = (lengths - SEEslack - xfull(iLce).* LCEopt)./SEEslack;

    t = t + tstep;
end
simtime = toc;

% report computation time on the screen
fprintf('CPU time per time step: %8.3f ms\n', 1000*simtime/nsteps);
fprintf('Simulation speed is %8.3f times real time\n',tend/simtime);
    
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};
make_osimm([filename '_forward'],dofnames,xout(:,1:ndof),tout);
save([filename '_forward.mat'], 'xout','u','SEEelong');

