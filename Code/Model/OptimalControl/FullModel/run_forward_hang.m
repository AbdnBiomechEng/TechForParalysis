function run_forward_hang(filename)
% Runs a forward dynamic simulation from the equilibrium posture 

% Load results and model
res_sim = load(filename);
model = res_sim.Result.model;

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
tend = 8;
tstep = .001;
nsteps = round((tend-t)/tstep);

% reserve space to store results
tout = tstep*(0:nsteps)';
xout = zeros(nsteps+1, nstates);
tout(1) = t;

% set state x and excitations u
x = res_sim.Result.X(1:nstates);
%u = res_sim.Result.u; % from the input file
u = zeros(nmus,1); % zero for passive equilibrium

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
make_osimm('hang_forward',dofnames,xout(:,1:ndof),tout);
save hang_forward.mat xout

% save final position as initial for simulation
LCEopt = das3('LCEopt');
SEEslack = das3('SEEslack');

X = xout(end,:)';
lengths = das3('Musclelengths',X);
% SEE elongation 
X(nstates+(1:nmus)) = (lengths - X(2*ndof + (1:nmus)).*LCEopt - SEEslack)./SEEslack;
Result.X = X;
Result.model = model;
save hang_forward_init Result


