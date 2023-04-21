modelparams = load('simplemus_model_struct.mat'); 
model = modelparams.model;

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% Lock thorax dofs
lockeddofs = [1,2,3];

% Precompute the indices for excitations
iexc = nstates+(1:nmus);

% Initialize the model
das3('Initialize',model);

init_pop = [];

for i=1:633
    one_file = load(['C:\Users\s04db9\OneDrive - University of Aberdeen\Research\Tech4Paralysis\ElbowOnly\pareto_con\random_feas_',num2str(i),'.mat']);
    X = one_file.Result.X;
    
    % Evaluate dynamics
    f = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iexc));
    f(lockeddofs) = zeros(length(lockeddofs),1);
    f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
    ceq = f';  
    if sum(ceq)>0.001, disp('Constraints not satisfied');
    else
        init_pop = [init_pop; X'];
    end
    if size(init_pop,1)==200, break; end
end