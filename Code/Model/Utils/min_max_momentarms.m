function min_max_momentarms(model_struct,posture)

% Prints out the maximum and minimum moment arms in a posture of interest
% The posture is a 1x14 vector of angles in radians

if nargin<2 
    posture = [0 0 0 -0.38 0.19 0.23 0.64 -0.38 -0.29 0.42 0.30 0.27 0.39 1.49];
end

if ~isequal(size(posture), [1 14])
    error('The posture should be a 1x14 vector of angles (in radians)');
end

modelstr = load(model_struct);
model = modelstr.model;

% Get muscle names
nmus = model.nMus;
musclenames = cell(nmus,1);
for imus=1:nmus
    musclenames{imus} = model.muscles{imus}.name;
end

% Get dof names
ndof = model.nDofs;
dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% Initialize the model
das3('Initialize',model);

% Place it in a posture of interest
nstates = 2*ndof + 2*nmus;
x = zeros(nstates,1);
x(1:14) = posture;

% Get moment arms (in sparse matrix)
MA = das3('Momentarms', x);
MA = full(MA); 

% for each DOF, report the largest positive and negative moment arms
fprintf('DOF   Largest positive moment arm          Largest negative moment arm\n')
fprintf('----  -----------------------------        -----------------------------------\n');
MA = full(MA);
for i = 1:ndof
    [dmin, imin] = min(MA(:,i));
    [dmax, imax] = max(MA(:,i));
    name_min = musclenames{imin};
    name_max = musclenames{imax};
    if dmin==0 && dmax==0
        fprintf('%-15s\n',dofnames{i});			% there are no muscles spanning this joint
    else
        fprintf('%-15s  %6.1f mm (%-16s)         %6.1f mm (%-16s)\n',dofnames{i},1000*dmax,name_max,1000*dmin,name_min);
    end
end
