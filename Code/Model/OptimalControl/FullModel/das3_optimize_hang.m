function Result = das3_optimize_hang(initialguess,out_filename)
% This program optimizes das3 position to hang with minimum activation

tic;

% optimization settings
solver = 'IPOPT'; 
MaxIterations = 10000;
OptimalityTolerance = 1e-3;
FeasibilityTolerance = 1e-3;

Weffort = 1;        % weight for the energy consumption term in the cost function
Wscap = 0.1;        % weight for scapulo-thoracic gliding plane 
Whum = 0.01;        % weight for glenohumeral stability 

modelparams = load('model_struct_clav_scap_orig.mat'); 
model = modelparams.model;

% Initialize the model
das3('Initialize',model);

ndof = model.nDofs;
nmus = model.nMus;
nstates = 2*ndof + 2*nmus;

% some other constants
nvar = nstates + nmus;        % number of unknowns 
ncon_eq = 2*ndof + nmus;	  % number of constraints due to discretized dynamics
ncon = ncon_eq + nmus;        % we add SEE elongation

% get DOF names
dofnames = cell(ndof,1);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.osim_name;
end

% find muscle names
musclenames = cell(nmus,1);
PEEslack = zeros(nmus,1);
fmax = zeros(nmus,1);

for imus=1:nmus
    musclenames{imus} = model.muscles{imus}.osim_name;
    PEEslack(imus) = model.muscles{imus}.PEEslack;
    fmax(imus) = model.muscles{imus}.fmax;
end

LCEopt = das3('LCEopt');
SEEslack = das3('SEEslack');

% lock thorax dofs
lockeddofs = 1:3;
lockeddofvalues = [0;0;0];

% these are the range of motion limits 
xlim = das3('Limits')';

% Lower and upper bounds for the optimization variables X
% Bounds are:
%   joint angles 	xlim
%   angular velocities 0
%   CE lengths		0.2 to 1.8
%	active states	0 to 1
%   normalised SEE elongation: from -inf to 0.06 -> max force is
%   ~2.25 times the maximum isometric force

max_elong = 0.06;

L = [xlim(1:ndof,1);         % q
    zeros(ndof,1);           % qdot
    zeros(nmus,1) + 0.2;     % Lce
    zeros(nmus,1);           % active states
    -0.01*ones(nmus,1)];     % normalised SEE elongation

U = [xlim(1:ndof,2);         % q
    zeros(ndof,1);           % qdot
    zeros(nmus,1) + 1.8;     % Lce
    1.0*ones(nmus,1);        % active states
    max_elong*ones(nmus,1)]; % normalised SEE elongation

L(lockeddofs) = lockeddofvalues;
U(lockeddofs) = lockeddofvalues;
L(ndof+lockeddofs) = zeros(length(lockeddofs),1);
U(ndof+lockeddofs) = zeros(length(lockeddofs),1);

objeval = 0;
coneval = 0;
objlog = [];
conlog = [];
print_flag = 0;

% define indices to the state variables within the state vector x
iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iact = max(iLce) + (1:nmus);

% precompute the indices for muscles that cross GH
iGH = [];
for imus=1:nmus
    if(das3('crossGH', imus))
        iGH=[iGH imus];
    end
end

% make an initial guess
if strcmp(initialguess, 'mid')
    X0 = (L + U)/2;						% halfway between upper and lower bound
elseif strcmp(initialguess, 'random')
    X0 = L + (U - L).*rand(size(L));	% random between upper and lower bound
elseif strcmp(initialguess, 'dsem')
    X0 = zeros(nstates,1);					
    %X0(4:14) = [-0.3802;0.11;0;0.808;0.0855;-0.0206;0;0;0;0;1.3963]; % initial DSEM position
    X0(4:14) = [-0.3802;0.11;0;0.8727;0.0855;-0.0206;0;0;0;0;1.3963]; % initial DSEM position - increased ACz
    X0(iLce) = 1;
    init_lengths = das3('Musclelengths',X0);
    X0(nstates+(1:nmus)) = 0; % SEE elongation
    X0(iLce) = (init_lengths - SEEslack)./LCEopt;
else
    % load a previous solution, initialguess contains file name
    initg = load(initialguess);
    X0 = initg.Result.X;
end

% for glenohumeral constraint:
aphi=tand(38.55);
atheta=tand(44.37);
Rgt = glenoid_scap;

% run optimization
if (MaxIterations > 0)
    print_flag=0;    
    % determine sparsity structure of Jacobian (midpoint discretization)
    % I have verified that we always get same structure by testing with random X
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
        X = L + (U-L).*rand(size(L));
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
        figure; spy(cjac_num);
        cjac_full = full(cjac);
        % find the max difference in constraint jacobian and objective gradient
        [maxerr,irow] = max(abs(cjac-cjac_num));
        [maxerr,icol] = max(maxerr);
        fprintf('Max.error in constraint jacobian: %8.5f at %d %d, %8.8f percent error\n', maxerr, irow(icol), icol, maxerr/cjac_full(irow(icol),icol));
        [maxerr,irow] = max(abs(grad_num-grad));
        fprintf('Max.error in objective gradient:  %8.5f at %d\n', maxerr, irow);
        keyboard
    end
    
    % evaluate initial guess
    print_flag=1;
    objfun(X0);
    confun(X0);
    print_flag=0;
    
    % set up constraints
    cl = zeros(ncon,1);
    cu = zeros(ncon,1);

    funcs.objective = @objfun;
    funcs.gradient  = @objgrad;
    
    options.cl = cl;
    options.cu = cu;
    funcs.constraints = @confun;
    funcs.jacobian    = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    
    options.lb = L;
    options.ub = U;

    % IPOPT options
    options.ipopt.print_level = 0;
    options.ipopt.max_iter = MaxIterations;
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy = 'adaptive';		% worked better than 'monotone'
    options.ipopt.bound_frac = 0.001;			% worked better than 0.01 or 0.0001
    options.ipopt.bound_push = options.ipopt.bound_frac;
    options.ipopt.tol = OptimalityTolerance;
    options.ipopt.acceptable_constr_viol_tol = FeasibilityTolerance;
    options.ipopt.acceptable_tol = FeasibilityTolerance;
    options.ipopt.constr_viol_tol = FeasibilityTolerance;
    options.ipopt.print_info_string = 'yes';
    options.ipopt.limited_memory_max_history = 12;	% 6 is default, 12 converges better, but may cause "insufficient memory" error when N is large
    options.ipopt.dual_inf_tol = 1;

    [X, info] = ipopt(X0,funcs,options);
    Result.status = info.status;
    fprintf('IPOPT info status: %d\n', info.status);
    Result.status = info.status;
    Result.info = info;
    Result.X = X;
    Result.X0 = X0;
    Result.solver = solver;

else		% skip optimization
    Result.X = X0;
end

Result.duration = toc;
disp(['Time elapsted: ' num2str(Result.duration) ' seconds']);
beep;

obj_str = '';
print_flag=1;
objfun(Result.X);
confun(Result.X);
print_flag=0;

% save this result on a file
Result.ndof = ndof;
Result.nmus = nmus;
Result.model = model;
Result.obj_str = obj_str;
Result.x = Result.X(1:nstates);
Result.MaxIterations = MaxIterations;
Result.OptimalityTolerance = OptimalityTolerance;
Result.FeasibilityTolerance = FeasibilityTolerance;
Result.times = 0.1;
Result.u = Result.X(iact);
Result.elong = Result.X(nstates+(1:nmus));
angles = Result.X(1:ndof);

% Calculate muscle lengths, forces, GH and scapula stability
Result.mus_lengths = das3('Musclelengths', Result.x);
Result.mus_forces = das3('Muscleforces', Result.x);
[f, ~, ~, ~, FGH, ~, Result.thor_hum] = das3('Dynamics', Result.x, 0*Result.x,zeros(nmus,1));
Result.f = f;
Result.Fscap = das3('Scapulacontact', Result.x)';
Result.FGHcontact = calculate_FGH(FGH);

fprintf('\n\nDOF               angle(deg)    limits (deg)            ang.vel(deg/s)   moment(Nm)  ext moment(Nm) \n');
fprintf('--------------- --------------  ---------------------   -------------- -----------  --------------\n');
moments = das3('Jointmoments',Result.x);
ext_moments = f(ndof+1:2*ndof) - moments;
for i=1:ndof
    fprintf('%-15s %9.3f      %9.3f   %9.3f      %9.3f    %9.3f    %9.3f\n',dofnames{i}, 180/pi*Result.x(i), 180/pi*xlim(i,1), 180/pi*xlim(i,2), 180/pi*Result.x(ndof+i), moments(i), ext_moments(i));
end


fprintf('\n\nMuscle           Lce/Lceopt   PEEslack    SEE elong     Muscletendon length    Activation    Force(N)  \n');
fprintf('--------------- ------------ ----------- -------------- ---------------------  ----------   ----------\n');
for i=1:nmus
    fprintf('%-15s %9.3f    %9.3f    %9.3f      %9.3f           %9.3f     %9.3f\n',musclenames{i}, Result.x(2*ndof+i), PEEslack(i), Result.elong(i), Result.mus_lengths(i), Result.x(2*ndof+nmus+i), Result.mus_forces(i));
end

save([out_filename '.mat'],'Result');
make_osimm(out_filename,dofnames,[angles,angles]); % two rows so Opensim 4.x can open it

% objective function
    function f = objfun(X)
        
        % First term is mean squared muscle activations
        wf1 = mean(X(iact).^2);
        f1 = Weffort * wf1;
        f = f1;
        
        if Wscap
            % Second term is thorax-scapula constraint
            Fscap = das3('Scapulacontact', X(1:nstates));

            wf2 = Fscap(1).^2+Fscap(2).^2;
            f2 = Wscap * wf2;
            f = f + f2;
        end
                
        if Whum
            % Third term is glenohumeral stability constraint
            [~, ~, ~, ~, FGH] = das3('Dynamics', X(1:nstates), zeros(nstates,1),X(iact));
            FGHcontact = calculate_FGH(FGH);

            wf3 = mean(FGHcontact.^2);
            f3 = Whum * wf3;
            f = f + f3;
        end
        
        obj_str = sprintf('Objfun (weighted): %9.5f (effort)', f);
        wobj_str = sprintf('Objfun (unweighted): %9.5f (effort)', wf1);

        if print_flag
            fprintf(obj_str);
            fprintf('\n');
            fprintf(wobj_str);
            fprintf('\n');
        end
            
    end


% gradient of objective function
    function g = objgrad(X)
        
        % Perturbation to estimate derivative
        dx = 1e-7;

        g = zeros(nvar,1);
        
        % First term is mean squared muscle activation and force
        g(iact) = g(iact) + 2*Weffort*(X(iact))/nmus;
        
        if Wscap
            % Second is the scapula term
            Fscap = das3('Scapulacontact', X(1:nstates));
            % Scapulo-thoracic glinding plane only depends on dofs 4-9
            % (SC_y,SC_z,SC_x,AC_y,AC_z,AC_x)
            for istate = 4:9
                xisave = X(istate);
                X(istate) = X(istate) + dx;
                Fscap_dx = das3('Scapulacontact', X(1:nstates));
                dFscap = (Fscap_dx-Fscap)/dx;
                g(istate) = g(istate) + Wscap*(2*Fscap(1)*dFscap(1) + 2*Fscap(2)*dFscap(2));
                X(istate) = xisave;
            end
        end
        
        if Whum
            % Third is the humerus term
            [~, ~, ~, ~, FGH] = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iact));

            % calculate GH constraint
            FGHcontact = calculate_FGH(FGH);

            % Perturb x1 and x2 by dx to estimate derivative
            for istate = [1:2*ndof 2*ndof+iGH 2*ndof+nmus+iGH]
                xisave = X(istate);
                X(istate) = X(istate) + dx;

                % estimate delta(GHconstraint)/delta(x)
                [~, ~, ~, ~, FGH] = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iact));
                FGHcontact_dx = calculate_FGH(FGH);
                dFGHcontact = (FGHcontact_dx-FGHcontact)/dx;
                g(istate) = g(istate) + Whum * 2 * FGHcontact * dFGHcontact;
                X(istate) = xisave;
            end     
        end
    end


% constraints
    function allcon = confun(X)
 
        % evaluate dynamics
        f = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iact));
        f(lockeddofs) = zeros(length(lockeddofs),1);
        f(ndof+lockeddofs) = zeros(length(lockeddofs),1);
        ceq(1:ncon_eq) = f(1:2*ndof+nmus)';        
                
        % calculate SEE elongation 
        lengths = das3('Musclelengths',X(1:nstates));
        ceq(ncon_eq+(1:nmus)) = X(nstates+(1:nmus)) - (lengths - X(iLce).*LCEopt - SEEslack)./SEEslack;            
                    
        allcon = ceq;
        
        if (print_flag)
            fprintf('Norm(ceq): %9.5f  \n', norm(ceq));
        end
    end

    function J = conjacstructure()
        % returns structure of constraint Jacobian matrix
        J = Jpattern;
    end


% Jacobian of constraints
    function J = conjac(X)
        
        % Perturbation to estimate derivative
        dx = 1e-7;
        
        allrows = zeros(Jnnz,1);
        allcols = zeros(Jnnz,1);
        allvals = zeros(Jnnz,1);
        index=1;
        
        % evaluate dynamics         
        [~, dfdx, ~, ~] = das3('Dynamics',X(1:nstates),zeros(nstates,1),X(iact));
        dfdx(lockeddofs,:) = zeros(length(lockeddofs),nstates);
        dfdx(ndof+lockeddofs,:) = zeros(length(lockeddofs),nstates);
                      
        % which generates one block in the Jacobian:
        [r,c,v] = find(dfdx(1:2*ndof+nmus,:));
        datal = length(v);
        allrows(index:index+datal-1) = r;
        allcols(index:index+datal-1) = c;
        allvals(index:index+datal-1) = v;
        index = index+datal;
                    
        % SEE elongation
        iceq_elong = ncon_eq+(1:nmus);

        elong = X(nstates+(1:nmus));
        lengths = das3('Musclelengths',X(1:nstates));
        elong_diff = elong - (lengths - X(iLce).*LCEopt - SEEslack)./SEEslack;

        % Perturb x1 to estimate derivative
        for istate = 1:nstates
            xisave = X(istate);
            X(istate) = X(istate) + dx;

            % estimate delta(elong_diff)/delta(x)
            lengths = das3('Musclelengths',X(1:nstates));
            elong_diff_dx = elong - (lengths - X(iLce).*LCEopt - SEEslack)./SEEslack;

            allrows(index:index+nmus-1) = iceq_elong;
            allcols(index:index+nmus-1) = istate;
            allvals(index:index+nmus-1) = (elong_diff_dx-elong_diff)/dx;
            index=index+nmus;
            X(istate) = xisave;
        end

        allrows(index:index+nmus-1) = iceq_elong;
        allcols(index:index+nmus-1) = nstates+(1:nmus);
        allvals(index:index+nmus-1) = ones(nmus,1);
        index=index+nmus;
                              
        if exist('Jpattern','var')
            [r,c,~] = find(Jpattern);
            nz_rows = ismember([allrows(1:index-1),allcols(1:index-1)],[r,c], 'rows');
            J = sparse(allrows(nz_rows),allcols(nz_rows),allvals(nz_rows),ncon, nvar);
        else
            J = sparse(allrows(1:index-1),allcols(1:index-1),allvals(1:index-1),ncon, nvar);
        end
               
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


end		% end of function das3_optimize_hang
