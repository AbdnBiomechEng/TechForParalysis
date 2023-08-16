function opt_output_analysis(filename)
% Function to plot various outputs of an optimal control simulation
% Input is the filename (with each relative path)

% Load results and model
res_sim = load(filename);
res_model = load(res_sim.Result.OptSetup.model);
model = res_model.model;

% Initialise model to calculate additional parameters
das3('Initialize',model);

nmus = model.nMus;
ndof = model.nDofs;
N = res_sim.Result.OptSetup.N;

% Get muscle names and parameters
musclenames = cell(nmus,1);
PEEslack = zeros(nmus,1);
SEEslack = zeros(nmus,1);
LCEopt = zeros(nmus,1);
for imus=1:nmus
    musclenames{imus} = model.muscles{imus}.name;
    PEEslack(imus) = model.muscles{imus}.PEEslack;
    SEEslack(imus) = model.muscles{imus}.lslack;
    LCEopt(imus) = model.muscles{imus}.lceopt;
end

% Get DOF names and range
dofnames = cell(ndof,1);
range = zeros(ndof,2);
for idof=1:ndof
    dofnames{idof} = model.dofs{idof}.name;
    range(idof,:) = model.dofs{idof}.range;
end

time = res_sim.Result.times;
lce_sim = res_sim.Result.x(2*ndof+1:2*ndof+nmus,:);
act_sim = res_sim.Result.x(2*ndof+nmus+1:end,:);
fmus_sim = res_sim.Result.mus_forces;

% Calculate SEE elongation, moments and momentarms for each node
SEE_elong = zeros(nmus,N);
moments = zeros(ndof,N);
momentarms = zeros(nmus,ndof,N);
for inode = 1:N
    SEE_elong(:,inode) = res_sim.Result.mus_lengths(:,inode) - res_sim.Result.x(2*ndof+(1:nmus),inode).*LCEopt - SEEslack;
    moments(:,inode) = das3('Jointmoments', res_sim.Result.x(:,inode));
    momentarms(:,:,inode) = full(das3('Momentarms', res_sim.Result.x(:,inode)));
end

%% Plot muscle variables
fh1 = figure;
subplot(4,1,1); plot(time,act_sim','o-'); 
title('Muscle activations');
legend(musclenames); legend off
subplot(4,1,2); plot(time,lce_sim','o-');
title('Normalised fibre lengths');
legend(musclenames); legend off
subplot(4,1,3); plot(time,SEE_elong','o-');
title('SEE elongation');
legend(musclenames); legend off
subplot(4,1,4); plot(time,fmus_sim','o-'); 
title('Muscle forces (N)'); xlabel('time (s)');
legend(musclenames); legend off
datacursormode on;
dcm = datacursormode(fh1);
set(dcm,'UpdateFcn',@customdatatip)

%% Plot angles (two different plots for input and output angles)
% figure;
% subplot(5,1,1); plot(time,res_sim.Result.thor_hum*180/pi,'o-') % thoraco-humeral angles
% ylabel('angles (degrees)'); legend('plane of elevation','angle of elevation','axial rotation');
% title('Optimisation output angles');
% subplot(5,1,2); plot(time,res_sim.Result.x(4:6,:)*180/pi,'o-') % sc joint
% ylabel('angles (degrees)'); legend('SC protraction','SC elevation','SC rotation');
% subplot(5,1,3); plot(time,res_sim.Result.x(7:9,:)*180/pi,'o-') % ac joint
% ylabel('angles (degrees)'); legend('AC protraction','AC elevation','AC rotation');
% subplot(5,1,4); plot(time,res_sim.Result.x(10:12,:)*180/pi,'o-') % gh joint
% ylabel('angles (degrees)'); legend('GH plane','GH elevation','GH rotation');
% subplot(5,1,5); plot(time,res_sim.Result.x(13:14,:)*180/pi,'o-') % elbow
% ylabel('angles (degrees)'); xlabel('time (s)'); legend('flexion', 'supination');
% 
input_angles = res_sim.Result.resampled_data';
% figure;
% subplot(4,1,1); plot(time,input_angles(1:3,:)*180/pi,'o-') % sc joint
% ylabel('angles (degrees)'); legend('SC protraction','SC elevation','SC rotation');
% title('Optimisation input angles');
% subplot(4,1,2); plot(time,input_angles(4:6,:)*180/pi,'o-') % ac joint
% ylabel('angles (degrees)'); legend('AC protraction','AC elevation','AC rotation');
% subplot(4,1,3); plot(time,input_angles(7:9,:)*180/pi,'o-') % gh joint
% ylabel('angles (degrees)'); legend('GH plane','GH elevation','GH rotation');
% subplot(4,1,4); plot(time,input_angles(10:11,:)*180/pi,'o-') % elbow
% ylabel('angles (degrees)'); xlabel('time (s)'); legend('flexion', 'supination');

%% Compare shoulder angles on same plot
fig = figure; colororder(fig,[0.8 0 0.5; 0 0.6 0; 0 0 0.6]);
subplot(4,1,1);
    plot(time,input_angles(1:3,:)*180/pi,'o-'); hold on; % sc joint
    plot(time,res_sim.Result.x(4:6,:)*180/pi,'x--');
    ylabel('angles (degrees)'); legend('SC protraction','SC elevation','SC rotation');
subplot(4,1,2);
    plot(time,input_angles(4:6,:)*180/pi,'o-'); hold on % ac joint
    plot(time,res_sim.Result.x(7:9,:)*180/pi,'x--');
    ylabel('angles (degrees)'); legend('AC protraction','AC elevation','AC rotation');
subplot(4,1,3);
    plot(time,input_angles(7:9,:)*180/pi,'o-'); hold on; %gh joint
    plot(time,res_sim.Result.x(10:12,:)*180/pi,'x--');
    ylabel('angles (degrees)'); legend('GH plane','GH elevation','GH rotation');
    
fig = figure; colororder(fig,[0.8 0 0.5; 0 0 0.6]);
    plot(time,input_angles(10:11,:)*180/pi,'o-'); hold on % elbow
    plot(time,res_sim.Result.x(13:14,:)*180/pi,'x--')
    ylabel('angles (degrees)'); legend('flexion', 'supination');
    xlabel('time (s)');

%% Plot muscle moments
fh2 = figure;
for idof=4:9
    subplot(6,1,idof-3);
    momarms_dof = squeeze(momentarms(:,idof,:));
    moments_dof = momarms_dof.*fmus_sim;
    plot(moments_dof','o-');
    ylabel(dofnames{idof});
    if idof==4, title('Moment (Nm)'); end
    legend(musclenames);
    legend off
end
datacursormode on;
dcm = datacursormode(fh2);
set(dcm,'UpdateFcn',@customdatatip)

fh3 = figure; % plot moments for humeral and forearm dofs
for idof=10:14
    subplot(5,1,idof-9);
    momarms_dof = squeeze(momentarms(:,idof,:));
    moments_dof = momarms_dof.*fmus_sim;
    plot(moments_dof','o-');
    ylabel(dofnames{idof});
    if idof==10, title('Moment (Nm)'); end
    legend(musclenames);
    legend off
end
datacursormode on;
dcm = datacursormode(fh3);
set(dcm,'UpdateFcn',@customdatatip)


function output_txt = customdatatip(obj,event_obj,str)
pos = get(event_obj, 'Position');
output_txt = {...
    ['X: ', num2str(pos(1),4)]...
    ['Y: ', num2str(pos(2),4)] ...
    ['muscle: ', event_obj.Target.DisplayName]...
};

