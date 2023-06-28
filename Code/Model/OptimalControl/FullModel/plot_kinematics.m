function plot_kinematics(input_file)

%% load data

res = load(input_file);
%modelstruct = load('simplified_model_struct.mat');

nmus = res.Result.nmus;
ndof = res.Result.ndof;

iq = 1:ndof;
iqdot = max(iq) + (1:ndof);
iLce = max(iqdot) + (1:nmus);
iAct = max(iLce) + (1:nmus);

time = res.Result.times;

%% plot kinematics

figure
subplot(5,1,1); plot(time,res.Result.thor_hum*180/pi) % thoraco-humeral angles
ylabel('angles (degrees)'); legend('plane of elevation','angle of elevation','axial rotation');
subplot(5,1,2); plot(time,res.Result.x(4:6,:)*180/pi) % sc joint
ylabel('angles (degrees)'); legend('SC protraction','SC elevation','SC rotation');
subplot(5,1,3); plot(time,res.Result.x(7:9,:)*180/pi) % ac joint
ylabel('angles (degrees)'); legend('AC protraction','AC elevation','AC rotation');
subplot(5,1,4); plot(time,res.Result.x(10:12,:)*180/pi) % gh joint
ylabel('angles (degrees)'); legend('GH plane','GH elevation','GH rotation');
subplot(5,1,5); plot(time,res.Result.x(13:14,:)*180/pi) % elbow
ylabel('angles (degrees)'); xlabel('time (s)'); legend('flexion', 'supination');

%% plot muscles

musclenames = cell(nmus,1);

for imus=1:nmus
    musclenames{imus} = res.Result.muscles{imus}.name;
end

lce = res.Result.x(2*ndof+1:2*ndof+nmus,:);
act = res.Result.x(2*ndof+nmus+1:end,:);
fmus = res.Result.mus_forces;

figure;
subplot(3,1,1)
p1 = plot(time,act);
ylabel('activation');
% dtRow = dataTipTextRow('Muscle', mus_name)
% p1.DataTipTemplate.DataTipRows(end+1) = dtRow;
subplot(3,1,2)
plot(time,fmus)
ylabel('muscle force (N)')
subplot(3,1,3)
plot(time,lce)
ylabel('fibre length (norm)'); xlabel('time (s)');
