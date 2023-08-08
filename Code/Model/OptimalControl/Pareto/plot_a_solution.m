function plot_a_solution(filename)
% Create Opensim motion file from one of the solutions in the Pareto front

fileload = load(filename);
one_solution = fileload.Result.x(1,:);
x = reshape(one_solution,[],Result.OptSetup.N);
angles = x(1:14,:)';
dofnames = {'TH_x','TH_y','TH_z','SC_y','SC_z','SC_x','AC_y','AC_z','AC_x','GH_y','GH_z','GH_yy','EL_x','PS_y'};
make_osimm(filename,dofnames,angles);
