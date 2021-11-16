%% Muscle parameter optimisation
% The elbow is actively flexed by the (1) biceps, (2) brachialis and (3) brachioradialis 
% muscles. The resultant moment is measured and we want to scale the muscle forces 
% to the measured data.
% 
% Script uses Matlab's built-in constrained optimisation routine 'fmincon'. 
% Type 'help fmincon' in the command window for details.
% 
% Define the measured moments (you can change this):

M_vol = 8.75; % measured voluntary moment in Nm
M_stim = 8.75; % moment achieved under stimulation of 1 muscle (biceps)
%% 
% Moment arms of the muscles are as follows:

r_mus = [0.05 0.03 0.07]; % moment arms in metres
Fmax  = [842 999 162]; % the maximum isometric force of each muscle in N.
%% 
% Assume for our individual, max muscle activation is limited to mc (max_control) 
% [0...1], and the value of Fmax is reduced by the same amount again due to atrophy. 
% There is an additional scaling factor applied uniformly to all the muscles simply 
% as a result of the size of the person: this is 'scale'. 
% 
% Therefore:
% 
% $$M_{vol} = \sum{{mc_i}^2.r_{musi}.F_{maxi}.scale$$
% 
% $$M_{stim} = mc_1.r_{mus1}.F_{max1}.scale$$
% 
% Hence $scale = \frac{M_{stim}}{mc_1.r_{mus1}.F_{max1}}$
% 
% and $M_{vol} = \sum{\frac{{{mc_i}^2}.r_{musi}.F_{maxi}.M_{stim}}{mc_1.r_{mus1}.F_{max1}}}$     
% (this is our non-linear equality constraint)
% 
% The goal is to minimise the function 'fun', which is:
% 
% $${(mc_1-mc_2)}^2 + {(mc_2-mc_3)}^2$$

fun = @(mc) ((mc(1)-mc(2))^2 + (mc(2)-mc(3))^2);  % anonymous function describing objective function
nonlcon = @(mc) nonlinfun(mc,r_mus,Fmax,M_vol,M_stim); % anonymous function for non-linear constraints with additional arguments

[mc, J, exitflag, output] = fmincon(fun, ... % 'fun' is the name of the objective function
    [0.5 0.5 0.5], ...  % Initial guess, or starting point
    [], [], ... % there are no linear inequality contraints
    [], [], ... % these are A and B in A[x] = B, the linear equality contraint
    [0 0 0], [1 1 1], ...  % lower and upper bounds on the values
    nonlcon,... % non-linear constraints
    []); % options
%% 
% Check results:
scale = 2*M_stim/(r_mus(1)*mc(1)*Fmax(1))
M_vol_opt = scale*(mc(1)^2*r_mus(1)*Fmax(1) + mc(2)^2*r_mus(2)*Fmax(2) + mc(3)^2*r_mus(3)*Fmax(3)); % M_vol
M_stim_opt = 0.5*scale*(r_mus(1)*mc(1)*Fmax(1)); % M_stim
% Error = r_mus * F_mus' - M_mus
output.message
fprintf('mc = %f\n', mc);
fprintf('M_vol_opt = %f\n', M_vol_opt);
fprintf('M_stim_opt = %f\n', M_stim_opt);