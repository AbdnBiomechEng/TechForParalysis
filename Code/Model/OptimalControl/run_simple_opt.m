% In this example, we have a known mass and we try to estimate the
% stiffness of a spring. We attach the mass to the spring in a vacuum (to 
% ensure no energy loss), and give it a push in the horizontal direction.
%
% We record (noisy) kinematic data (we can't get too close because of the
% whole vacuum setup)
%
% Here I am giving you the answer: the frequency of the kinematic data 
% is such that mass = stiffness

mass = 10; % kg <-- change this value to see if the optimisation can find the correct stiffness!
                % Make sure it is within the range of stiffness we search
                % in (line 18 below)

Result.objfun = 10000;  % initial value of the objective function: a very high number

% We go through a range of stiffnesses...
for stiffness = 1:20
    newRes = simple_opt_ipopt('random',200,mass,stiffness); % ... and try to match the measured kinematic data
    if newRes.objfun < Result.objfun  % if the new value of the objective function is lower, this is a good solution, so we keep it
        Result = newRes;
    end
end

% And finally, we plot the results
plot(Result.input_t,Result.input); hold on;
plot(Result.times,Result.x(1,:),'LineWidth',2);
title(['Optimisation result: Spring stiffness = ' num2str(Result.stiffness) ' N/m']);
xlabel('time')
legend('Measured','Optimised');
ylabel('Mass displacement (m)');
