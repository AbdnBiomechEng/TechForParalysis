%% Compares the das3_simplified.osim forearm moment arms to the results of the paper https://doi.org/10.1016/j.clinbiomech.2006.03.002

% Load model
das3_modelfile = 'das3_simplified.osim';
das3 = org.opensim.modeling.Model(das3_modelfile);
state = das3.initSystem();
nDofs = das3.getCoordinateSet().getSize();
das3_flex = das3.getCoordinateSet().get('EL_x');
das3_pron = das3.getCoordinateSet().get('PS_y');

% Muscles of interest
musclesinterest = {'supinator_3', 'bic_l', 'bic_b_1', 'brachiorad_2','pron_quad_2','pron_teres_1', 'pron_teres_2'};
nMuscles = numel(musclesinterest);

% Define joint angle combinations to test
angles_flex = [0 45 90];
legend_flex = cellstr(num2str(angles_flex', 'flexion=%-d'));

angles_pron = 5:5:160;
[flex, pron] = meshgrid(angles_flex, angles_pron);
posturecombo = [flex(:), pron(:)];
nPostures = size(posturecombo,1);
nAngles = sqrt(nPostures);

% Allocate moment arm matrix (nPostures x nMuscles)
momentarms = zeros(nPostures, nMuscles);

% Iterate model through posture combo and find moment arms
for posturenum = 1:nPostures
    % Forearm posture
    EL_x = posturecombo(posturenum, 1);
    PS_y = posturecombo(posturenum, 2);

    % Set posture
    das3.updCoordinateSet().get('EL_x').setValue(state, EL_x * pi / 180);
    das3.updCoordinateSet().get('PS_y').setValue(state, PS_y * pi / 180);

    % Get moment arms
    for musclenum = 1:nMuscles
        muscle = char(musclesinterest(musclenum));
        activemuscle = das3.getMuscles().get(muscle);
        momentarms(posturenum, musclenum) = activemuscle.computeMomentArm(state, das3_pron);
    end
end

% Converting from m to mm
momentarms = 1000 * momentarms;

% Visualizing elbow pronation-supination moment arms (+ = pronation)
figure;
for iflex=angles_flex
    indeces = posturecombo(:,1)==iflex;
    xdata = posturecombo(indeces, 2);
    for musclenum = 1:nMuscles
        muscle = char(musclesinterest(musclenum));
        ax = subplot(4, 2, musclenum);
        ydata = momentarms(indeces, musclenum);
        plot(xdata,ydata); hold on;
        xlabel(ax, 'Pronation (deg)', 'FontSize', 8);
        ylabel(ax, 'Moment arm (cm)', 'FontSize', 8);
        title(ax, muscle, 'FontSize', 12);
        ax.FontSize = 8;
    end
end
legend(legend_flex);

return;


% Identify muscles of interest
musclesinterest = {'bic_l', 'bic_b_1', 'tric_long_2', 'tric_med_4', 'brachialis_4', 'brachiorad_2', ...
                   'pron_teres_1', 'pron_teres_2', 'supinator_3', 'pron_quad_2', 'tric_lat_3', 'anconeus_2'};
nMuscles = numel(musclesinterest);

% Allocate moment arm matrix (nPostures x nMuscles x nCoordinates)
momentarms = zeros(nPostures, nMuscles, 2);

% Iterate model through posture combo and find moment arms
for posturenum = 1:nPostures
    test = 1;
    % Forearm posture
    EL_x = posturecombo(posturenum, 1);
    PS_y = posturecombo(posturenum, 2);

    % Set posture
    das3.updCoordinateSet().get('EL_x').setValue(state, EL_x * pi / 180);
    das3.updCoordinateSet().get('PS_y').setValue(state, PS_y * pi / 180);

    % Get moment arms
    for musclenum = 1:nMuscles
        muscle = char(musclesinterest(musclenum));
        activemuscle = das3.getMuscles().get(muscle);
        momentarms(posturenum, musclenum, 1) = activemuscle.computeMomentArm(state, das3_flex);
        momentarms(posturenum, musclenum, 2) = activemuscle.computeMomentArm(state, das3_pron);
    end
end

% Converting from m to cm
momentarms = 100 * momentarms;

% Visualizing elbow flexion-extension moment arms (+ = elbow flexion)
fig = figure;
xdata = posturecombo(:, 1);
ydata = posturecombo(:, 2);
for musclenum = 1:nMuscles
    muscle = char(musclesinterest(musclenum));
    ax = subplot(3, 4, musclenum, 'Parent', fig, 'Projection', 'perspective');

    zdata = momentarms(:, musclenum, 1);
    surf(ax, reshape(xdata, nAngles, nAngles)', reshape(ydata, nAngles, nAngles)', ...
        reshape(zdata, nAngles, nAngles)', 'EdgeColor', 'none', 'FaceColor', 'interp');
    xlabel(ax, 'Flexion (deg)', 'FontSize', 8);
    ylabel(ax, 'Pronation (deg)', 'FontSize', 8);
    zlabel(ax, 'Moment arm (cm)', 'FontSize', 8);
    zlim(ax, [-4, 8]);
    title(ax, muscle, 'FontSize', 12);
    ax.FontSize = 8;
    view(ax, 20, -120);
end

% Visualizing elbow pronation-supination moment arms (+ = pronation)
fig = figure;
xdata = posturecombo(:, 1);
ydata = posturecombo(:, 2);
for musclenum = 1:nMuscles
    muscle = char(musclesinterest(musclenum));
    ax = subplot(3, 4, musclenum, 'Parent', fig, 'Projection', 'perspective');

    zdata = momentarms(:, musclenum, 2);
    surf(ax, reshape(xdata, nAngles, nAngles)', reshape(ydata, nAngles, nAngles)', ...
        reshape(zdata, nAngles, nAngles)', 'EdgeColor', 'none', 'FaceColor', 'interp');
    xlabel(ax, 'Flexion (deg)', 'FontSize', 8);
    ylabel(ax, 'Pronation (deg)', 'FontSize', 8);
    zlabel(ax, 'Moment arm (cm)', 'FontSize', 8);
    zlim(ax, [-4, 8]);
    title(ax, muscle, 'FontSize', 12);
    ax.FontSize = 8;
    view(ax, 20, -120);
end

