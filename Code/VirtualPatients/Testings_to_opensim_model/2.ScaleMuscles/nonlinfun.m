function[c, ceq_new] = nonlinfun(mc, r_mus, Fmax, M_vol, M_stim)

    c = [];
    
    % original one
    ceq = (2*      (mc(1)^2*r_mus(1)*Fmax(1)+mc(2)^2*r_mus(2)*Fmax(2)+mc(3)^2*r_mus(3)*Fmax(3))     *M_stim)/   (   r_mus(1)*mc(1)*Fmax(1)   ) - M_vol;
    
    % Do not restrict to 3 muscles
    momentOfAllMuscles = 0;
    for x = 1:numel(mc); % matlab starts the loops at 1 and not 0!
        momentOfAllMuscles = momentOfAllMuscles + (mc(x)^2*r_mus(x)*Fmax(x));
    end

    momentOfBicepsMuscles = r_mus(1)*mc(1)*Fmax(1) + r_mus(2)*mc(2)*Fmax(2) + r_mus(3)*mc(3)*Fmax(3)

    ceq_new = (2 * momentOfAllMuscles * M_stim) / momentOfBicepsMuscles - M_vol;

end

