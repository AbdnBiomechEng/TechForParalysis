function[c, ceq] = nonlinfun(mc, r_mus, Fmax, M_vol, M_stim)
    c = [];
    ceq = (2*(mc(1)^2*r_mus(1)*Fmax(1)+mc(2)^2*r_mus(2)*Fmax(2)+mc(3)^2*r_mus(3)*Fmax(3))*M_stim)/(r_mus(1)*mc(1)*Fmax(1)) - M_vol;
end
