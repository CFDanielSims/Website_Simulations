function mach = M_from_statstag_ratio(r)
    gamma = 1.4;
    syms M;
    ratio = (1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1))==r;

    newmboth = solve(ratio, M, "Real", true);
    mach = double(newmboth(newmboth>0));
end

