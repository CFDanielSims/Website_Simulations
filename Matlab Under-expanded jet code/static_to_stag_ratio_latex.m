function show = static_to_stag_ratio_latex()
    gamma = 1.4;
    syms M;
    ratio = (1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1));

    show = latex(ratio);