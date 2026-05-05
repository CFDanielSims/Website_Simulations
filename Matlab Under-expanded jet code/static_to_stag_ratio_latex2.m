function static_to_stag_ratio_latex2()
    M = sym('M'); 
    gamma = sym('gamma');
    expr = (1+((gamma-1)/2)*M.^2).^(-gamma/(gamma-1));
    str  = latex(expr);

    figure('Color','w');      % white background
    ax = axes;
    text(0.05, 0.5, ['$' str '$'], ...
        'Interpreter','latex', ...
        'FontSize', 50, ...
        'Color','k', ...       % <-- force black text
        'Parent', ax);
    axis off
end
