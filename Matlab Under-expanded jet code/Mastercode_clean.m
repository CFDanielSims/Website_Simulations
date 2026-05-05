%%Mastercode clean version
%%Now use symmetry so you only have to compute one expansion fan
xlimits = [0 8];     % set your x range
ylimits = [0 1.5];   % set your y range (this fixes plot height)
machplot     = false;   % true = show Mach heatmap
pressureplot = false;   % true = show pressure heatmap
onlychars    = true;    % true = only show characteristic lines


M_exit = 2.0;         %Starting M_exit values and slope calculations
patm_pexit_ratio = 0.5;
pb_p0_ratio = 1.0;
gamma = 1.4;
strang_exit = 0;

M_last = M_exit_to_M_last(M_exit, patm_pexit_ratio, gamma);
v_exit = prandtlmeyer_angled(M_exit);       %Prandtl Meyer angles in rad
v_last = prandtlmeyer_angled(M_last);
Mangle_exit = asind(1/M_exit);     %Mach angles
Mangle_last = asind(1/M_last);
symmetry_axis = 0;

stepsize_x = 0.01;      %Stepsize of x, smaller than 0.01 takes too long
radius = 0.5;           %Setting for starting points of nozzle
simlength = 10;        %Length of simulation (x axis)
s = simlength;                       
x = 0:stepsize_x:s;     %x values
Nx = length(x);         %Number of x values
N = 20;                  %Number of characteristic lines total
ymax = 1.5;             %Height of the mach distribution plot 1.2-1.7 good
gridpoints = 100;       %Number of gridpoints for the mach number mesh

y_axis_of_symmetry = symmetry_axis.* (1+zeros(1, Nx));  %Horizontal "wall"

%Stream angle at the start (Jet boundary angle in degrees)
strangle_last = moc_strang(v_exit, v_last, strang_exit, 'positive');

Mi = linspace(M_exit, M_last, N);
strangd = moc_strang_from_M(Mi, M_exit, strang_exit, 'pos');
strangd = strangd(:);
h = expfanslope_from_M(Mi, M_exit, strang_exit, 'positive');
h = h(:);

init_char = radius + h .* x;                            

[char_lines, yjetb, s, s2, s3] = solver_loop1clean(Mi, strangd, strang_exit, 'neg', x, init_char, y_axis_of_symmetry, gamma, v_last, radius);

[staticpressure, M_interp, X, Y] = machdistribution(gamma, x, char_lines, Mi, yjetb, s, s2, s3, simlength, v_last, N, ymax, gridpoints, M_exit);

if machplot
    figure
    pcolor(X, Y, M_interp)
    shading interp
    axis equal tight
    c = colorbar; ylabel(c, 'Mach number')
    hold on
    plot(x, char_lines, 'w', x, yjetb, 'r', x, y_axis_of_symmetry, 'b', 'LineWidth', 0.8);
    hold off
    xlabel('x'); ylabel('y')
    title('Jet stream plot – Mach number')
    xlim(xlimits)
    ylim(ylimits)
end

if pressureplot
    figure
    pcolor(X, Y, staticpressure)
    shading interp
    axis equal tight
    c = colorbar; ylabel(c, 'Static pressure ratio')
    hold on
    plot(x, char_lines, 'w', x, yjetb, 'r', x, y_axis_of_symmetry, 'b', 'LineWidth', 0.8);
    hold off
    xlabel('x'); ylabel('y')
    title('Jet stream plot – Static pressure ratio')
    xlim(xlimits)
    ylim(ylimits)
end

if onlychars && ~machplot && ~pressureplot
    figure
    hold on
    plot(x, char_lines, 'w', x, yjetb, 'r', x, y_axis_of_symmetry, 'b', 'LineWidth', 0.8);
    hold off
    axis equal tight
    xlabel('x'); ylabel('y')
    title('Jet stream plot – Characteristics only')
    xlim(xlimits)
    ylim(ylimits)
end