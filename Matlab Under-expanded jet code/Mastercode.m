%%Mastercode
%%Now use symmetry so you only have to compute one expansion fan
doplot = true;                                         %Show plot or not (False)

M_exit = 2;                                             %Starting M_exit values and slope calculations
patm_pexit_ratio = 0.5;
pb_p0_ratio = 1;
gamma = 1.4;
strang_exit = 0;

M_last = M_exit_to_M_last(M_exit, patm_pexit_ratio, gamma);
v_exit = prandtlmeyer_angled(M_exit);                    %Prandtl Meyer angles in rad
v_last = prandtlmeyer_angled(M_last);
Mangle_exit = asind(1/M_exit);                           %Mach angles
Mangle_last = asind(1/M_last);
symmetry_axis = 0;

stepsize_x = 0.01;                                           %Stepsize of x
radius = 0.5;                                           %Setting for starting points of nozzle
simlength = 8;                                         %Length of simulation (x axis)
s = simlength;                       
x = 0:stepsize_x:s;                                           %x values
Nx = length(x);                                         %Number of x values
N = 20;                                                  %Number of characteristic lines total

y_axis_of_symmetry = symmetry_axis.* (1+zeros(1, Nx));  %Horizontal "wall"

strangle_last = moc_strang(v_exit, v_last, strang_exit, 'positive'); %Stream angle at the start (Jet boundary angle in degrees)
hjetb = tand(strangle_last);                            %Slope of upper jet boundary
hjetbdeg = atand(hjetb);
yjetb = radius + hjetb.*x;                              %Upper jet boundary

Mi = linspace(M_exit, M_last, N);
strangd = moc_strang_from_M(Mi, M_exit, strang_exit, 'pos');
strangd = strangd(:);
h = expfanslope_from_M(Mi, M_exit, strang_exit, 'positive');
h = h(:);

init_char = radius + h .* x;                            

[char_lines, yjetb] = solver_loop1(Mi, strangd, strang_exit, 'neg', x, init_char, y_axis_of_symmetry, gamma, v_last, radius);

if doplot
    plot(x, char_lines, 'w', x, yjetb, 'r', x, y_axis_of_symmetry, 'b');
end














% size(characterstic_lines)
% view_data_in_table_array(characteristic_lines)


% 
% h = linspace(h_exit, h_last, N+2);                      %Slopes of starting characteristic lines of expansion fan


% h = h(:)                                                %Reshapes array into column vector

% h_exit = tand(strangle_exit - Mangle_exit);              %Slope first characteristic (At exit)
% h_last = tand(strangle_last - Mangle_last);              %Slope last characteristic initial

% view_data_in_table(Mi, h, Strangd)

% characteristic_lines = zeros(Nx, N);                   %All y values of characteristic lines, maybe not necessary

% pm_strangle_char = zeros(Nx, Nh, Nv);

% char_lines1 = char_lines(:, :, 1);
% char_lines2 = char_lines(:, :, 2);
% char_lines3 = char_lines(:, :, 3);
% char_lines4 = char_lines(:, :, 4);
% char_lines5 = char_lines(:, :, 5);
% char_lines5 = char_lines(:, :, 6);

% view_data_in_table_array(char1)

% size(y_axis_of_symmetry);