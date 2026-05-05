function [staticpressure, M_interp, X, Y] = machdistribution(gamma, xreal, char_lines, Mi, yjetb, s, s2, s3, simlength, vjetb, N, ymax, gridpoints, M_exit)

%FOR S
%1st index is Prandtl Meyer angle, 
%2nd index is Stream Angle
%3rd index is x coordinate intersection
%4th index is y coordinate intersection
%5th index is gamma plus slope characteristic
%6th index is gamma minus slope characteristic
%7th index is index of x at which x_int coordinate is found closest to x

Nx = gridpoints; %Number of gridpoints
xmax = simlength;%Approximated from graph;
x = 0:xmax/Nx:xmax;
y = 0:ymax/Nx:ymax;
[X, Y] = meshgrid(x, y);

Nxreal = length(xreal); %Number of x points in xreal

ylist = zeros(N, s(N, 1, 7));
xlist = zeros(N, s(N, 1, 7));
Mlist = zeros(N, s(N, 1, 7));

ylist2 = zeros(N, s2(N, 1, 7)-s(N, 1, 7));
xlist2 = zeros(N, s2(N, 1, 7)-s(N, 1, 7));
Mlist2 = zeros(N, s2(N, 1, 7)-s(N, 1, 7));

ylist3 = zeros(N, s3(N, 1, 7)-s2(N, 1, 7));
xlist3 = zeros(N, s3(N, 1, 7)-s2(N, 1, 7));
Mlist3 = zeros(N, s3(N, 1, 7)-s2(N, 1, 7));

% y2list = zeros()

for j = 1:N
    
    l1 = char_lines(j, 1:s(j, 1, 7));
    m1 = xreal(1:s(j, 1, 7));
    n1 = Mi(j).*ones(1, s(j, 1, 7));

    l2 = char_lines(j, s(N, j, 7)+1:s2(j, 1, 7));
    m2 = xreal(s(N, j, 7)+1:s2(j, 1, 7));
    n2 = M_from_pm(s(N, j, 1)).*ones(1, s2(j, 1, 7)-s(N, j, 7));

    l3 = char_lines(j, s2(N, j, 7)+1:s3(j, 1, 7));
    m3 = xreal(s2(N, j, 7)+1:s3(j, 1, 7));
    n3 = M_from_pm(s2(N, j, 1)).*ones(1, s3(j, 1, 7)-s2(N, j, 7));
    % 
    % t3v(j) = s2(N, j, 1);
    % t3(j) = M_from_pm(s2(N, j, 1));

    ylist(j, 1:length(l1)) = l1;
    xlist(j, 1:length(m1)) = m1;
    Mlist(j, 1:length(n1)) = n1;

    ylist2(j, 1:length(l2)) = l2;
    xlist2(j, 1:length(m2)) = m2;
    Mlist2(j, 1:length(n2)) = n2;

    ylist3(j, 1:length(l3)) = l3;
    xlist3(j, 1:length(m3)) = m3;
    Mlist3(j, 1:length(n3)) = n3;
end

% viewarray(t3v, t3);
ylistflat = [ylist(:);ylist2(:);ylist3(:)];
xlistflat = [xlist(:);xlist2(:);xlist3(:)];
Mlistflat = [Mlist(:);Mlist2(:);Mlist3(:)];
% viewarray(ylistflat, xlistflat, Mlistflat);

% viewarray(ylistflat, xlistflat, Mlistflat);


%Reverse calculate the mach number based on Prandtl Meyer angle

%M_from_pm works only with scalars, so we have to use forloop

%For more accurate mach number distribution
gridexit = 10;
b = gridexit;

M_exit = [M_exit*ones(1, 2*b), M_exit, M_exit];
x_exit = [0.06.*ones(1, b), 0.4.*ones(1, b), 0, 0.16];
y_exit = [linspace(0.45, 0, b), linspace(0.27, 0, b), 0.499, 0.405];

M_exitflat = M_exit(:);
x_exitflat = x_exit(:);
y_exitflat = y_exit(:);

%Solving Mach number according to prandtl meyer number at each node
Mjetb = zeros(1, length(xreal));
[Ms, Ms2, Ms3] = deal(zeros(N, N));

for j = 1:Nxreal

        Mjetb(j) = M_from_pm(vjetb, gamma);
end

for j = 1:N
    
    for i = j:N

    Ms(i, j) = M_from_pm(s(i, j, 1), gamma);              
    Ms2(i, j) = M_from_pm(s2(i, j, 1), gamma);
    Ms3(i, j) = M_from_pm(s3(i, j, 1), gamma);

    end
end

Mjetbflat = Mjetb(:);
Msflat = Ms(:);
Ms2flat = Ms2(:);
Ms3flat = Ms3(:);

xs = s(:, :, 3);
xs2 = s2(:, :, 3);
xs3 = s3(:, :, 3);

xjetbflat = xreal(:);
xsflat = xs(:);
xs2flat = xs2(:);
xs3flat = xs3(:);

ys = s(:, :, 4);
ys2 = s2(:, :, 4);
ys3 = s3(:, :, 4);

yjetbflat = yjetb(:);
ysflat = ys(:);
ys2flat = ys2(:);
ys3flat = ys3(:);

Mcolumn = [Mlistflat; M_exitflat; Mjetbflat; Msflat; Ms2flat; Ms3flat];
x_column = [xlistflat; x_exitflat; xjetbflat; xsflat; xs2flat; xs3flat];
y_column = [ylistflat; y_exitflat; yjetbflat; ysflat; ys2flat; ys3flat];

mask = ~(Mcolumn==0 & x_column==0 & y_column==0);
Mcolumn = Mcolumn(mask);
x_column = x_column(mask);
y_column = y_column(mask);

% Interpolate across entire grid
F = scatteredInterpolant(x_column, y_column, Mcolumn, 'linear', 'nearest');
M_interp = F(X, Y);

% viewarray(Mcolumn, x_column, y_column);

% Known data (scattered points)
% x_known = [1 3 5];
% y_known = [0.1 0.5 0.6];
% T_known = [30 50 70];

% Pcolumn = static_to_stag_ratio([Mlistflat; M_exitflat; Mjetbflat; Msflat; Ms2flat; Ms3flat], gamma)
% 
% Pcolumn(length(Mlistflat + M_exitflat):length(Mlistflat + M_exitflat+ Mjetbflat)) = 1

staticpressure = static_to_stag_ratio(M_interp, gamma);


