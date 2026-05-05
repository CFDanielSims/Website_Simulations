clear; clf; close all
M1 = 2;                             %Starting M1 values and slope calculations
pbp0ratio = 1;                          %Ratio between pressure of boundary and outside pressure
% pa_p0_ratio =; 

doplot = true;                      %Show plot or not (False)

%%External functions
pa_p0 = 0.5*static_to_stag_ratio(M1);
M2 = M_from_statstag_ratio(pa_p0);  %M2 Value based on M1

%%Prandtl Meyer angle
v1 = deg2rad(prandtlmeyer_angle(M1));
v2 = deg2rad(prandtlmeyer_angle(M2));

phi2 = v1-v2;                       %Streamline angle for last characteristic

%%Mach angle
Mangle1 = asin(1/M1);
Mangle2 = asin(1/M2);

h1 = tan(Mangle1);                  %Slope for equations
h2 = tan(Mangle2 - phi2);

rad2deg(Mangle2 + phi2);

yradius = 0.5;                      %Setting for starting points of nozzle
ylower = -yradius;
yupper = yradius;
simlength = 10;                     %Length of simulation
s = simlength;                       
x = 0:0.01:s;                       %x values
Nx = length(x);                     %Nx number of x values
N = 0;                              %Number of characteristics
h = linspace(h1, h2, N+2);          %Slopes of all characteristic lines
h = h(:);                           %Reshapes array into column vector

Ybottom = ylower + h .* x;          %N by Nx array of coordinates of characteristics
Ytop = -Ybottom;                    %Mirror, top side

%%reflection

reflectSlope = @(m, mb) tan(2*atan(mb) - atan(m)); %Gives slope of reflected line, m = incoming line, mb = slope jet boundary

%%Jet boundary

pm1 = prandtlmeyer_angle(M1);       %Prandtl-meyer angle of M1  
pm2 = prandtlmeyer_angle(M2);       %Prandtl-meyer angle of M2
hj1 = tand(pm2-pm1)                %Slope of upper jet boundary
hj1deg = atand(hj1)
jetb1 = yupper + hj1.*x;            %Upper jet boundary
jetb2 = -jetb1;                     %Lower jet boundary

% view_data_in_table(x, jetb1, Ybottom(1, :), Ybottom(2, :), h(1), h(2), hj1 )

%%Refleciton
[y1ref, x_int1, y_int1] = reflection1(x, jetb1, Ybottom(1, :), h(1), hj1);
[y2ref, x_int2, y_int2] = reflection1(x, jetb1, Ybottom(2, :), h(2), hj1);

if ~isequal(x_int1, false)
    % do something when an intersection was found
    [~, n1] = min(abs(x_int1-x));
    Ybottomcut1 = Ybottom(1, 1:n1);
    xbottomcut1 = x(1:n1);
    y1refcut = y1ref(n1:end);
    x1cut = x(n1:end);
else
    Ybottomcut1 = Ybottom(1, :); xbottomcut1 = x; y1refcut = []; x1cut = [];
end

if ~isequal(x_int2, false)
    [~, n2] = min(abs(x_int2-x));
    Ybottomcut2 = Ybottom(2, 1:n2);
    xbottomcut2 = x(1:n2);
    y2refcut = y2ref(n2:end);
    x2cut = x(n2:end);
else
    Ybottomcut2 = Ybottom(2, :); xbottomcut2 = x; y2refcut = []; x2cut = [];
end

% x_int2
%%First solver for wave characteristic interaction with jet boundary
% jet1b1cut = upperjetboundary2(x_int2, y_int2, x, phi2, pm2, pbp0ratio, n2);

%%PLOT
if doplot
    figure
    plot( xbottomcut1, Ybottomcut1.', 'm:', ...
          xbottomcut2, Ybottomcut2, 'c:',                          ...      % transpose 2×N -> N×2
          x, Ytop.',    'w:', ...
          x, jetb1,     'g-', ...
          x, jetb2,     'w-', ...
          x1cut, y1refcut,     'y-', ...
          x2cut, y2refcut,     'r-' );
    grid on
    xlabel('x'); ylabel('y');
    title('Mach wave expansion fans');
    legend('Ybottom','Ytop','jetb1','jetb2','y1ref','y2ref', 'Ybottomcut1', 'Ybottomcut2');
end

% angle1 = asind(1/M1)           %%for angles
% angle2 = asind(1/M2)


integercounter = 0:1:length(x);
% view_data_in_table(x, x_int1, integercounter)