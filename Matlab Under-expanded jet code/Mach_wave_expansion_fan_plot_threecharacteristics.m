M1 = 2;  %Starting M1 values and slope calculations
h1 = tan(asin(1/M1))

pa_p0 = 0.5*static_to_stag_ratio(M1);

M2 = M_from_statstag_ratio(pa_p0); 

h2 = tan(asin(1/M2)); 

% angle1 = asind(1/M1); %%for angles
% angle2 = asind(1/M2); 

yradius = 0.5; %Setting for lower starting point of nozzle

ylower = -yradius;
yupper = yradius;

x = 0:0.1:5;

hbetween = linspace(h1, h2, 5);
h11 = hbetween(2)
h12 = hbetween(3)
h13 = hbetween(4)

y1 = ylower + h1*x;
y2 = ylower + h2*x;
y11 = ylower + h11*x;
y12 = ylower + h12*x;
y13 = ylower + h13*x;

ny1 = -y1
ny2 = -y2
ny11 = -y11
ny12 = -y12
ny13 = -y13

plot(x, y1, 'r--', x, y2, 'b--', x, y11, 'w', x, y12, 'w', x, y13, 'w', ...
    x, ny1, 'r--', x, ny2, 'b--', x, ny11, 'w', x, ny12, 'w', x, ny13, 'w')
xlabel('x')
ylabel('y')
title('Mach wave expansion fans')



