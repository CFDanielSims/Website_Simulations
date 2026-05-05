M1 = 2;
h1 = tan(asin(1/M1))

angle1 = asind(1/M1)

pa_p0 = 0.5*static_to_stag_ratio(M1);

M2 = M_from_statstag_ratio(pa_p0)

h2 = tan(asin(1/M2))

angle2 = asind(1/M2)

x = 0:0.1:20;

y1 = h1*x;
y2 = h2*x;

plot(x, y1, 'r--', x, y2, 'b--')
xlabel('x')
ylabel('y')
title('Mach wave expansion fans')