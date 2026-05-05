gamma = 1.4
M = 1:0.1:10

gammaratio = sqrt((gamma+1)/(gamma-1))
machsquared = sqrt(M.^2-1)

v = gammaratio*atand(gammaratio.^-1*machsquared)-atand(machsquared)

plot(M, v, 'r--')
xlabel('Mach number')
ylabel('Prandtl-Meyer angle');
title('Velocity vs. Prandtl Meyer angle');
grid on;