function v = prandtlmeyer_angled(M)

gamma = 1.4;

gammaratio = sqrt((gamma+1)/(gamma-1));
machsquared = sqrt(M.^2-1);

v = gammaratio*atand(gammaratio.^-1*machsquared)-atand(machsquared);

