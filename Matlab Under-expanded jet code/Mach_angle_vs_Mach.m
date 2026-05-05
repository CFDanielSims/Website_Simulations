M = 0.1:0.01:10

mu = asind(1./M)

plot(M, mu, 'r--')
xlabel('Mach speed')
ylabel('Mach angle')
title('Mach angle vs. Mach speed')