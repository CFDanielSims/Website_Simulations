clear all
close all
clc

counter = 0;

for i = 1:5
    counter = counter + 1
    disp(counter)
end

while counter >= 5
    counter = counter-1;
    disp(counter)
end

% Plotting

x = 0:0.1:5;
y = x.^2;

plot(x,y,'r+')
title('My first plot')

xlabel('x')
ylabel('y')

grid on

hold

y2 = x.^3
y3 = x.^4

plot(x, y2, 'g*')
plot(x, y3)
hold off
legend('Plot 1','Plot 2', 'Plot 3')