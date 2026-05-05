x = -5:0.01:5;

y1 = 1 - x.^2;

y2 = -2 + x.^2;

[x_int, y_int] = intersection_finder(x, y1, y2)