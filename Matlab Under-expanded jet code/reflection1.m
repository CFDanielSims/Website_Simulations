function [f, x_int, y_int] = reflection1(x, y1, y2, sy1, sy2)

reflectSlope = @(m, mb) tan(2*atan(mb) - atan(m)); %Gives slope of reflected line, m = incoming line, mb = slope jet boundary

[x_int, y_int] = intersection_finder(x, y1, y2);
    
First_outgoing_slope = reflectSlope(sy1, sy2);

[bref, fref] = line_through_point(x_int, y_int, First_outgoing_slope); 

f = fref(x);