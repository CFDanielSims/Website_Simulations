function [b, f] = line_through_point(x0, y0, m)
% Returns y-intercept b for y = m*x + b so the line passes through (x0,y0).
% Also returns a function handle f(x) = m*x + b.
b = y0 - m*x0;
if nargout > 1
    f = @(x) m*x + b;
end
end
