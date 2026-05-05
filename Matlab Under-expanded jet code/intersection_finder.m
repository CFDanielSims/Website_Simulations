%%Intersection point finder
function [x_int, y_int] = intersection_finder(x, y1, y2);

doplot = false;                %With or without plot

%Checks
if nargin < 3, error('Not enough arguments'); end
if ~isvector(x) || ~isequal(size(x), size(y1), size(y2))
    error('x, y1, y2 must be equal-size vectors');
end
if any(~isfinite(x)) || any(~isfinite(y1)) || any(~isfinite(y2))
    error('Inputs must be finite');
end
if ~issorted(x) || any(diff(x)==0)
    error('x must be strictly increasing');
end

dy = y1 - y2;

idx = find(dy(1:end-1).*dy(2:end)<=0);

x_int = zeros(size(idx));
y_int = zeros(size(idx));

for k = 1:length(idx)
    x1 = x(idx(k)); x2 = x(idx(k)+1);
    f1 = @(xq) interp1(x, y1, xq, 'linear');                
    f2 = @(xq) interp1(x, y2, xq, 'linear');                %Linearily interpolates for a function y with as variable x at xq
    f  = @(xq) f1(xq) - f2(xq);                             %Difference y between intersection points at both equations
    x_int(k) = fzero(f, [x1 x2]);                           %Finds when xq is equal to 0 (When f1 = f2, linear approximation)
    y_int(k) = f1(x_int(k));                                %Y value for xq
end

if isempty(idx)
    x_int = false;
    y_int = false;
    return
end

%fprintf('Intersection near (%.4f, %.4f)\n', x_int, y_int)  %Uncomment for
%more detail

% Plot check
if doplot;
    plot(x,y1,'r',x,y2,'b--'); hold on; grid on;
    scatter(x_int,y_int,50,'k','filled');
    legend('y1','y2','intersection');
end

