%% Intersection point finder (returns nearest x-index too)
function [x_int, y_int, idx_near] = intersection_finder_index(x, xrange, y1, y2)

doplot = false;

% Checks
if nargin < 3, error('Not enough arguments'); end
if ~isvector(xrange) || ~isequal(size(xrange), size(y1), size(y2))
    error('x, y1, y2 must be equal-size vectors');
end
if any(~isfinite(xrange)) || any(~isfinite(y1)) || any(~isfinite(y2))
    error('Inputs must be finite');
end
if ~issorted(xrange) || any(diff(xrange)==0)
    error('x must be strictly increasing');
end

dy  = y1 - y2;
idx = find(dy(1:end-1).*dy(2:end) <= 0);

x_int   = zeros(size(idx));
y_int   = zeros(size(idx));
idx_near = zeros(size(idx));             % << add: preallocate nearest-index

for k = 1:numel(idx)
    x1 = xrange(idx(k)); x2 = xrange(idx(k)+1);
    f1 = @(xq) interp1(xrange, y1, xq, 'linear');
    f2 = @(xq) interp1(xrange, y2, xq, 'linear');
    f  = @(xq) f1(xq) - f2(xq);
    x_int(k) = fzero(f, [x1 x2]);
    y_int(k) = f1(x_int(k));
    [~, idx_near(k)] = min(abs(x - x_int(k)));   % map to full x

end

if isempty(idx)
    x_int = false;
    y_int = false;
    idx_near = [];                          % << add: no intersections
    return
end

if doplot
    plot(xrange,y1,'r',xrange,y2,'b--'); hold on; grid on;
    scatter(x_int,y_int,50,'k','filled');
    legend('y1','y2','intersection');
end
