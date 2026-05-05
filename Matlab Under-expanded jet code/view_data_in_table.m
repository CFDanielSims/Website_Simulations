function T = view_data_in_table(varargin)
% Opens selected data side-by-side in a UI table (no Command Window spam).
% Usage:
%   view_data_in_table(x, jetb1, Ybottom(1,:))
%   view_data_in_table('x','jetb1','Ybottom(1,:)')

narginchk(1, inf);

cols  = cell(1, nargin);
names = cell(1, nargin);
maxlen = 0;

for k = 1:nargin
    arg = varargin{k};
    if ischar(arg) || isstring(arg)
        expr = char(arg);
        data = evalin('base', expr);
        nm   = matlab.lang.makeValidName(expr);
    else
        data = arg;
        nm   = inputname(k); if isempty(nm), nm = sprintf('Var%d', k); end
        nm   = matlab.lang.makeValidName(nm);
    end
    data = data(:);
    cols{k}  = data;
    names{k} = nm;
    maxlen   = max(maxlen, numel(data));
end

for k = 1:nargin
    if numel(cols{k}) < maxlen
        cols{k}(end+1:maxlen,1) = NaN;
    end
end

T = table(cols{:}, 'VariableNames', names);

% Always open a dedicated spreadsheet-like window (works in desktop & Online)
uif = uifigure('Name','ViewData','Position',[100 100 900 500]);
uit = uitable(uif, 'Data', T, 'Position', [0 0 900 500]);
drawnow;  %#ok<*NASGU>
end
