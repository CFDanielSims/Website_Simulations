function T = viewarray(varargin)
narginchk(1, inf);
cols = {}; names = {}; maxlen = 0;

for k = 1:nargin
    arg = varargin{k};

    % resolve value + base name
    if ischar(arg) || isstring(arg)
        expr = char(arg);
        val  = evalin('base', expr);
        base = matlab.lang.makeValidName(expr);
    else
        val  = arg;
        base = inputname(k); if isempty(base), base = sprintf('Var%d',k); end
        base = matlab.lang.makeValidName(base);
    end

    % normalize shapes -> expand matrices into multiple variables (by column)
    if isvector(val)
        v = val(:);
        cols{end+1}  = v; %#ok<*AGROW>
        names{end+1} = base;
        maxlen = max(maxlen, numel(v));
    elseif ismatrix(val)
        [m,n] = size(val);
        for j = 1:n
            v = val(:,j);
            cols{end+1}  = v;
            names{end+1} = sprintf('%s_%d', base, j);
        end
        maxlen = max(maxlen, m);
    else
        error('Only vectors or 2-D matrices supported.');
    end
end

% pad with NaN to common length
for i = 1:numel(cols)
    if numel(cols{i}) < maxlen
        cols{i}(end+1:maxlen,1) = NaN;
    end
end

T = table(cols{:}, 'VariableNames', names);

uif = uifigure('Name','ViewData','Position',[100 100 900 500]);
uitable(uif, 'Data', T, 'Position', [0 0 900 500]);
drawnow;
end
