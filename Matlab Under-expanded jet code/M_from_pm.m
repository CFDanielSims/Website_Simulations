function M = M_from_pm(nu_deg, gamma)
% Returns Mach number M from Prandtl-Meyer angle (nu_deg) and gamma.
% nu_deg: Prandtl-Meyer angle in degrees
% gamma : Specific heat ratio (default 1.4)

if nargin < 2
    gamma = 1.4;
end

nu = deg2rad(nu_deg);

% Define Prandtl-Meyer function
f = @(M) sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M.^2-1))) - atan(sqrt(M.^2-1)) - nu;

% Solve numerically for M
M = fzero(f, [1, 50]); % valid range for supersonic flow
end
