% underexpanded_jet_step1.m
% Step 1: compute ambient state and plot the first simple expansion fan (ABC)
% Gamma=1.4, Me=2, pe=2*pa  -> pa = pe/2
clear; clc;

% ---------------- Inputs ----------------
g   = 1.4;         % ratio of specific heats
Me  = 2.0;         % exit Mach number
N    = 25;         % number of characteristics in the initial fan (>=10 later)
H    = 1.0;        % half-height scale (sets y-coordinate of the lip at +H)
Lseg = 1.5*H;      % length to draw each characteristic segment

% --------------- Helpers ----------------
% Prandtl-Meyer nu(M) [rad]
nu = @(M) sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(M.^2-1))) - atan(sqrt(M.^2-1));
% Inverse: M from nu (scalar or array) via Newton-bisection
function M = M_from_nu(nu_t,g)
    M = zeros(size(nu_t));
    for k=1:numel(nu_t)
        f  = @(M) sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(M.^2-1))) - atan(sqrt(M.^2-1)) - nu_t(k);
        a=1+1e-8; b=50; fa=f(a); fb=f(b);
        for it=1:80
            m = 0.5*(a+b); fm=f(m);
            if fa*fm<=0, b=m; fb=fm; else, a=m; fa=fm; end
            if abs(b-a)<1e-10, break; end
        end
        M(k) = 0.5*(a+b);
    end
end
% Isentropic p/p0
p_over_p0 = @(M) (1 + 0.5*(g-1).*M.^2).^(-g/(g-1));
% Invert isentropic relation for M given p/p0 (supersonic branch)
function M = M_from_p_over_p0(pp0,g)
    f = @(M) (1 + 0.5*(g-1).*M.^2).^(-g/(g-1)) - pp0;
    a=1+1e-8; b=50; fa=f(a); fb=f(b);
    for it=1:80
        m = 0.5*(a+b); fm=f(m);
        if fa*fm<=0, b=m; fb=fm; else, a=m; fa=fm; end
        if abs(b-a)<1e-12, break; end
    end
    M = 0.5*(a+b);
end
% Mach angle mu(M) [rad]
mu = @(M) asin(1./M);

% ------------- Exit & ambient -----------
pe_p0 = p_over_p0(Me);      % exit static-to-stagnation
pa_p0 = 0.5*pe_p0;          % from pe = 2*pa  -> pa/pe = 1/2
Ma    = M_from_p_over_p0(pa_p0,g); % ambient Mach (supersonic)
nu_e  = nu(Me);
nu_a  = nu(Ma);
phi_e = 0.0;                % parallel exit flow
phi_a = nu_a - nu_e;        % from V+ invariant (phi - nu) const across ABC

% ------------- Build initial fan ABC -------------
% Discretize flow turning from phi_e -> phi_a (equally spaced in phi)
phi_i = linspace(phi_e, phi_a, N);                % [rad]
nu_i  = (phi_e - nu_e) + phi_i;                   % V+ const: nu_i = phi_e - nu_e + phi_i
M_i   = M_from_nu(nu_i, g);
mu_i  = mu(M_i);

% Geometry: treat Γ− lines as straight in a simple wave when V+ is const.
% Slope of Γ−: dy/dx = tan(phi - mu). All emanate from the nozzle lip (0,H).
x0 = 0.0; y0 = H;
m_minus = tan(phi_i - mu_i);

% Pick an x-extent and draw straight segments
x1 = Lseg*ones(1,N);
y1 = y0 + m_minus.*(x1 - x0);

% ------------- Plot ----------------
figure(1); clf; hold on; box on; axis equal;
% nozzle lip and centerline
plot([x0 x0],[0 2*H],'k-','LineWidth',1.0);                 % exit plane (for reference)
y_center = 0; plot([0, max(x1)], [y_center, y_center], 'k--'); % symmetry line

% plot characteristics of the initial expansion fan
for i=1:N
    plot([x0 x1(i)], [y0 y1(i)], 'b-');
end

% annotate
xlabel('x'); ylabel('y'); title(sprintf('Initial expansion fan (ABC), N=%d',N));
legend({'Exit plane','Centerline','\Gamma^- in ABC'},'Location','northeast');
grid on;

% quick text dump of key numbers
fprintf('Exit:  M_e=%.4f, nu_e=%.4f deg\n', Me, nu_e*180/pi);
fprintf('Amb.:  M_a=%.4f, nu_a=%.4f deg, phi_a=%.4f deg\n', Ma, nu_a*180/pi, phi_a*180/pi);
