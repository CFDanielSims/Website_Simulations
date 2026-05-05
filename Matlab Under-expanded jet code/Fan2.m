% jet_boundary_simple.m
% Jet boundary of an underexpanded, inviscid, 2D supersonic jet
% First segment (ABC): boundary is the streamline from the nozzle lip
% that reaches ambient pressure in the initial simple expansion fan.
% Assumes: isentropic expansion, no shocks before pa is reached.
% Inputs you can change: Me, pe_over_pa, g, H, Xmax.

clear; clc;

% ---- Inputs ----
g           = 1.4;      % ratio of specific heats
Me          = 2.0;      % nozzle exit Mach number
pe_over_pa  = 2.0;      % underexpanded: pe = pe_over_pa * pa  (e.g., 2 -> pe = 2*pa)
H           = 1.0;      % lip height (plot scale)
Xmax        = 100.0*H;    % plot extent in x

% ---- Helper functions ----
nu = @(M) sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1)*(M.^2-1))) - atan(sqrt(M.^2-1));  % [rad]
mu = @(M) asin(1./M);                                                                % [rad]
p_over_p0 = @(M) (1 + 0.5*(g-1).*M.^2).^(-g/(g-1));                                  % isentropic

% Invert p/p0 -> M on supersonic branch by bisection
function M = M_from_pp0(pp0,g)
    f = @(M) (1 + 0.5*(g-1).*M.^2).^(-g/(g-1)) - pp0;
    a = 1+1e-12; b = 50; fa=f(a); fb=f(b);
    for it=1:100
        m = 0.5*(a+b); fm = f(m);
        if fa*fm <= 0, b=m; fb=fm; else, a=m; fa=fm; end
        if (b-a) < 1e-12, break; end
    end
    M = 0.5*(a+b);
end

% ---- States at exit and ambient (relative to the jet's p0) ----
pe_p0 = p_over_p0(Me);         % exit static-to-stagnation (jet total pressure)
pa_p0 = pe_p0 / pe_over_pa;    % because pe/pa = pe_over_pa -> pa/p0 = (pe/p0)/(pe/pa)

Ma    = M_from_pp0(pa_p0, g);  % ambient-pressure state reached by isentropic expansion

% ---- Turning angles ----
nu_e  = nu(Me);
nu_a  = nu(Ma);

% In a simple expansion from a straight lip, streamlines are straight.
% The streamline that reaches pa has constant flow angle:
phi_e = 0.0;                    % exit flow parallel to x-axis
phi_a = (nu_a - nu_e);          % from V+ = phi - nu = const  -> phi_a = nu_a - nu_e  (rad)

% ---- Jet boundary (first segment) as a straight streamline from the lip ----
x0 = 0; y0 = H;                 % nozzle lip
slope = tan(phi_a);
x = linspace(0, Xmax, 200);
y = y0 + slope*(x - x0);

% ---- Plot ----
figure(1); clf; hold on; axis equal; box on;
plot([0 0],[0 2*H],'k-','LineWidth',1);           % exit plane
plot([0 Xmax],[0 0],'k--');                       % centerline
plot(x,y,'r-','LineWidth',1.5);                   % jet boundary (first segment)
xlabel('x'); ylabel('y');
title('Underexpanded jet: initial jet boundary (simple expansion)');
legend({'Exit plane','Centerline','Jet boundary'},'Location','northeast');
grid on;

% ---- Print key numbers ----
fprintf('Me = %.4f, Ma(at p=pa) = %.4f\n', Me, Ma);
fprintf('nu_e = %.3f deg, nu_a = %.3f deg\n', nu_e*180/pi, nu_a*180/pi);
fprintf('Boundary streamline angle phi_a = %.3f deg\n', phi_a*180/pi);
fprintf('Mach angle at boundary state mu(Ma) = %.3f deg\n', asin(1/Ma)*180/pi);
