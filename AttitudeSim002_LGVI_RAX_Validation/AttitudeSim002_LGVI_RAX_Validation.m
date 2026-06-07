%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Replication of Park et al.'s SO(3) Lie Group Variational Integrator (LGVI) for RAX (Scenario 4)
%% SO(3) LGVI with IGRF-11 and Continuous Flatley-Henretty Hysteresis
%% Reference Paper: "A Dynamic Model of a Passive Magnetic Attitude Control System for
%% the RAX Nanosatellite"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAX PARAMETERS SCENARIO 4 USED IN THIS SIMULATION
% Orbit type:        circular
% altitude:          650 km
% timespan:          6 orbital periods
% field model:       IGRF
% config:            3U
% inertia:           J = diag(2.91058, 2.91058, 0.59261) * 10^(-2) kg * m^2
% magnet B_p:        1.28 T
% magnet axis:       z axis
% magnet V_p:        3.0137 * 10^(-6) m^3
% rod V_h:           7.15 * 10^(-8) m^3
% B_m:               0.73 T
% H_c:               1.59 A/m
% H_r:               1.696 A/m
% Integrator:        SO(3) Lie Group Variational Integrator (2nd Order, fixed timestep)
% Step size:         0.05 s fixed timestep
% Total sim time:    6 orbital periods
% Hysteresis update: Forward Euler
% Initial velocity:  [0.05, 0.05, 0.05] rad/s
% Initial rotation:  [0.9947 -0.1008 -0.0221
%                     0.0102  0.1172 -0.9931
%                     0.1027  0.9880  0.1155]

%% This script is designed to reproduce the paper's Scenario 4 assumptions.
% - Circular orbit: altitude 650 km, inclination 72 deg
% - Total time: 6 orbital periods
% - The IGRF-11 precompute may take several minutes for dtB = 0.05 s.
% - The paper specifies that an initial orbital position is used that
% provides zero initial pointing error. A searching algorithm is used to
% find an initial orbital geometry that satisfies the condition of zero
% pointing error.
% - Earth field: IGRF-11 (ECEF), transformed using the paper's R_I^E(t) rotation
% - Integrator: Lie Group Variational Integrator (2nd-order, fixed timestep)
% - Hysteresis: Continuous Flatley-Henretty ODE, forward Euler update

%%%%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
warning('off', 'all');

% Optional global matrix for logging the number of iterations needed for
% Newton solver convergence.
global NEWTON_ITERS 
NEWTON_ITERS = [];

%% -------------------- USER SETTINGS -------------------------------------
h = 0.05;                  % [s] fixed timestep, paper value
dtB = h;                   % [s] coarse IGRF precompute spacing
igrf_decimal_year = 2010;  % 2010-era RAX paper
igrf_generation  = 11;     % use if supported by MATLAB; fallback included
print_every_seconds_wall = 10;

% These will be auto-updated by the search below
Omega0 = 0;                % [rad] initial RAAN
u0     = 0;                % [rad] initial argument of latitude

%% -------------------- CONSTANTS -----------------------------------------
mu_earth = 3.986004418e14;   % [m^3/s^2]
Re_mean  = 6371e3;           % [m]
mu0      = 4*pi*1e-7;        % [N/A^2]
omega0   = 7.292115e-5;      % [rad/s] Earth rotation rate

%% -------------------- RAX SCENARIO 4 PARAMETERS -------------------------
altitude = 650e3;            % [m]
incl     = deg2rad(72.0);    % [rad]

r_orbit  = Re_mean + altitude;
n_orbit  = sqrt(mu_earth / r_orbit^3);
T_orbit  = 2*pi / n_orbit;
t_end    = 6 * T_orbit;
steps    = floor(t_end / h);
t        = (0:steps)' * h;
N        = numel(t);

fprintf('RAX Scenario 4: h = %.3f s, N = %d, sim time = %.2f hr (%.2f orbits)\n', ...
    h, N, t(end)/3600, t(end)/T_orbit);

% Inertia matrix
J    = diag([2.91058, 2.91058, 0.59261]) * 1e-2;   % [kg*m^2]
invJ = diag(1 ./ diag(J));

% Permanent magnet
Bp = 1.28;              % [T]
Vp = 3.0137e-6;         % [m^3]

% Hysteresis rods
Vh = 7.15e-8;           % [m^3]
Bm = 0.73;              % [T]
Hc = 1.59;              % [A/m]
Hr = 1.696;             % [A/m]

% Body principal axes
b1 = [1;0;0];
b2 = [0;1;0];
b3 = [0;0;1];

% Initial orientation
R_BI = [ 0.9947 -0.1008 -0.0221;
         0.0102  0.1172 -0.9931;
         0.1027  0.9880  0.1155 ];

% Initial angular velocity
omega = [0.05; 0.05; 0.05];    % [rad/s], Scenario 4

% Initial angular momentum
Pi    = J * omega;

% Initial magnetic flux density within hysteresis rods
Brod = [0; 0];                 % [B1; B2]

%%%%%%%%%%% AUTO-FIND ORBIT GEOMETRY FOR ZERO INITIAL POINTING ERROR %%%%%%
% Park et al. specifies that their simulation is configured to an initial
% state in which the pointing error is zero. The authors do not specify the
% initial location of the spacecraft in orbit nor the IGRF model generation
% or deployment date. In this replication, IGRF-11 is evaluated at the 2010.0
% epoch. A searching algorithm is used to find an initial orbital geometry 
% to satisfy the condition of zero pointing error.
fprintf('Searching for orbit geometry (Omega0, u0) that minimizes initial pointing error...\n');

% ---------- COARSE SEARCH ------------------------------------------------
Omega_candidates = linspace(0, 2*pi, 181);
Omega_candidates(end) = [];

u0_candidates = linspace(0, 2*pi, 181);
u0_candidates(end) = [];

best_cost   = inf;
best_theta0 = inf;
best_H10    = inf;
best_H20    = inf;

REI0_search = R_I2E(omega0, 0);

for iO = 1:numel(Omega_candidates)
    Omega_test = Omega_candidates(iO);

    for iu = 1:numel(u0_candidates)
        u0_test = u0_candidates(iu);

        r_eci0  = circularOrbitECI(r_orbit, incl, Omega_test, u0_test);
        r_ecef0 = REI0_search * r_eci0;

        B_ecef0 = igrfECEF_Tesla(r_ecef0, igrf_decimal_year, igrf_generation);
        H_E0    = B_ecef0 / mu0;
        H_body0 = R_BI * (REI0_search.') * H_E0;

        H10 = b1.' * H_body0;
        H20 = b2.' * H_body0;

        theta0_test = acosd(max(-1, min(1, (b3.' * H_body0) / norm(H_body0))));
        cost = theta0_test^2;

        if cost < best_cost
            best_cost   = cost;
            best_theta0 = theta0_test;
            best_H10    = H10;
            best_H20    = H20;
            Omega0      = Omega_test;
            u0          = u0_test;
        end
    end
end

fprintf('Coarse best Omega0 = %.6f rad\n', Omega0);
fprintf('Coarse best u0     = %.6f rad\n', u0);
fprintf('Coarse best theta0 = %.6e deg\n', best_theta0);

% ---------- FINE SEARCH --------------------------------------------------
deg = pi/180;
Omega_center = Omega0;
u0_center    = u0;

Omega_candidates_fine = linspace(Omega_center - 4*deg, Omega_center + 4*deg, 81); % 0.1 deg
u0_candidates_fine    = linspace(u0_center    - 4*deg, u0_center    + 4*deg, 81); % 0.1 deg

best_cost_fine   = inf;
best_theta0_fine = inf;
best_H10_fine    = inf;
best_H20_fine    = inf;

for iO = 1:numel(Omega_candidates_fine)
    Omega_test = mod(Omega_candidates_fine(iO), 2*pi);

    for iu = 1:numel(u0_candidates_fine)
        u0_test = mod(u0_candidates_fine(iu), 2*pi);

        r_eci0  = circularOrbitECI(r_orbit, incl, Omega_test, u0_test);
        r_ecef0 = REI0_search * r_eci0;

        B_ecef0 = igrfECEF_Tesla(r_ecef0, igrf_decimal_year, igrf_generation);
        H_E0    = B_ecef0 / mu0;
        H_body0 = R_BI * (REI0_search.') * H_E0;

        H10 = b1.' * H_body0;
        H20 = b2.' * H_body0;

        theta0_test = acosd(max(-1, min(1, (b3.' * H_body0) / norm(H_body0))));
        cost = theta0_test^2;

        if cost < best_cost_fine
            best_cost_fine   = cost;
            best_theta0_fine = theta0_test;
            best_H10_fine    = H10;
            best_H20_fine    = H20;
            Omega0           = Omega_test;
            u0               = u0_test;
        end
    end
end

fprintf('Fine best Omega0 = %.6f rad\n', Omega0);
fprintf('Fine best u0     = %.6f rad\n', u0);
fprintf('Fine best theta0 = %.6e deg\n', best_theta0_fine);
fprintf('Fine best H1     = %.6e A/m\n', best_H10_fine);
fprintf('Fine best H2     = %.6e A/m\n', best_H20_fine);

%% -------------- PRECOMPUTE ORBIT + IGRF IN EARTH FRAME ------------------
t_coarse = (0:dtB:t(end))';
if abs(t_coarse(end) - t(end)) > eps(t(end))
    t_coarse = [t_coarse; t(end)];
end

cache_filename = sprintf( ...
    'RAX_S4_IGRFcache_dtB_%s_year_%g_gen_%d_Omega_%s_u0_%s.mat', ...
    strrep(num2str(dtB, '%.6g'), '.', 'p'), ...
    igrf_decimal_year, ...
    igrf_generation, ...
    strrep(num2str(Omega0, '%.6g'), '.', 'p'), ...
    strrep(num2str(u0, '%.6g'), '.', 'p'));

if exist(cache_filename, 'file')
    fprintf('Loading IGRF cache: %s\n', cache_filename);
    S = load(cache_filename, 't_coarse', 'H_e_coarse', 'Hdot_e_coarse', 'H_e_hist', 'Hdot_e_hist');
    H_e_hist    = S.H_e_hist;
    Hdot_e_hist = S.Hdot_e_hist;

    if size(H_e_hist,1) ~= N || size(Hdot_e_hist,1) ~= N
        error('Cached IGRF history size does not match current simulation grid. Delete cache and rerun.');
    end
else
    fprintf('Precomputing IGRF values along orbit (coarse grid, dtB = %.2f s)...\n', dtB);
    tic;

    NC = numel(t_coarse);
    H_e_coarse = zeros(NC, 3);

    for k = 1:NC
        tk = t_coarse(k);
        u  = u0 + n_orbit * tk;

        % Position along circular orbit in inertial frame
        r_eci = circularOrbitECI(r_orbit, incl, Omega0, u);

        % Inertial -> Earth-fixed
        REI = R_I2E(omega0, tk);
        r_ecef = REI * r_eci;

        % IGRF field in ECEF, then convert B -> H
        B_ecef = igrfECEF_Tesla(r_ecef, igrf_decimal_year, igrf_generation);
        H_e_coarse(k, :) = (B_ecef / mu0).';
    end

    % Time derivative of H_E on the coarse grid
    Hdot_e_coarse = zeros(NC, 3);
    if NC >= 2
        Hdot_e_coarse(1,  :) = (H_e_coarse(2,  :) - H_e_coarse(1,    :)) / (t_coarse(2)   - t_coarse(1));
        Hdot_e_coarse(end,:) = (H_e_coarse(end,:) - H_e_coarse(end-1,:)) / (t_coarse(end) - t_coarse(end-1));
    end
    for k = 2:NC-1
        Hdot_e_coarse(k, :) = (H_e_coarse(k+1, :) - H_e_coarse(k-1, :)) / (t_coarse(k+1) - t_coarse(k-1));
    end
    
    % Uncomment these lines if dtB ~= h.
    % Interpolate coarse data onto full simulation grid
    %H_e_hist    = interp1(t_coarse, H_e_coarse,    t, 'linear', 'extrap');
    %Hdot_e_hist = interp1(t_coarse, Hdot_e_coarse, t, 'linear', 'extrap');
    
    H_e_hist    = H_e_coarse;
    Hdot_e_hist = Hdot_e_coarse;

    fprintf('IGRF precompute finished in %.1f s. Saving cache...\n', toc);
    save(cache_filename, 't_coarse', 'H_e_coarse', 'Hdot_e_coarse', 'H_e_hist', 'Hdot_e_hist', '-v7.3');
end

%% -------------------- INITIAL GEOMETRY DIAGNOSTIC -----------------------
REI0 = R_I2E(omega0, t(1));
H_body0 = R_BI * (REI0.') * H_e_hist(1, :).';
theta0 = acosd(max(-1, min(1, (b3.' * H_body0) / norm(H_body0))));
H10 = b1.' * H_body0;
H20 = b2.' * H_body0;

fprintf('Initial pointing error after cache/load = %.6e deg\n', theta0);
fprintf('Initial H1 after cache/load = %.6e A/m\n', H10);
fprintf('Initial H2 after cache/load = %.6e A/m\n', H20);

%% -------------------- DATA STORAGE --------------------------------------
theta_deg  = zeros(N, 1);
omega_hist = zeros(N, 3);
Pi_hist    = zeros(N, 3);
Brod_hist  = zeros(N, 2);
Hrod_hist  = zeros(N, 2);

ineq_resid_upper = zeros(N, 2);
ineq_resid_lower = zeros(N, 2);

sat_guard_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SIMULATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Starting LGVI integration...\n');
tic;
last_print = 0;

for k = 1:N
    tk     = t(k); % current simulation time
    REI    = R_I2E(omega0, tk); % inertial-to-Earth-fixed rotation matrix
    REIdot = R_I2E_dot(omega0, tk); % time derivative of REI

    H_E    = H_e_hist(k, :).'; % magnetic field vector in Earth-fixed frame
    Hdot_E = Hdot_e_hist(k, :).'; % time derivative of H_E

    H_body = R_BI * (REI.') * H_E; % provides field in body-frame

    H1 = b1.' * H_body; % projected field onto hysteresis rod axis
    H2 = b2.' * H_body;

    % Pointing error relative to b3
    Hmag = norm(H_body); % magnitude of magnetic field
    if Hmag > 0
        theta_deg(k) = acosd(max(-1, min(1, (b3.' * H_body) / Hmag))); % pointing error
    else
        theta_deg(k) = 0;
    end

    omega = invJ * Pi; % angular velocity recovery
    
    % Compute field-strength derivatives used by the hysteresis model.
    % These terms account for spacecraft's rotation through the field,
    % orbital translation through the spatially varying field, and
    % Earth-frame rotation of the magnetic field.
    common_mat = -skew(omega) * R_BI * (REI.') + R_BI * (REIdot.');
    Hdot1 = b1.' * (common_mat * H_E + R_BI * (REI.') * Hdot_E); % time derivative of magnetic field along rod axis
    Hdot2 = b2.' * (common_mat * H_E + R_BI * (REI.') * Hdot_E);

    % Store
    omega_hist(k, :) = omega.';
    Pi_hist(k, :)    = Pi.';
    Brod_hist(k, :)  = Brod.';
    Hrod_hist(k, :)  = [H1, H2];

    % Hysteresis invariant region check
    xi1 = H1 - Hr * tan(pi * Brod(1) / (2 * Bm));
    xi2 = H2 - Hr * tan(pi * Brod(2) / (2 * Bm));
    % Store upper/lower inequality residuals
    ineq_resid_upper(k, :) = [xi1 - Hc, xi2 - Hc];
    ineq_resid_lower(k, :) = [-Hc - xi1, -Hc - xi2];

    if k == N
        break;
    end

    % Hysteresis update: forward Euler, with a singularity guard
    [B1_dot, sat1] = duhemHenretty_Bdot_guarded(H1, Hdot1, Brod(1), Bm, Hc, Hr);
    [B2_dot, sat2] = duhemHenretty_Bdot_guarded(H2, Hdot2, Brod(2), Bm, Hc, Hr);
    sat_guard_count = sat_guard_count + sat1 + sat2;
    Brod_next = Brod + h * [B1_dot; B2_dot];

    % Momentum is calculated via a second-order splitting scheme (kick-
    % drift-kick). 
    % Current torques calculation
    M_k = (Bp * Vp) * cross(b3, H_body) + ...
          (Brod(1) * Vh) * cross(b1, H_body) + ...
          (Brod(2) * Vh) * cross(b2, H_body);

    % Newton solver used to find rotation increment F
    Pi_half = Pi + 0.5 * h * M_k; % first half-kick
    F = lgviSolveF(J, h, Pi_half);

    % For R_BI satisfying Rdot + S(w)R = 0
    R_BI_next = F.' * R_BI; % drift

    % Next-step torque calculation
    tkp1     = t(k+1);
    REI_next = R_I2E(omega0, tkp1);
    H_E_next = H_e_hist(k+1, :).';
    H_body_next = R_BI_next * (REI_next.') * H_E_next;

    M_kp1 = (Bp * Vp) * cross(b3, H_body_next) + ...
            (Brod_next(1) * Vh) * cross(b1, H_body_next) + ...
            (Brod_next(2) * Vh) * cross(b2, H_body_next);

    % Momentum update
    Pi_next = F.' * Pi_half + 0.5 * h * M_kp1; % second half-kick

    % Commit
    R_BI = R_BI_next;
    Pi   = Pi_next;
    Brod = Brod_next;

    wall = toc;
    if wall - last_print > print_every_seconds_wall
        fprintf(' %5.1f%% | sim t = %.2f hr | wall = %.0f s\n', ...
            100 * k / N, tk / 3600, wall);
        last_print = wall;
    end
end
fprintf('Done. Wall time: %.1f s\n', toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write logging matrix to text file. This is the number of Newton solver iterations
% needed for residual convergence for every Newton solver call. This is
% optional.
global NEWTON_ITERS
% writematrix(NEWTON_ITERS','newton_iterations.txt');

%% -------------------- REPORT HYSTERESIS REGION CHECK --------------------
max_upper_violation = max(max(ineq_resid_upper));
max_lower_violation = max(max(ineq_resid_lower));

fprintf('\nHysteresis region check (paper Eq. 11):\n');
fprintf(' max upper violation = %.3e\n', max_upper_violation);
fprintf(' max lower violation = %.3e\n', max_lower_violation);
if max_upper_violation > 0 || max_lower_violation > 0
    fprintf(' WARNING: hysteresis state left the invariant region at some step.\n');
else
    fprintf(' All simulated states remained inside the invariant region.\n');
end
fprintf('Saturation guard activations = %d\n', sat_guard_count);

%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orbits = t / T_orbit;

results_filename = fullfile(pwd, 'AttitudeSim002_Results.mat');
save(results_filename, ...
    't', 'orbits', 'theta_deg', 'omega_hist', 'Pi_hist', ...
    'Brod_hist', 'Hrod_hist', ...
    'ineq_resid_upper', 'ineq_resid_lower', ...
    'NEWTON_ITERS');
fprintf('Saved results to: %s\n', results_filename);

% Global formatting: Times New Roman
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultLegendFontName','Times New Roman');

% Consistent Sizing for Side-by-Side in Word
fontAxis  = 12; 
fontLabel = 14; 
figWidth  = 4.5; % Inches
figHeight = 3.5; % Inches

% Mission Dashboard Color Scheme (from your dashboard snippet)
c_wx   = [0 0.25 0.85];   % blue
c_wy   = [0.85 0.1 0.1];  % red
c_wz   = [0.93 0.69 0.13];% yellow
c_beta = [0 0 0];         % black

%% (a) Pointing Error
fig1 = figure('Color','w', 'Units','inches', ...
              'Position',[1 1 figWidth figHeight], ...
              'Name','Pointing Error');

% Solid line, reduced thickness to keep the dense data clear
plot(orbits, theta_deg, '-', 'Color', c_beta, 'LineWidth', 0.8)
grid on
box on

xlabel('Orbit Number', 'FontSize', fontLabel)
ylabel('Pointing Error [deg]', 'FontSize', fontLabel)
xticks(0:1:6)
xlim([0 6])
yticks(0:15:90)
ylim([0 90])

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

% Export as high-res PNG for easy "Insert" or "Paste" into Word
% exportgraphics(fig1, 'RAX_Scenario4_pointing_error.png', 'Resolution', 600)

%% (b) Angular Velocity
fig2 = figure('Color','w', 'Units','inches', ...
              'Position',[6 1 figWidth figHeight], ... 
              'Name','Angular Velocity');

% Solid lines only, using the RGB dashboard colors
plot(orbits, omega_hist(:,1), '-', 'Color', c_wx, 'LineWidth', 1.2); hold on
plot(orbits, omega_hist(:,2), '-', 'Color', c_wy, 'LineWidth', 1.2); 
plot(orbits, omega_hist(:,3), '-', 'Color', c_wz, 'LineWidth', 1.2);  

grid on
box on

xlabel('Orbit Number', 'FontSize', fontLabel)
ylabel('Angular Velocity [rad/s]', 'FontSize', fontLabel)
xticks(0:1:6)
xlim([0 6])

% Legend in bottom-right (southeast) with LaTeX-style Greek letters
lgd = legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'southeast', 'FontSize', fontAxis);
set(lgd, 'Interpreter', 'tex', 'Box', 'on');

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

% Export as high-res PNG for easy "Insert" or "Paste" into Word
% exportgraphics(fig2, 'RAX_Scenario4_angular_velocity.png', 'Resolution', 600)


%% (c) Hysteresis Rod 1
fig3 = figure('Color','w', 'Units','inches', ...
              'Position',[1 1 figWidth figHeight], ...
              'Name','Hysteresis Rod 1');

% Black solid line for the B-H loop
plot(Hrod_hist(:,1), Brod_hist(:,1), '-k', 'LineWidth', 1.2)
grid on
box on

xlabel('Magnetic Field Strength H [A/m]', 'FontSize', fontLabel)
ylabel('Induced Flux Density B [T]', 'FontSize', fontLabel)

% Consistent axis formatting (Times New Roman, 12pt, Inward Ticks)
set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

% Export as high-res PNG for Word
% exportgraphics(fig3, 'RAX_Scenario4_hysteresis_rod1.png', 'Resolution', 600)

%% (d) Hysteresis Rod 2
fig4 = figure('Color','w', 'Units','inches', ...
              'Position',[6 1 figWidth figHeight], ...
              'Name','Hysteresis Rod 2');

% Black solid line for the B-H loop
plot(Hrod_hist(:,2), Brod_hist(:,2), '-k', 'LineWidth', 1.2)
grid on
box on

xlabel('Magnetic Field Strength H [A/m]', 'FontSize', fontLabel)
ylabel('Induced Flux Density B [T]', 'FontSize', fontLabel)

% Consistent axis formatting
set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

% Export as high-res PNG for Word
% exportgraphics(fig4, 'RAX_Scenario4_hysteresis_rod2.png', 'Resolution', 600)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CIRCULAR ORBITAL POSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r_eci = circularOrbitECI(r_orbit, incl, Omega0, u)
    cO = cos(Omega0); sO = sin(Omega0);
    ci = cos(incl);   si = sin(incl);
    cu = cos(u);      su = sin(u);

    r_eci = r_orbit * [ cO * cu - sO * su * ci;
                        sO * cu + cO * su * ci;
                        su * si ];
end
%%%%%%%%%%%%%%%%%%%% ECI-TO-ECEF DCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function REI = R_I2E(omega0, t)
    c = cos(omega0 * t);
    s = sin(omega0 * t);
    REI = [ c, -s, 0;
            s,  c, 0;
            0,  0, 1 ];
end
%%%%%%%%%%%%%%%%%%%% DERIVATIVE OF ECI-TO-ECEF DCM %%%%%%%%%%%%%%%%%%%%%%%%
function REIdot = R_I2E_dot(omega0, t)
    c = cos(omega0 * t);
    s = sin(omega0 * t);
    REIdot = [ -omega0 * s, -omega0 * c, 0;
                omega0 * c, -omega0 * s, 0;
                0,          0,           0 ];
end
%%%%%%%%%%%%%%%%%%%% EARTH MAGNETIC FIELD MODEL (IGRF) %%%%%%%%%%%%%%%%%%%%
function B_ecef_T = igrfECEF_Tesla(r_ecef_m, dyear, generation)
    if exist('igrfmagm', 'file') ~= 2
        error('igrfmagm not found. Aerospace Toolbox required.');
    end
    if exist('ecef2lla', 'file') ~= 2
        error('ecef2lla not found. Aerospace Toolbox required.');
    end

    lla = ecef2lla(r_ecef_m.');
    lat = lla(1);
    lon = lla(2);
    alt = lla(3);

    try
        XYZ_ned_nT = igrfmagm(alt, lat, lon, dyear, generation);
    catch
        XYZ_ned_nT = igrfmagm(alt, lat, lon, dyear);
    end

    B_ned_nT = XYZ_ned_nT(:);

    slat = sind(lat);  clat = cosd(lat);
    slon = sind(lon);  clon = cosd(lon);

    % NED -> ECEF
    R_ned2ecef = [ -slat * clon, -slon,        -clat * clon;
                   -slat * slon,  clon,        -clat * slon;
                    clat,         0,           -slat       ];

    B_ecef_nT = R_ned2ecef * B_ned_nT;
    B_ecef_T  = B_ecef_nT * 1e-9;
end
%%%%%%%%%%%%%%%%%%%% MAGNETIC HYSTERESIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bdot, guard_used] = duhemHenretty_Bdot_guarded(H, Hdot, B, Bm, Hc, Hr)
    % Tiny singularity guard only; no projection back to the major loop
    epsB = 1e-9;
    B_guard = min(max(B, -Bm + epsB), Bm - epsB);
    guard_used = abs(B_guard - B) > 0;

    x = (pi * B_guard) / (2 * Bm);
    c = cos(x);
    s = sin(x);

    if Hdot >= 0
        term = ((H + Hc) * c - Hr * s) / (2 * Hc); % ascending-field branch
    else
        term = ((Hc - H) * c + Hr * s) / (2 * Hc); % descending-field branch
    end

    % Differential Flatley/Henretty update law: dB/dt = (dB/dH) * dH/dt
    Bdot = (2 * Bm) / (Hr * pi) * (term ^ 2) * Hdot;
end
%%%%%%%%%%%%%%%% NEWTON SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the Newton solver used in the LGVI. This solver finds a solution
% to the discrete rotational momentum balance (Eq. 25) from Lee et al.'s "A 
% Lie Group Variational Integrator for the Attitude Dynamics of a
% Rigid Body with Applications to the 3D Pendulum." It iteratively 
% drives the residual equation toward zero. The Newton iteration converged in
% two iterations at every timestep for h = 0.05 s.
function F = lgviSolveF(J, h, Pi_half)
    invJ = diag(1 ./ diag(J));
    f = h * (invJ * Pi_half); % initial guess

    tol = 1e-12; % Newton solver settings
    maxIter = 12;
    
    % Theta equals magnitude of rotation increment vector (f)
    for it = 1:maxIter % Newton solver loop
        theta = norm(f);
        theta2 = theta^2;

        if theta < 1e-8
            % Approximations are used when theta is close to zero
            c1 = 1 - theta2 / 6 + theta2^2 / 120; % Taylor series expansion of sin(theta)/theta
            c2 = 1/2 - theta2 / 24 + theta2^2 / 720; % Taylor series expansion of (1-cos(theta))/theta^2
            dc1_th = -1/3 + theta2 / 30; % derivative term for c1 Jacobian
            dc2_th = -1/12 + theta2 / 180; % derivative term for c2 Jacobian
        else
            % Exact trig formulas are used when theta is not close to zero
            s = sin(theta);
            c = cos(theta);
            c1 = s / theta; % Rodrigues coefficient
            c2 = (1 - c) / theta2; % Rodrigues coefficient

            dc1_val = (theta * c - s) / theta2; % derivative of c1 with respect to theta
            dc2_val = (theta * s - 2 * (1 - c)) / (theta2 * theta); % derivative of c2 with respect to theta
            dc1_th = dc1_val / theta; % computes derivative terms used in Jacobian
            dc2_th = dc2_val / theta;
        end

        % Terms used in residual equation
        Jf = J * f;
        fxJf = cross(f, Jf);

        % Residual equation
        G = c1 * Jf + c2 * fxJf - h * Pi_half;

        if norm(G) < tol
            break;
        end
    
        % Builds skew matrices
        S_f  = skew(f);
        S_Jf = skew(Jf);

        % Jacobian construction
        Jac_A = c1 * J + dc1_th * (Jf * f.');
        Jac_B = c2 * (S_f * J - S_Jf) + dc2_th * (fxJf * f.');
        Jac   = Jac_A + Jac_B;

        % Rotation increment update
        f = f - Jac \ G;
    end

    theta = norm(f);
    S = skew(f);

    if theta < 1e-10
        c1 = 1 - theta^2 / 6; % Taylor approximation
        c2 = 1/2 - theta^2 / 24; % Taylor approximation
    else
        c1 = sin(theta) / theta; % Exact Rodrigues coefficient
        c2 = (1 - cos(theta)) / theta^2; % Exact Rodrigues coefficient
    end
    
    % Rodrigues formula used to convert rotation increment vector (f) into 
    % matrix (F).
    F = eye(3) + c1 * S + c2 * (S * S);

    % Optional data logging. This stores the number of iterations needed for
    % convergence. It is written to a text file at the end of the
    % integration loop.
    global NEWTON_ITERS
    NEWTON_ITERS(end+1) = it;
end
%%%%%%%%%%%%%%%%%%%% SKEW-SYMMETRIC MATRIX OPERATOR %%%%%%%%%%%%%%%%%%%%%%%
function S = skew(v)
    S = [   0,   -v(3),  v(2);
          v(3),    0,   -v(1);
         -v(2),  v(1),    0  ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%