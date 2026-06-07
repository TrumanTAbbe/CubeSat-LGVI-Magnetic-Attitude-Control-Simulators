%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% SO(3) Lie Group Variational Integrator (LGVI) GASRATS-3U Attitude Simulation
%% Uses IGRF-14 and Continuous Flatley-Henretty Hysteresis 

% NOTE:
% The RK45 GASRATS case uses m_bar = 0.3 A*m^2.
% The LGVI GASRATS case uses m_bar = 0.03 A*m^2.
% These values are intentionally different because the continuous
% Flatley-Henretty hysteresis model and LGVI framework did not converge
% acceptably for the 0.3 A*m^2 candidate. The 0.03 A*m^2 case is retained
% as a stable candidate for the LGVI formulation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GASRATS PARAMETERS USED IN THIS SIMULATION
% Configuration 3U
% Starting Angular Velocity: [10 5 5] deg/s 
% Rods per axis: 2 (4 total)
% Rod diameter: 0.0015875 m
% Rod length: 0.070485 m
% Apparent B_s: 0.0128 T
% Apparent B_r: 0.004 T
% Apparent H_c: 12 A/m
% Permanent magnet: 0.03 A*m^2
% Timestep: 0.05 s
% Field model: IGRF
% Altitude: 400 km
% Inclination: 51.6 deg
% Time Date of Deployment: April 1st 2027
% m 1U   = 1.40 kg
% m 1.5U = 1.64 kg
% m 2U   = 1.76 kg
% m 3U   = 2.00 kg

%%%%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Set working directory to the current script's location
if ~isdeployed
    scriptPath = fileparts(mfilename('fullpath'));
    cd(scriptPath);
end

fprintf('Working directory set to: %s\n', pwd);

% Base Constants
m_core = 1.4;   % Base mass for the center 1U cube (kg)
s      = 0.1;   % Side length of the cube face (meters)

% -------------------------------------------------------------------------
% Configuration 4: 3U CubeSat
% -------------------------------------------------------------------------
% Logic: 1.4 kg Base + 0.60 kg Extra (10 panels @ 60g each)
% Geometry: Total length 30cm -> Extensions are 10cm (0.1m) each side
m_extra_total = 0.60;
m_ext_side    = m_extra_total / 2; % 0.30 kg per side
L_ext         = 0.1;               % Length of extension (meters)
d             = (s/2) + (L_ext/2);

I_mid = (1/6) * m_core * s^2;
I_ext_z  = (1/6) * m_ext_side * s^2;
I_ext_xy = (1/12) * m_ext_side * (s^2 + L_ext^2);

I_zz = I_mid + 2 * I_ext_z;
I_xx = I_mid + 2 * (I_ext_xy + m_ext_side * d^2);
I_yy = I_xx;

J_3U = diag([I_xx, I_yy, I_zz]);

%% -------------------- USER SETTINGS -------------------------------------
h = 0.05;                  % [s] fixed step, paper value
% A coarser field is precomputed for 14 days simulation duration (compared 
% to dtb = 0.05 s) to save runtime.
dtB = 1;                   % [s] coarse IGRF precompute spacing
igrf_decimal_year = 2027;  % GASRATS deployment year
igrf_generation  = 14;     % use if supported by MATLAB; fallback included
print_every_seconds_wall = 5;

%% -------------------- CONSTANTS -----------------------------------------
mu_earth = 3.986004418e14;   % [m^3/s^2]
Re_mean  = 6371e3;           % [m]
mu0      = 4*pi*1e-7;        % [N/A^2]
omega0   = 7.292115e-5;      % [rad/s] Earth rotation rate

%% -------------------- GASRATS-3U PARAMETERS -----------------------------
altitude = 400e3;            % [m] GASRATS
incl     = deg2rad(51.6);    % [rad] GASRATS

r_orbit  = Re_mean + altitude;
n_orbit  = sqrt(mu_earth / r_orbit^3);
T_orbit  = 2*pi / n_orbit;

t_end = 1 * 14 * 24 * 3600;  % 14 days in seconds
steps    = floor(t_end / h);
t        = (0:steps)' * h;
N        = numel(t);

fprintf('GASRATS: h = %.3f s, N = %d, sim time = %.2f hr (%.2f days)\n', ...
    h, N, t(end)/3600, t(end)/86400);

% Inertia matrix
J    = J_3U;
invJ = diag(1 ./ diag(J));

% Permanent magnet: 0.03 A*m^2 is candidate for LGVI framework
m_bar = 0.03;           % [A*m^2]
mag_term = m_bar * mu0; % Equivalent to Bp*Vp in the torque equation

% Hysteresis rods
rod_L = 0.070485; %length
rod_r = 0.0015875/2; %radius
Vh = pi * rod_r^2 * rod_L;

% Apparent B_s: 0.0128 T
% Apparent B_r: 0.004 T
% Apparent H_c: 12 A/m
Bm = 0.0128;              % [T]
Hc = 12;              % [A/m]
Br = 0.004;
Hr = Hc * tan(pi * Br / (2 * Bm)); % Result: ~6.414 A/m

% Body principal axes
b1 = [1;0;0];
b2 = [0;1;0];
b3 = [0;0;1];

%% -------------------- GASRATS INITIAL CONDITIONS ------------------------
R_BI = eye(3);      
omega = deg2rad([10; 5; 5]);   % [rad/s], GASRATS tumble
Pi    = J * omega;
Brod = [0; 0];                 % [B1; B2]

%% -------------------- FIXED INITIAL ORBIT GEOMETRY ----------------------
Omega0 = 0;   % [rad] initial RAAN
u0     = 0;   % [rad] initial argument of latitude
fprintf('Using fixed initial orbit geometry: Omega0 = %.3f rad, u0 = %.3f rad\n', Omega0, u0);

%% -------------------- PRECOMPUTE ORBIT + IGRF IN EARTH FRAME ------------
t_coarse = (0:dtB:t(end))';
if abs(t_coarse(end) - t(end)) > eps(t(end))
    t_coarse = [t_coarse; t(end)];
end
cache_filename = fullfile(pwd, sprintf('GASRATS_IGRFcache_dtB_%s_year_%g_gen_%d_Omega_%s_u0_%s.mat', ...
    strrep(num2str(dtB, '%.6g'), '.', 'p'), ...
    igrf_decimal_year, ...
    igrf_generation, ...
    strrep(num2str(Omega0, '%.6g'), '.', 'p'), ...
    strrep(num2str(u0, '%.6g'), '.', 'p')));
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

    igrf_wall_tic = tic;
    last_progress_wall = 0;

    NC = numel(t_coarse);
    H_e_coarse = zeros(NC, 3);

    % Print progress about every 10 wall-clock seconds
    progress_print_interval_wall = 10;

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

        % Progress marker
        wall_elapsed = toc(igrf_wall_tic);
        if (wall_elapsed - last_progress_wall >= progress_print_interval_wall) || (k == NC)
            frac_done = k / NC;
            if frac_done > 0
                est_total_wall = wall_elapsed / frac_done;
                est_remaining_wall = est_total_wall - wall_elapsed;
            else
                est_remaining_wall = NaN;
            end

            fprintf('  IGRF precompute: %6.2f%% | k = %d / %d | sim t = %.2f hr | wall = %.0f s | ETA = %.0f s\n', ...
                100 * frac_done, k, NC, tk / 3600, wall_elapsed, est_remaining_wall);

            last_progress_wall = wall_elapsed;
        end
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

    % Interpolate coarse data onto full simulation grid
    H_e_hist    = interp1(t_coarse, H_e_coarse,    t, 'linear', 'extrap');
    Hdot_e_hist = interp1(t_coarse, Hdot_e_coarse, t, 'linear', 'extrap');

    total_igrf_wall = toc(igrf_wall_tic);
    fprintf('IGRF precompute finished in %.1f s (%.2f min). Saving cache...\n', ...
        total_igrf_wall, total_igrf_wall / 60);

    save(cache_filename, 't_coarse', 'H_e_coarse', 'Hdot_e_coarse', 'H_e_hist', 'Hdot_e_hist', '-v7.3');
end

%% -------------------- INITIAL GEOMETRY DIAGNOSTIC --------------------
REI0 = R_I2E(omega0, t(1));
H_body0 = R_BI * (REI0.') * H_e_hist(1, :).';
theta0 = acosd(max(-1, min(1, (b3.' * H_body0) / norm(H_body0))));
H10 = b1.' * H_body0;
H20 = b2.' * H_body0;

fprintf('Initial pointing error after cache/load = %.6e deg\n', theta0);
fprintf('Initial H1 after cache/load = %.6e A/m\n', H10);
fprintf('Initial H2 after cache/load = %.6e A/m\n', H20);

%% -------------------- DATA STORAGE --------------------
theta_deg  = zeros(N, 1);
omega_hist = zeros(N, 3);
Pi_hist    = zeros(N, 3);
Brod_hist  = zeros(N, 2);
Hrod_hist  = zeros(N, 2);

ineq_resid_upper = zeros(N, 2);
ineq_resid_lower = zeros(N, 2);
sat_guard_count = 0;

%%%%%%%%%%%%%%%%%%%%%%% STEADY-STATE TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steady_state_threshold_w = 0.75;   % [deg/s]
steady_state_threshold_theta = 15; % [deg]

steady_state_time = NaN;
steady_state_index = NaN;

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
    M_k = mag_term * cross(b3, H_body) + ...
          (2 * Brod(1) * Vh) * cross(b1, H_body) + ...
          (2 * Brod(2) * Vh) * cross(b2, H_body);

    % Newton solver used to find rotation increment F
    Pi_half = Pi + 0.5 * h * M_k; % first half-kick
    F = lgviSolveF(J, h, Pi_half);

    % For R_BI satisfying Rdot + S(w)R = 0
    R_BI_next = F.' * R_BI; % drift

    tkp1     = t(k+1);
    REI_next = R_I2E(omega0, tkp1);
    H_E_next = H_e_hist(k+1, :).';
    H_body_next = R_BI_next * (REI_next.') * H_E_next;
    
    % Next-step torque calculation
    M_kp1 = mag_term * cross(b3, H_body_next) + ...
            (2 * Brod_next(1) * Vh) * cross(b1, H_body_next) + ...
            (2 * Brod_next(2) * Vh) * cross(b2, H_body_next);

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

%%%%%%%%%%%%%%%%%%%%%%% STEADY-STATE DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega_mag_deg = vecnorm(omega_hist, 2, 2) * (180/pi);

violation_idx = find(theta_deg >= steady_state_threshold_theta | ...
                     omega_mag_deg >= steady_state_threshold_w);

if isempty(violation_idx)
    steady_state_index = 1;
    steady_state_time = t(1);
else
    steady_state_index = violation_idx(end) + 1;

    if steady_state_index <= N
        steady_state_time = t(steady_state_index);
    else
        steady_state_index = NaN;
        steady_state_time = NaN;
    end
end

fprintf('\n================ STEADY-STATE REPORT ================\n');

if ~isnan(steady_state_time)
    fprintf('Steady-state achieved at:\n');
    fprintf('  t = %.2f seconds\n', steady_state_time);
    fprintf('  t = %.4f hours\n', steady_state_time / 3600);
    fprintf('  t = %.4f days\n', steady_state_time / 86400);
else
    fprintf('Steady-state was NOT achieved during simulation.\n');
end

fprintf('=====================================================\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir = pwd;
filename = 'AttitudeSim004_Results.mat';
fullPath = fullfile(saveDir, filename);
% Downsample saved data to reduce .mat file size
save_stride = 20;   % h = 0.05 s, so 20 steps = 1 second

idx_save = 1:save_stride:N;

t       = t(idx_save);
omega_hist    = omega_hist(idx_save,:);
theta_deg    = theta_deg(idx_save);
Brod_hist     = Brod_hist(idx_save,:);
Hrod_hist     = Hrod_hist(idx_save,:);

save(fullPath, ...
     't', 'omega_hist', 'theta_deg', ...
     'Brod_hist', 'Hrod_hist', ...
     'J', 'h', 'save_stride', ...
     '-v7.3');

fprintf('--- DATA SECURED TO DISK: %s ---\n', fullPath);

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

N_plot = length(t);
days = t / 86400;

% Downsample dense history arrays for plotting only
plot_stride = 200;
idx_plot = 1:plot_stride:N_plot;

% Global formatting
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
fontAxis  = 14;
fontTitle = 15;
fontLabel = 15;

%% (a) Pointing Error
fig1 = figure('Color','w', ...
              'Units','normalized', ...
              'Position',[0.05 0.08 0.9 0.85], ...
              'Name','Pointing Error');

plot(days, theta_deg, 'k', 'LineWidth', 1.6)

grid on
box on
xlabel('Time [days]','FontSize',fontLabel)
ylabel('Pointing Error [deg]','FontSize',fontLabel)
xticks(0:1:14)
xlim([0 14])
set(gca,'FontSize',fontAxis,'LineWidth',1.0)
%exportgraphics(fig1,'GASRATS_from_RAX_pointing_error.pdf','ContentType','vector','Resolution',300)

%% (b) Angular Velocity
omega_hist_deg = omega_hist * (180/pi);

fig2 = figure('Color','w', ...
              'Units','normalized', ...
              'Position',[0.05 0.08 0.9 0.85], ...
              'Name','Angular Velocity');

plot(days, omega_hist_deg(:,1), 'b', 'LineWidth', 1.6); hold on
plot(days, omega_hist_deg(:,2), 'r', 'LineWidth', 1.6)
plot(days, omega_hist_deg(:,3), 'Color', [0.93 0.69 0.13], 'LineWidth', 1.6)

grid on
box on
xlabel('Time [days]','FontSize',fontLabel)
ylabel('Angular Velocity [deg/s]','FontSize',fontLabel)
xticks(0:1:14)
xlim([0 14])
legend('wx','wy','wz','Location','best','FontSize',11)
set(gca,'FontSize',fontAxis,'LineWidth',1.0)
%exportgraphics(fig2,'GASRATS_from_RAX_angular_velocity.pdf','ContentType','vector','Resolution',300)

%% (c) Hysteresis Rod 1
fig3 = figure('Color','w', ...
              'Units','normalized', ...
              'Position',[0.05 0.08 0.9 0.85], ...
              'Name','Hysteresis Rod 1');

plot(Hrod_hist(idx_plot,1), Brod_hist(idx_plot,1), 'k', 'LineWidth', 1.0)

grid on
box on
xlabel('Magnetic Field Strength H [A/m]','FontSize',fontLabel)
ylabel('Induced Flux Density B [T]','FontSize',fontLabel)
set(gca,'FontSize',fontAxis,'LineWidth',1.0)
%exportgraphics(fig3,'GASRATS_from_RAX_hysteresis_rod1.pdf','ContentType','vector','Resolution',300)

%% (d) Hysteresis Rod 2
fig4 = figure('Color','w', ...
              'Units','normalized', ...
              'Position',[0.05 0.08 0.9 0.85], ...
              'Name','Hysteresis Rod 2');

plot(Hrod_hist(idx_plot,2), Brod_hist(idx_plot,2), 'k', 'LineWidth', 1.0)

grid on
box on
xlabel('H [A/m]','FontSize',fontLabel)
ylabel('Induced Flux Density B [T]','FontSize',fontLabel)
set(gca,'FontSize',fontAxis,'LineWidth',1.0)
%exportgraphics(fig4,'GASRATS_from_RAX_rod2.pdf','ContentType','vector','Resolution',300)
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
    % Converts ECEF position to Geodetic (Lat/Lon/Alt) for IGRF compatibility
    lla = ecef2lla(r_ecef_m.');
    lat = lla(1);
    lon = lla(2);
    alt = lla(3);
    % Query IGRF model for magnetic field in North-East-Down (NED) coordinates
    try
        XYZ_ned_nT = igrfmagm(alt, lat, lon, dyear, generation);
    catch
        XYZ_ned_nT = igrfmagm(alt, lat, lon, dyear);
    end
    B_ned_nT = XYZ_ned_nT(:);
    % Construct rotation matrix to transform B-field from local NED to global ECEF
    slat = sind(lat);  clat = cosd(lat);
    slon = sind(lon);  clon = cosd(lon);
    % NED -> ECEF
    R_ned2ecef = [ -slat * clon, -slon,        -clat * clon;
                   -slat * slon,  clon,        -clat * slon;
                    clat,         0,           -slat       ];
    % Rotate vector and convert from nanoTesla to Tesla for SI consistency
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
% drives the residual equation toward zero and provides a valid rotation
% increment (F). The Newton iteration converged in two iterations at every 
% timestep for h = 0.05 s.
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
end
%%%%%%%%%%%%%%%%%%%% SKEW-SYMMETRIC MATRIX OPERATOR %%%%%%%%%%%%%%%%%%%%%%%
function S = skew(v)
    S = [   0,   -v(3),  v(2);
          v(3),    0,   -v(1);
         -v(2),  v(1),    0  ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%