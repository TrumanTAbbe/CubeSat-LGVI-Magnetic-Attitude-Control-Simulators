%% Truman Abbe | Utah State University | truman.abbe23@gmail.com
%% 801-Trial Monte Carlo Study
%% SO(3) Lie Group Variational Integrator (LGVI) GASRATS Attitude Simulation
%% Uses IGRF-14 and the Flatley Differential Hysteresis 

%% This provides modular .mat ouput for to allow for any settling criteria 
%% evaluation, and provides a plot with evaluation of 0.75 deg/s and 15 deg requirement.

%% This runs 801 trials using Gerhardt's open circuit hysteres parameters 
%% and magnet strength of 0.09 A*m^2.
%% The goal is to confirm whether the 1U anomaly from Sim005 was an artifact 
%% of the hysteresis parameter mismatch. The trials are run for 30 day durations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GASRATS PARAMETERS USED IN THIS SIMULATION
% Baseline Angular Velocity: [10 5 5] deg/s 
% Rods per axis: 2 (4 total)
% Rod diameter: 0.0015875 m 
% Rod length: 0.070485 m 
% Apparent rod properties
% Timestep: 0.05 s
% Field model: IGRF
% Altitude: 400 km
% Inclination: 51.6 deg
% Time Date of Deployment: April 1st 2027
% m 1U   = 1.40 kg
% m 1.5U = 1.64 kg
% m 2U   = 1.76 kg
% m 3U   = 2.00 kg

%% %%%%%%%%%%% SIMULATION CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

assert(license('test','Aerospace_Toolbox') == 1, ...
    'Aerospace Toolbox required.');

assert(license('test','Distrib_Computing_Toolbox') == 1, ...
    'Parallel Computing Toolbox required.');

%%%%%%%%%%%%%% CALCULATE INERTIAL MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global J_1U J_1_5U J_2U J_3U
% Base Constants
m_core = 1.4;   
s      = 0.1;   
% 1U
I_core_val = (1/6) * m_core * s^2;
J_1U = diag([I_core_val, I_core_val, I_core_val]);
% 1.5U
m_extra_total = 0.24;         
m_ext_side    = m_extra_total / 2; 
L_ext         = 0.025;             
d             = (s/2) + (L_ext/2); 
I_mid = (1/6) * m_core * s^2;
I_ext_z  = (1/6) * m_ext_side * s^2;                  
I_ext_xy = (1/12) * m_ext_side * (s^2 + L_ext^2);     
I_zz = I_mid + 2 * I_ext_z;
I_xx = I_mid + 2 * (I_ext_xy + m_ext_side * d^2);
I_yy = I_xx;
J_1_5U = diag([I_xx, I_yy, I_zz]);
% 2U
m_extra_total = 0.36;
m_ext_side    = m_extra_total / 2; 
L_ext         = 0.05;              
d             = (s/2) + (L_ext/2);
I_mid = (1/6) * m_core * s^2;
I_ext_z  = (1/6) * m_ext_side * s^2;
I_ext_xy = (1/12) * m_ext_side * (s^2 + L_ext^2);
I_zz = I_mid + 2 * I_ext_z;
I_xx = I_mid + 2 * (I_ext_xy + m_ext_side * d^2);
I_yy = I_xx;
J_2U = diag([I_xx, I_yy, I_zz]);
% 3U
m_extra_total = 0.60;
m_ext_side    = m_extra_total / 2; 
L_ext         = 0.1;               
d             = (s/2) + (L_ext/2);
I_mid = (1/6) * m_core * s^2;
I_ext_z  = (1/6) * m_ext_side * s^2;
I_ext_xy = (1/12) * m_ext_side * (s^2 + L_ext^2);
I_zz = I_mid + 2 * I_ext_z;
I_xx = I_mid + 2 * (I_ext_xy + m_ext_side * d^2);
I_yy = I_xx;
J_3U = diag([I_xx, I_yy, I_zz]);
disp("This is J_3U");
disp(J_3U);

%%%%%%%%%%%% MONTE CARLO MASTER CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section handles the setup for 801 trials. It initializes the
% parallel computing environment, maps out the configurations, and
% generates random initial angular velocities and orientations for each
% trial.
if isempty(gcp('nocreate')), parpool; end
total_trials = 801; 
J_list = {J_1U, J_1_5U, J_2U, J_3U};
config_names = {'1U', '1.5U', '2U', '3U'};

% Generate the 801 random trial satellite configs (200 of each 1, 2, 3, 4)
base_sat_map = repelem([1, 2, 3, 4], 200); 

% Baseline 3U trial followed by the 800 MC trials
sat_config_map = [4, base_sat_map];

% Every trial uses 2 rods/axis
rod_config_map = 2*ones(1,total_trials);

results_lock_time = zeros(1, total_trials);
results_nominal   = zeros(total_trials, 1); 

Master_Inputs = struct(); 
w_all_matrix = zeros(total_trials, 3); 
q_all_matrix = zeros(total_trials, 4);
max_mag = 35; nominal_threshold = 8.660254038; % norm(5, 5, 5) deg/s
rng('shuffle', 'twister');
run_seed = rng;

for i = 1:total_trials
   if i == 1
       w_rand_deg = [10; 5; 5]; 
       q_rand = [1; 0; 0; 0];   % Identity for baseline
   else
       % rng(rng_state.Seed + i, 'twister'); 
       % Random initial velocity
       rand_mag = rand() * max_mag; 
       rand_dir = randn(3,1); rand_dir = rand_dir / norm(rand_dir);
       w_rand_deg = rand_mag * rand_dir;
       
       % Random initial orientation (uniform quaternion sampling)
       u = rand(3,1);
       q_rand = [sqrt(u(1))*cos(2*pi*u(3));
                sqrt(1-u(1))*sin(2*pi*u(2));
                sqrt(1-u(1))*cos(2*pi*u(2));
                sqrt(u(1))*sin(2*pi*u(3))];
   end
   
   % Trial configurations are stored in the Master_Inputs structure array
   results_nominal(i) = (norm(w_rand_deg) <= nominal_threshold);
   c_idx = sat_config_map(i);
   Master_Inputs(i).trialID = i;
   Master_Inputs(i).config = config_names{c_idx};
   Master_Inputs(i).rod_config = rod_config_map(i);
   Master_Inputs(i).w0_deg = w_rand_deg;
   Master_Inputs(i).q0 = q_rand; % Save q0 to struct
   w_all_matrix(i, :) = w_rand_deg';
   q_all_matrix(i, :) = q_rand';  % Save q0 to matrix
end
save('AttitudeSim007_Inputs.mat', 'Master_Inputs', 'w_all_matrix', 'q_all_matrix', 'run_seed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------- USER SETTINGS -------------------------------------
h = 0.05;                  
dtB = 1;                 
igrf_decimal_year = 2027;  
igrf_generation  = 14;     

%% -------------------- CONSTANTS -----------------------------------------
mu_earth = 3.986004418e14;   
Re_mean  = 6371e3;           
mu0      = 4*pi*1e-7;        
omega0   = 7.292115e-5;      

%% -------------------- SCENARIO PARAMETERS -------------------------------
altitude = 400e3;            
incl     = deg2rad(51.6);    
r_orbit  = Re_mean + altitude;
n_orbit  = sqrt(mu_earth / r_orbit^3);
T_orbit  = 2*pi / n_orbit;
t_end    = 30 * 86400; % 30 days        
steps    = floor(t_end / h);
t        = (0:steps)' * h;
N        = numel(t);
fprintf('GASRATS: h = %.3f s, N = %d, sim time = %.2f hr (%.2f days)\n', ...
    h, N, t(end)/3600, t(end)/86400);

% Permanent magnet
m_bar = 0.09; % [A*m^2] This candidate strength was found to be roughly near 
% a magnet strength convergence in the isolated test
mag_term = m_bar * mu0;

% Hysteresis rods for GASRATS
rod_L = 0.070485; 
rod_r = 0.0015875/2; 
Vh = pi * rod_r^2 * rod_L; 

% Open circuit parameters (Gerhardt's measurements)
Bs = 0.3; % [T]             
Br = 6.0618 * 1e-4; % [T]
Hc = 0.3381; % [A/m]
q0 = 0;
p = 2;

% Body principal axes
b1 = [1;0;0]; b2 = [0;1;0]; b3 = [0;0;1];
Omega0 = 0; u0 = 0;   

%% -------------- PRECOMPUTE ORBIT + IGRF IN EARTH FRAME ------------------
% To achieve a 801-trial Monte Carlo within a reasonable wall-clock time, 
% the IGRF-14 magnetic field is precomputed along the orbit once. 

% Create a coarse time vector for initial field sampling
t_coarse = (0:dtB:t(end))';
if abs(t_coarse(end) - t(end)) > eps(t(end))
    t_coarse = [t_coarse; t(end)];
end

% Generate a unique filename based on orbital and model parameters.
% This ensures a cache isn't loaded from a different mission profile 
% (e.g., different year or inclination).
cache_filename = sprintf( ...
    'GASRATS_IGRFcache_dtB_%s_year_%g_gen_%d_Omega_%s_u0_%s_days_%g.mat', ...
    strrep(num2str(dtB, '%.6g'), '.', 'p'), ...
    igrf_decimal_year, igrf_generation, ...
    strrep(num2str(Omega0, '%.6g'), '.', 'p'), ...
    strrep(num2str(u0, '%.6g'), '.', 'p'), ...
    t_end/86400);
cache_fullpath = fullfile(pwd, cache_filename);
recompute_igrf = false;
if exist(cache_fullpath, 'file')
    fprintf('Loading IGRF cache: %s\n', cache_fullpath);
    S = load(cache_fullpath, 't_coarse', 'H_e_coarse', 'Hdot_e_coarse', 'H_e_hist', 'Hdot_e_hist');
    
    % Validation: Ensure the cached data matches the current simulation's 
    % required number of steps (N). If N has changed (due to 'h' or 't_end'), 
    % the field is recomputed to avoid interpolation errors.
    if size(S.H_e_hist, 1) == N
        H_e_hist    = S.H_e_hist;
        Hdot_e_hist = S.Hdot_e_hist;
    else
        fprintf('WARNING: Cache file contains mismatched data (Found %d rows, Expected %d).\n', size(S.H_e_hist, 1), N);
        fprintf('Forcing a fresh IGRF recompute...\n');
        recompute_igrf = true;
    end
else
    % No existing cache found for this configuration
    recompute_igrf = true;
end

if recompute_igrf
    fprintf('Precomputing IGRF values along orbit (dtB = %.2f s)...\n', dtB);
    igrf_wall_tic = tic;
    last_progress_wall = 0;
    progress_print_interval_wall = 10;
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

            fprintf(['  IGRF precompute: %6.2f%% | k = %d / %d | ' ...
                     'sim t = %.2f hr | wall = %.0f s | ETA = %.0f s\n'], ...
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
    save(cache_fullpath, 't_coarse', 'H_e_coarse', 'Hdot_e_coarse', 'H_e_hist', 'Hdot_e_hist', '-v7.3');
        fprintf('IGRF cache saved to: %s\n', cache_fullpath);
end

%%%%%%%%%%%%%%%%%%%%% PARALLEL LGVI MONTE CARLO LOOP %%%%%%%%%%%%%%%%%%%%%%
% This section executes the high-volume simulation across available CPU cores.
% The architecture focuses on memory efficiency and real-time status reporting.

fprintf('\nStarting Parallel LGVI Integration for %d Trials...\n', total_trials);

% --- PARALLEL STATUS MESSAGING ---
% DataQueue handles asynchronous communication between workers and the main 
% thread to provide real-time console updates without slowing the integration.
dq = parallel.pool.DataQueue;
afterEach(dq, @(msg) updateMCProgress(msg));

% Initialize results array with NaN to distinguish between failures and crashes
results_lock_time = NaN(1, total_trials);
global_start_time = tic;

% --- MEMORY MANAGEMENT: DATA DOWNSAMPLING ---
save_cadence_sec = 10; 
save_stride = round(save_cadence_sec / h);
idx_save = 1:save_stride:N;
t_save = t(idx_save);

% Pre-allocate downsampled matrices for parallel efficiency
results_ortho_hist = zeros(total_trials, length(idx_save));
results_w_mag_hist = zeros(total_trials, length(idx_save));
results_theta_hist = zeros(total_trials, length(idx_save));
results_sat_guard = zeros(total_trials, 1);

parfor i = 1:total_trials
    % --- 1. EXTRACT TRIAL SPECIFICS ---
    % Retrieve hardware configuration and starting conditions for current worker
    c_idx = sat_config_map(i);
    J_current = J_list{c_idx}; % Load current inertia matrix
    invJ_current = diag(1 ./ diag(J_current)); % Invert J for momentum recovery
    num_rods_per_axis = rod_config_map(i); % Set damping rod count
    w_init_rad = deg2rad(Master_Inputs(i).w0_deg); % Convert velocity to radians
    
    % Trigger asynchronous worker-to-manager status message
    send(dq, sprintf('Trial %d START: %s, %d total rods, W0=[%.1f, %.1f, %.1f]', ...
        i, config_names{c_idx}, num_rods_per_axis*2, Master_Inputs(i).w0_deg(1), Master_Inputs(i).w0_deg(2), Master_Inputs(i).w0_deg(3)));
        
    % --- 2. INITIALIZE STATE ---
    q0 = Master_Inputs(i).q0; % Load random initial quaternion
    
    % Quaternion to Direction Cosine Matrix (R_BI)
    q0 = q0 / norm(q0); % Ensure normalization
    q0_s = q0(1); q0_v = q0(2:4);
    R_BI = eye(3) + 2*q0_s*skew(q0_v) + 2*skew(q0_v)^2; % Initial Body-to-Inertial matrix
    
    Pi   = J_current * w_init_rad; % Initial angular momentum
    Brod = [0; 0];                 % Initialize rod induced flux density [T]
    sat_guard_count = 0; 

    % --- WORKER-LOCAL DIAGNOSTICS ---
    % These arrays exist only in the worker's RAM to ensure thread safety.
    theta_worker = zeros(N, 1); % Pointing error log
    w_mag_worker = zeros(N, 1); % Velocity magnitude log
    ortho_worker = zeros(N, 1); % Orthogonality log
    
    worker_last_print_time = tic;
    
    % --- 3. THE LGVI INTEGRATION ---
    for k = 1:N
        tk     = t(k);

        % Field transformations and magnetic field evaluation:
        REI    = R_I2E(omega0, tk);
        REIdot = R_I2E_dot(omega0, tk);
        H_E    = H_e_hist(k, :).';
        Hdot_E = Hdot_e_hist(k, :).';
        H_body = R_BI * (REI.') * H_E;
        H1 = b1.' * H_body;
        H2 = b2.' * H_body;
        
        % Attitude and energy diagnostic quantities:
        Hmag = norm(H_body);
        if Hmag > 0
            theta_worker(k) = acosd(max(-1, min(1, (b3.' * H_body) / Hmag)));
        else
            theta_worker(k) = 0;
        end
        omega = invJ_current * Pi;
        w_mag_worker(k) = norm(omega) * (180/pi);
        ortho_worker(k)  = max(norm(R_BI.' * R_BI - eye(3), 'fro'), 1e-16);
        
        if k == N, break; end
        
        % Magnetic field time derivative evaluation:
        common_mat = -skew(omega) * R_BI * (REI.') + R_BI * (REIdot.');
        Hdot1 = b1.' * (common_mat * H_E + R_BI * (REI.') * Hdot_E);
        Hdot2 = b2.' * (common_mat * H_E + R_BI * (REI.') * Hdot_E);
        
        % Hysteresis rod magnetization dynamics:
        [Brod_next, sat1, sat2] = RK4_Flatley_Henretty(H1, Hdot1, H2, Hdot2, Brod, Bs, Hc, Br, q0, p, h);

        % Gerhardt's post-step clamp: check the state AFTER integration
        if Brod_next(1) > Bs
            Brod_next(1) = Bs;
            sat_guard_count = sat_guard_count + 1;
        elseif Brod_next(1) < -Bs
            Brod_next(1) = -Bs;
            sat_guard_count = sat_guard_count + 1;
        end
        
        % Same for rod 2
        if Brod_next(2) > Bs
            Brod_next(2) = Bs;
            sat_guard_count = sat_guard_count + 1;
        elseif Brod_next(2) < -Bs
            Brod_next(2) = -Bs;
            sat_guard_count = sat_guard_count + 1;
        end

        % Current magnetic torque computation:
        M_k = (mag_term) * cross(b3, H_body) + ...
              (num_rods_per_axis * Brod(1) * Vh) * cross(b1, H_body) + ...
              (num_rods_per_axis * Brod(2) * Vh) * cross(b2, H_body);
        
        % Half-step momentum update
        Pi_half = Pi + 0.5 * h * M_k;

        % LGVI rotation update
        F = lgviSolveF(J_current, h, Pi_half);
        R_BI_next = F.' * R_BI;

        if k >= N
            R_BI = R_BI_next;
            break; 
        end
        
        tkp1        = t(k+1);
        REI_next    = R_I2E(omega0, tkp1);
        H_E_next    = H_e_hist(k+1, :).';
        H_body_next = R_BI_next * (REI_next.') * H_E_next;
        
        % Future magnetic torque computation
        M_kp1 = (mag_term) * cross(b3, H_body_next) + ...
                (num_rods_per_axis * Brod_next(1) * Vh) * cross(b1, H_body_next) + ...
                (num_rods_per_axis * Brod_next(2) * Vh) * cross(b2, H_body_next);
        
        % Final momentum update
        Pi_next = F.' * Pi_half + 0.5 * h * M_kp1;
        
        % Commit
        R_BI = R_BI_next;
        Pi   = Pi_next;
        Brod = Brod_next;
        
        if toc(worker_last_print_time) > 10
            pct_done = (k / N) * 100;
            wall_now = toc(global_start_time); % Time since integration began
            send(dq, sprintf('   [Trial %d] Progress: %5.1f%% | Sim Day: %5.2f | Wall: %.1f s', ...
                i, pct_done, tk/86400, wall_now));
            worker_last_print_time = tic; 
        end
    end
    
% --- 4. LOCK DETECTION (2-Orbit Verification Window) ---
    % Find the indices where the satellite violates the capture requirements
    % Thresholds: 15 deg pointing error
    violation_idx = find(theta_worker >= 15 | w_mag_worker >= 0.75);
    
    % Requirement: 2 full orbits (~180 mins / 10800 seconds) of clean data
    % to prove the satellite is truly trapped in the field.
    required_clean_steps = 10800 / h; 
    
    if isempty(violation_idx)
        results_lock_time(i) = t(1) / 86400; % Never violated, locked from start
    elseif (N - violation_idx(end)) < required_clean_steps
        % It was still wobbling/tumbling within the last 2 orbits.
        % This identifies the "False Positives" you saw on the graph.
        results_lock_time(i) = NaN; 
    else
        % Lock occurs at the step immediately following the FINAL violation
        lock_idx = violation_idx(end) + 1;
        results_lock_time(i) = t(lock_idx) / 86400; 
    end
    
    % Send end message for progress tracking
    if isnan(results_lock_time(i))
        txt = sprintf('Trial %d END: FAILED (Wobble detected in final 2 orbits).', i);
    else
        txt = sprintf('Trial %d END: LOCKED at Day %.2f', i, results_lock_time(i));
    end
    msg_struct = struct('text', txt, 'isEnd', true, 'total', total_trials);
    send(dq, msg_struct);
    
    % Slice downsampled history out to the main workspace
    results_ortho_hist(i, :)  = ortho_worker(idx_save);
    results_w_mag_hist(i, :)  = w_mag_worker(idx_save);
    results_theta_hist(i, :) = theta_worker(idx_save);
    results_sat_guard(i) = sat_guard_count;
end % <--- End of the parfor loop

elapsedTotal = toc(global_start_time);
fprintf('\n======================================================\n');
fprintf('TOTAL INTEGRATION LOOP TIME: %.2f seconds (%.2f hours)\n', elapsedTotal, elapsedTotal/3600);
fprintf('======================================================\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveDir = pwd; 

filename = 'AttitudeSim007_Results.mat';
fullPath = fullfile(saveDir, filename);
% 1. --- SAVE MAT FILE ---
save(fullPath, 'Master_Inputs', 'results_lock_time', 'results_nominal', ...
    'w_all_matrix', 'q_all_matrix', 't_save', 'results_ortho_hist', ...
    'results_w_mag_hist', 'results_theta_hist', 'results_sat_guard', '-v7.3');
txtFilename = fullfile(saveDir, 'AttitudeSim007_Text_Report.txt');
% 2. --- GENERATE TEXT REPORT ---
fileID = fopen(txtFilename, 'w');
% Updated Header with Quaternion columns
fprintf(fileID, 'Trial,Sat_Config,Total_Rods,Nominal,Reached_Steady_State,Lock_Day,W_Mag(deg/s),Wx,Wy,Wz,q0,q1,q2,q3,Sat_Guard_Count\n');

for i = 1:total_trials
    w_i = Master_Inputs(i).w0_deg;
    q_i = Master_Inputs(i).q0;
    total_rods = Master_Inputs(i).rod_config * 2;
    
    l_day = results_lock_time(i);
    if isnan(l_day), l_day = 999.99; status = 'No'; else, status = 'Yes'; end
    if results_nominal(i) == 1, nom = 'Yes'; else, nom = 'No'; end
    
    % Print values including q0-q3
    fprintf(fileID, '%d,%s,%d,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
        i, Master_Inputs(i).config, total_rods, nom, status, l_day, ...
        norm(w_i), w_i(1), w_i(2), w_i(3), q_i(1), q_i(2), q_i(3), q_i(4), results_sat_guard(i));
end
fclose(fileID);
fprintf('Detailed Trial Log saved to: %s\n', txtFilename);

%% -------------------- MONTE CARLO MASTER PLOT ---------------------------
fprintf('Generating Master Plot (Invisible Batch Rendering)...\n');

% Create figure invisibly with the stable CPU renderer
fig = figure('Color', 'w', 'Name', 'GASRATS Monte Carlo Master Plot', ...
             'Visible', 'off', 'Renderer', 'painters');
ax = axes(fig, 'Position', [0.1 0.1 0.75 0.8]); 
hold(ax, 'on'); grid(ax, 'on');

yline(ax, 7, '-r', 'LineWidth', 2); 
xline(ax, 8.660254038, '--k', 'Nominal Limit', 'LabelVerticalAlignment', 'bottom'); 
markers = {'^', 'd', 'o', 's'}; 

% Vectorized grouping variables
w_mags = vecnorm(w_all_matrix, 2, 2);

for c = 1:4
    for r = 1:2
        if r == 1
            c_color = 'r'; 
        else
            c_color = 'b'; 
        end
        
        c_marker = markers{c};
        idx = find(sat_config_map == c & rod_config_map == r);
        if isempty(idx), continue; end
        
        l_times = results_lock_time(idx);
        success_mask = ~isnan(l_times);
        fail_mask = isnan(l_times);
        
        if any(success_mask)
            scatter(ax, w_mags(idx(success_mask)), l_times(success_mask), 50, ...
                'MarkerFaceColor', c_color, 'MarkerEdgeColor', 'k', ...
                'Marker', c_marker, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        
        if any(fail_mask)
            scatter(ax, w_mags(idx(fail_mask)), ones(sum(fail_mask),1)*30.5, 50, ...
                'Marker', 'x', 'MarkerEdgeColor', c_color, 'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
    end
end

xlabel(ax, 'Initial Total Tumble Rate (deg/s)', 'FontSize', 12);
ylabel(ax, 'Time to Stable Lock (Days)', 'FontSize', 12);
title(ax, '801 Trial LGVI Monte Carlo Verification (30-day)', 'FontSize', 14);
xlim(ax, [0 36]); ylim(ax, [0 31]);
% Clean up y-axis to match 30 days
yticks(ax, 0:5:30); 

hj = zeros(6, 1);
hj(1) = scatter(ax, NaN,NaN, 50, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', '1U (1.4kg)');
hj(2) = scatter(ax, NaN,NaN, 50, 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', '1.5U (1.64kg)');
hj(3) = scatter(ax, NaN,NaN, 50, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', '2U (1.76kg)');
hj(4) = scatter(ax, NaN,NaN, 50, 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'DisplayName', '3U (2.0kg)');
hj(5) = scatter(ax, NaN,NaN, 50, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', '2 Total Rods');
hj(6) = scatter(ax, NaN,NaN, 50, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'DisplayName', '4 Total Rods');
legend(ax, hj, 'Location', 'eastoutside');

% Export
plotPath_MC = fullfile(saveDir, '007_30Day_MC_801_LGVI_Plot.png');
exportgraphics(fig, plotPath_MC, 'Resolution', 300);
set(fig, 'Visible', 'on'); % Reveal window only after saving

fprintf('Master plot safely saved as %s\n', plotPath_MC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% WALL CLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateMCProgress(msg)
    persistent completed last_print_time;
    if isempty(completed), completed = 0; end
    if isempty(last_print_time), last_print_time = tic; end
    if ischar(msg)
        disp(msg);
        if toc(last_print_time) >= 10
            fprintf('>>> GLOBAL SUMMARY: %d / 801 Trials Fully Finished <<<\n', completed);
            last_print_time = tic; % Reset the manager's timer
        end
    elseif isstruct(msg)
        disp(msg.text);
        if isfield(msg, 'isEnd') && msg.isEnd
            completed = completed + 1;
            fprintf('\n------------------------------------------------------\n');
            fprintf('>>> OVERALL PROGRESS: %d / 801 Trials Completed (%.1f%%) <<<\n', ...
                completed, (completed/msg.total)*100);
            fprintf('------------------------------------------------------\n\n');
            last_print_time = tic;
        end
    end
end
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
    % Converts ECEF position to Geodetic (Lat/Lon/Alt) for IGRF compatibility
    lla = ecef2lla(r_ecef_m.');
    lat = lla(1); lon = lla(2); alt = lla(3);
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
    R_ned2ecef = [ -slat * clon, -slon,        -clat * clon;
                   -slat * slon,  clon,        -clat * slon;
                    clat,         0,           -slat       ];
    % Rotate vector and convert from nanoTesla to Tesla for SI consistency
    B_ecef_nT = R_ned2ecef * B_ned_nT;
    B_ecef_T  = B_ecef_nT * 1e-9;
end
%%%%%%%%%%%%%%%%%%%% MAGNETIC HYSTERESIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bdot, guard_used] = Flatley_Henretty_Bdot(H, Hdot, B, Bs, Hc, Br, q0, p)
    % Disabled in-model singularity guard; no projection back to the major loop
    % epsB = 1e-9;
    % B_guard = min(max(B, -Bs + epsB), Bs - epsB);
    % guard_used = abs(B_guard - B) > 0;

    guard_used = 0;

    k = (1 / Hc) * tan((pi * Br) / (2 * Bs));

    if Hdot >= 0
        % ascending-field branch
        Bdot1 = ...
            (q0 + (1 - q0) * ...
            ((1 / (2 * Hc)) * (H - (1 / k) * tan((pi * B) / (2 * Bs)) + Hc))^p) * ...
            ((2 * k * Bs) / pi) * cos((pi * B) / (2 * Bs))^2 * Hdot; 
        Bdot = Bdot1;
    else
        % descending-field branch
        Bdot2 = ...
            (q0 + (1 - q0) * ...
            ((1 / (2 * Hc)) * (H - (1 / k) * tan((pi * B) / (2 * Bs)) - Hc))^p) * ...
            ((2 * k * Bs) / pi) * cos((pi * B) / (2 * Bs))^2 * Hdot; 
        Bdot = Bdot2;
    end

    Bdot = Bdot(1);
end
%%%%%%%%%%%%%%%% NEWTON SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = lgviSolveF(J, h, Pi_half)
% This is the Newton solver used in the LGVI. This solver finds a solution
% to the discrete rotational momentum balance (Eq. 25) from Lee et al.'s "A 
% Lie Group Variational Integrator for the Attitude Dynamics of a
% Rigid Body with Applications to the 3D Pendulum." It iteratively 
% drives the residual equation toward zero and provides a valid rotation
% increment (F).
    invJ = diag(1 ./ diag(J));
    f = h * (invJ * Pi_half);
    tol = 1e-12; maxIter = 12;
    for it = 1:maxIter
        theta = norm(f); theta2 = theta^2;
        if theta < 1e-8
            c1 = 1 - theta2 / 6 + theta2^2 / 120;
            c2 = 1/2 - theta2 / 24 + theta2^2 / 720;
            dc1_th = -1/3 + theta2 / 30; dc2_th = -1/12 + theta2 / 180;
        else
            s = sin(theta); c = cos(theta);
            c1 = s / theta; c2 = (1 - c) / theta2;
            dc1_val = (theta * c - s) / theta2; dc2_val = (theta * s - 2 * (1 - c)) / (theta2 * theta);
            dc1_th = dc1_val / theta; dc2_th = dc2_val / theta;
        end
        Jf = J * f; fxJf = cross(f, Jf);
        G = c1 * Jf + c2 * fxJf - h * Pi_half; % Residual equation
        if norm(G) < tol, break; end
        S_f  = skew(f); S_Jf = skew(Jf);
        Jac_A = c1 * J + dc1_th * (Jf * f.'); Jac_B = c2 * (S_f * J - S_Jf) + dc2_th * (fxJf * f.');
        Jac   = Jac_A + Jac_B;
        f = f - Jac \ G; % Rotation increment update
    end
    theta = norm(f); S = skew(f);
    if theta < 1e-10
        c1 = 1 - theta^2 / 6; c2 = 1/2 - theta^2 / 24;
    else
        c1 = sin(theta) / theta; c2 = (1 - cos(theta)) / theta^2;
    end
    % Rodrigues formula used to convert rotation increment vector (f) into 
    % matrix (F).
    F = eye(3) + c1 * S + c2 * (S * S);
end
%%%%%%%%%%%%%%%%%%%% SKEW-SYMMETRIC MATRIX OPERATOR %%%%%%%%%%%%%%%%%%%%%%%
function S = skew(v)
    S = [ 0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0 ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Brod_next, sat1, sat2] = RK4_Flatley_Henretty(H1, Hdot1, H2, Hdot2, Brod, Bs, Hc, Br, q0, p, h)
    % RK4 integration for Flatley-Henretty B-dot controller
    %
    % Inputs:
    %   H1, Hdot1 - Magnetic field and derivative for axis 1
    %   H2, Hdot2 - Magnetic field and derivative for axis 2
    %   Brod      - Current B-dot values [2x1 vector]
    %   Bs, Hc, Br, q0, p - Controller parameters
    %   h         - Integration step size
    %
    % Outputs:
    %   Brod_next - Updated B-dot values [2x1 vector]
    %   sat1, sat2 - Saturation values (from final evaluation)

    Brod = Brod(:);
    
    Brod(1) = Brod(1);
    Brod(2) = Brod(2);

    % k1: slope at current point
    [B1_dot_k1, sat1_k1] = Flatley_Henretty_Bdot(H1, Hdot1, Brod(1), Bs, Hc, Br, q0, p);
    [B2_dot_k1, sat2_k1] = Flatley_Henretty_Bdot(H2, Hdot2, Brod(2), Bs, Hc, Br, q0, p);
    k1 = [B1_dot_k1(1); B2_dot_k1(1)];

    % k2: slope at midpoint using k1
    Brod_mid1 = Brod + (h/2) * k1;
    [B1_dot_k2, sat1_k2] = Flatley_Henretty_Bdot(H1, Hdot1, Brod_mid1(1), Bs, Hc, Br, q0, p);
    [B2_dot_k2, sat2_k2] = Flatley_Henretty_Bdot(H2, Hdot2, Brod_mid1(2), Bs, Hc, Br, q0, p);
    k2 = [B1_dot_k2(:); B2_dot_k2(:)];

    % k3: slope at midpoint using k2
    Brod_mid2 = Brod + (h/2) * k2;
    [B1_dot_k3, sat1_k3] = Flatley_Henretty_Bdot(H1, Hdot1, Brod_mid2(1), Bs, Hc, Br, q0, p);
    [B2_dot_k3, sat2_k3] = Flatley_Henretty_Bdot(H2, Hdot2, Brod_mid2(2), Bs, Hc, Br, q0, p);
    k3 = [B1_dot_k3(:); B2_dot_k3(:)];

    % k4: slope at end point using k3
    Brod_end = Brod + h * k3;
    [B1_dot_k4, sat1_k4] = Flatley_Henretty_Bdot(H1, Hdot1, Brod_end(1), Bs, Hc, Br, q0, p);
    [B2_dot_k4, sat2_k4] = Flatley_Henretty_Bdot(H2, Hdot2, Brod_end(2), Bs, Hc, Br, q0, p);
    k4 = [B1_dot_k4(:); B2_dot_k4(:)];

    % RK4 update: weighted average of slopes
    Brod_next = Brod + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % Return saturation values from final evaluation
    sat1 = sat1_k4;
    sat2 = sat2_k4;
end