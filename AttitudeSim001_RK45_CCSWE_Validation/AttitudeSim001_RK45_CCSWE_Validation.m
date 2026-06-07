%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Quaternion RK45 Validation Against CSSWE Euler Angle RK45
%% Uses Non-Tilted Dipole Field Model and Discontinuous Flatley-Henretty Hysteresis 
%% Reference Paper: "Passive Magnetic Attitude Control for CubeSat Spacecraft"

% This program will save a .mat data file with the results to your working
% directory.

% Reference Literature: Gerhardt's "Passive Magnetic Attitude Control for
% CubeSat Spacecraft" (2010)
% Simulation Specifications:
% Mass: 3.15 kg
% Configuration: 3U
% Rods per axis: 2 (4 total)
% Theta_0 = [0 0 0]
% W_0 = [10 5 5] deg/s
% Inertias: diag(0.00551 0.02552 0.02565)
% Relative Tolerance: 1E-7
% Absolute Tolerance: 1E-7
% Max timestep: none
% Rod Length: 0.095 m
% Rod Radius: 0.0005 m
% Rod Apparent Attributes: B_s = 0.027; H_c = 12; B_r = 0.004
% Altitude: 600 km
% Orbital inclination: 55 deg

clear; close all; clc;
%%%%%%%%%%%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial state y0 is given with scalar first quaternion:
y0 = [1 0 0 0 deg2rad(10) deg2rad(5) deg2rad(5)]';

% Define 14 days in seconds:
t_end = 1 * 14 * 24 * 3600;
t_span = [0, t_end];

fprintf('Starting 14-day simulation...\n');

progressTimer = tic;

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'OutputFcn', @(t,y,flag) wallClockProgress(t,y,flag,progressTimer,t_end));
[T, Y] = ode45(@myODE, t_span, y0, options);

% Post-integration normalization to guarantee unit quaternions
Y(:, 1:4) = Y(:, 1:4) ./ vecnorm(Y(:, 1:4), 2, 2);

elapsedTime = toc; 
fprintf('Simulation complete in %.2f seconds.\n', elapsedTime);
final_state = Y(end, :);
disp(T(end, :)); 
disp(final_state); 

% Store in-flight component-driven hysteresis history
mu_0 = 4*pi*1e-7;
B_hyst_history = zeros(3, length(T));
H_body_history = zeros(3, length(T));

for k = 1:length(T)
    % Extract state history at step k
    t_k = T(k);
    q_k = Y(k, 1:4)';
    w_k = Y(k, 5:7)';
    
    % Reconstruct kinematics
    B_eci_k = getGerhardtEarthField(t_k);
    B_eci_dot_k = getGerhardtEarthFieldDerivative(t_k);
    R_bi_k = quat2rotm_manual(q_k);
    B_body_k = R_bi_k' * B_eci_k;
    
    % Store the driving H-field parallel to each body axis
    H_body_history(:, k) = B_body_k / mu_0;
    
    % Drive Y-loop (axis [0;1;0])
    B_hyst_history(2,k) = calculate_B_hyst(B_body_k, B_eci_dot_k, R_bi_k, [0;1;0], w_k);
    
    % Drive Z-loop (axis [0;0;1])
    B_hyst_history(3,k) = calculate_B_hyst(B_body_k, B_eci_dot_k, R_bi_k, [0;0;1], w_k);
end


singlePos = [0.13 0.18 0.82 0.74];
topPos    = [0.13 0.58 0.82 0.34];
botPos    = [0.13 0.12 0.82 0.34];

%%%%%%%%%%%%% HYSTERESIS LOOP PLOT: Y ROD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_loop_y = figure('Name','Hysteresis Loop - Y Rod','NumberTitle','off','Color','w');
set(fig_loop_y,'Units','inches');
set(fig_loop_y,'Position',[1 1 7.2 3.8]);
set(fig_loop_y,'PaperPositionMode','auto');

axy = axes('Parent',fig_loop_y,'Position',singlePos);
plot(axy, H_body_history(2,:), B_hyst_history(2,:), 'LineWidth', 1.4);

xlabel(axy,'H [A/m]','FontSize',18);
ylabel(axy,'B_{hyst} [T]','FontSize',18);
title(axy,'Y Rod','FontSize',18);
grid(axy,'on');
box(axy,'on');

axy.LineWidth = 1.2;
axy.TickDir   = 'in';
axy.Layer     = 'top';
axy.FontSize  = 16;

%exportgraphics(fig_loop_y,'Hysteresis_Loop_Y_Rod.png','Resolution',600);

%%%%%%%%%%%%% HYSTERESIS LOOP PLOT: Z ROD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_loop_z = figure('Name','Hysteresis Loop - Z Rod','NumberTitle','off','Color','w');
set(fig_loop_z,'Units','inches');
set(fig_loop_z,'Position',[1 1 7.2 3.8]);
set(fig_loop_z,'PaperPositionMode','auto');

axz = axes('Parent',fig_loop_z,'Position',singlePos);
plot(axz, H_body_history(3,:), B_hyst_history(3,:), 'LineWidth', 1.4);
xlabel(axz,'H [A/m]','FontSize',18);
ylabel(axz,'B_{hyst} [T]','FontSize',18);
title(axz,'Z Rod','FontSize',18);
grid(axz,'on');
box(axz,'on');

axz.LineWidth = 1.2;
axz.TickDir   = 'in';
axz.Layer     = 'top';
axz.FontSize  = 16;

%exportgraphics(fig_loop_z,'Hysteresis_Loop_Z_Rod.png','Resolution',600);

%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_days = T(:) / 86400;   
Y = Y(:,:);              
idx = 1:size(Y,1);    

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);

c_wx   = [0 0.25 0.85];   % blue
c_wy   = [0.85 0.1 0.1];  % red
c_wz   = [0.93 0.69 0.13];% yellow
c_wmag = [0 0 0];         % black
c_beta = [0 0 0];         % pointing error (black)

x_ticks_days = 0:1:14;

% angular velocity plotting
fig_vel = figure('Name','Angular Velocity Damping Profile',...
                 'NumberTitle','off','Color','w');
set(fig_vel,'Units','inches');
set(fig_vel,'Position',[1 1 7.2 3.8]);
set(fig_vel,'PaperPositionMode','auto');

axv = axes('Parent',fig_vel,'Position',singlePos);
Y_deg = Y(:,5:7)*(180/pi);
w_total_deg = sqrt(sum(Y_deg.^2,2));

plot(axv,T_days(idx),Y_deg(idx,1),'Color',c_wx,'LineWidth',1.4); hold(axv,'on');
plot(axv,T_days(idx),Y_deg(idx,2),'Color',c_wy,'LineWidth',1.4);
plot(axv,T_days(idx),Y_deg(idx,3),'Color',c_wz,'LineWidth',1.4);
plot(axv,T_days(idx),w_total_deg(idx),'Color',c_wmag,'LineWidth',2.4);
yline(axv,0,'-k','HandleVisibility','off','LineWidth',1);

xlabel(axv,'Time [days]','FontSize',18);
ylabel(axv,'Angular Velocity [deg/sec]','FontSize',18);
legend(axv,{'$\omega_x$','$\omega_y$','$\omega_z$','$|\omega|$'},...
       'Interpreter','latex',...
       'Location','northeast',...
       'Box','on',...
       'FontSize',16);
grid(axv,'on');
box(axv,'on');
xlim(axv,[0 14]);
xticks(axv,x_ticks_days);
ylim(axv,[-10 15]);

axv.LineWidth = 1.2;
axv.TickDir   = 'in';
axv.Layer     = 'top';
axv.FontSize  = 16;

%exportgraphics(fig_vel,'Final_Angular_Velocity_Profile.png','Resolution',600);

% pointing error plotting
% For CSSWE, the permanent magnet is aligned with the Body X-axis [1; 0; 0].
% Pointing error is the angle between the X-axis and the local Earth magnetic field.
pointing_error = zeros(length(T), 1);
for k = 1:length(T)
    q_k = Y(k, 1:4)';
    %pos_k = circularPropagatorGerhardt(T(k));
    B_eci_k = getGerhardtEarthField(T(k));
    R_bi_k = quat2rotm_manual(q_k);
    B_body_k = R_bi_k' * B_eci_k;
    
    B_mag_k = norm(B_body_k);
    if B_mag_k > 0
        pointing_error(k) = acosd(max(-1, min(1, B_body_k(1) / B_mag_k)));
    else
        pointing_error(k) = 0;
    end
end

%% -------------------- SAVE RESULTS --------------------------------------
results_filename = fullfile(pwd, 'AttitudeSim001_Results.mat');

save(results_filename, ...
    'T', 'T_days', 'Y', ...
    'pointing_error', ...
    'B_hyst_history', 'H_body_history', ...
    'final_state', 'elapsedTime');

fprintf('Saved results to: %s\n', results_filename);
%% ------------------------------------------------------------------------
fig_point = figure('Name','Pointing Error',...
                   'NumberTitle','off','Color','w');
set(fig_point,'Units','inches');
set(fig_point,'Position',[1 1 7.2 3.8]);
set(fig_point,'PaperPositionMode','auto');

axp = axes('Parent',fig_point,'Position',singlePos);
plot(axp,T_days(idx),pointing_error(idx),'Color',c_beta,'LineWidth',1.4);

xlabel(axp,'Time [days]','FontSize',18);
ylabel(axp,'$\beta$ [degrees]','Interpreter','latex','FontSize',18);
grid(axp,'on');
box(axp,'on');
xlim(axp,[0 14]);
xticks(axp,x_ticks_days);
ylim(axp,[0 180]);

axp.LineWidth = 1.2;
axp.TickDir   = 'in';
axp.Layer     = 'top';
axp.FontSize  = 16;

%exportgraphics(fig_point,'Final_Pointing_Error.png','Resolution',600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% ODE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = myODE(t, y)

    q_in = y(1:4);

    %% Hard normalization for quaternion and converted to row vector.
    % Quaternion drift occurs in integrators like ODE45 because
    % the solver truncates quaternion values. A stabilizing method is
    % needed to prevent drift. In this function, a two-step stabilizing
    % method is implemented. A hard normalization at the start ensures that
    % a perfect unit quaternion is used for the current timestep's physics
    % calculations.
    q_norm = norm(q_in);
    q = (q_in / q_norm)'; 
    
    w = y(5:7); 

    %% Inertias.
    J = diag([0.00551 0.02552 0.02565]);
    
    %% Physical constants of hysteresis rods
    mu_0 = 4 * pi * 10^(-7);
    rod_L = 0.095; 
    rod_r = 0.0005; 
    V_hyst = pi * rod_r^2 * rod_L; % volume of one rod
    
    %% Magnetic flux density of Earth's field is generated and transformed
    %% to body frame.
    % A propagator function models a 600 km altitude orbit with inclination
    % of 55 deg. An Earth field function is used that implements a non-tilted
    % dipole Earth field model. A frame transformation function is used that 
    % converts from body to inertial frame, requiring transpose (') for 
    % inertial to body.

    B_eci = getGerhardtEarthField(t); 
    R_bi = quat2rotm_manual(q); 
    B_body = R_bi' * B_eci;
    B_eci_dot = getGerhardtEarthFieldDerivative(t);

    %% Hysteresis torque calculation.
    rod_axes = eye(3); 
    tau_hyst = zeros(3,1);
    for i = 1:3
        rod_dir = rod_axes(:,i);
        if i == 1 % Skip X-axis for rod placement
            num_rods = 0; % CSSWE used a Bar Magnet on X, no hysteresis rods.
        else
            num_rods = 2; % Two rods each for Y- and Z-axes.
        end
    
        if num_rods > 0
            % This calculates the magnetic flux density of hysteresis rods
            % using a discontinuous Henretty hysteresis model.
            B_rod = calculate_B_hyst(B_body, B_eci_dot, R_bi, rod_dir, w);
        
            % The magnetic moment is calculated using flux density and
            % rod volume.
            m_mag = (B_rod * (V_hyst * num_rods)) / mu_0;
            m_vec = m_mag * rod_dir;
        
            % Torque is calculated using T = m x B_external
            tau_hyst = tau_hyst + cross(m_vec, B_body);
        end
    end
    
    %% Kinematics: q_dot (Quaternion Rate).
    % This is the quaternion kinematic equation using scalar first
    % convention and hamilton product (body frame).
    qdot_raw = 0.5 * quatmultiply(q, [0, w']);
    
    %% Baumgarte Stabilization: Gently nudges the norm back to 1.0.
    % This is the second norm stabilizing technique.
    k_stabilize = 0.1; 
    norm_error = 1 - q_norm; % Used the original norm to calculate error.
    qdot = qdot_raw + k_stabilize * norm_error * q;
    
    %% Permanent Magnet: 0.3 Am^2 aligned with Body X-axis.
    m_bar = [0.3; 0; 0]; 
    
    %% Permanent Magnet Torque: tau = m x B.
    tau_bar = cross(m_bar, B_body);
    
    %% Dynamics: w_dot (Angular Acceleration).
    % External torques tau_hyst and tau_bar are included in calculation.
    w_col = y(5:7);
    wdot = J \ (tau_hyst + tau_bar - cross(w_col, J * w_col));
    
    %% Repack into state matrix.
    dydt = [qdot(:); wdot(:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% HAMILTON PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = quatmultiply(p, q)
    ps = p(1);   % Scalar part (p_s)
    pv = p(2:4); % Vector part (p_v)
    qs = q(1);   % Scalar part (q_s)
    qv = q(2:4); % Vector part (q_v)
    r_s = ps * qs - dot(pv, qv);
    r_v = ps * qv + qs * pv + cross(pv, qv);
    r = [r_s, r_v];
end
%%%%%%%%%%%%%%%%%%%%% QUAT2ROTM FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_bi = quat2rotm_manual(q)
    % QUAT2ROTM_MANUAL Converts a quaternion to a 3x3 rotation matrix.
    q = q / norm(q);
    w = q(1); 
    x = q(2); 
    y = q(3); 
    z = q(4);
    R_bi = [1 - 2*y^2 - 2*z^2,   2*x*y - 2*w*z,       2*x*z + 2*w*y;
            2*x*y + 2*w*z,       1 - 2*x^2 - 2*z^2,   2*y*z - 2*w*x;
            2*x*z - 2*w*y,       2*y*z + 2*w*x,       1 - 2*x^2 - 2*y^2];
end
%%%%%%%%%%%%%%%%% MAGNETIC FIELD FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B_eci = getGerhardtEarthField(t)
    % Non-tilted dipole model
    % Returns magnetic flux density vector in Tesla
    H_eq = 18.3;              % A/m
    i    = deg2rad(55.0);     % orbit inclination
    mu_0 = 4*pi*1e-7;         % H/m
    mu   = 3.986004418e14;    % m^3/s^2

    % Circular orbit mean motion at 600 km
    R  = (6371 + 600) * 1e3;      
    n  = sqrt(mu / R^3);

    u0 = 0;
    u  = n*t + u0;

    H1 = 3 * H_eq * sin(i) * cos(i) * (sin(u)^2);
    H2 = -3 * H_eq * sin(i) * sin(u) * cos(u);
    H3 = H_eq * (1 - 3 * (sin(i)^2) * (sin(u)^2));

    H_eci = [H1; H2; H3];     % A/m
    B_eci = mu_0 * H_eci;     % Tesla
end
%%%%%%%%%%%%%%%%%%% MAGNETIC HYSTERESIS FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
function B_hyst = calculate_B_hyst(B_body, B_eci_dot, R_bi, rod_axis, w)
    % B_eci_dot and R_bi input parameters are retained to provide option for 
    % full derivative testing (rotation and orbital field variation).

    B_s = 0.027; H_c = 12; B_r = 0.004; % CSSWE rod parameters
    p = (1 / H_c) * tan((pi * B_r) / (2 * B_s));
    mu_0 = 4 * pi * 10^(-7);

    H_parallel = dot(B_body / mu_0, rod_axis);
    
    % Rotation-only body-frame derivative used for Gerhardt result
    % replication:
    B_body_dot = -cross(w, B_body); 
    
    % Instead of also accounting for orbital field variation:
    %B_body_dot = R_bi' * B_eci_dot - cross(w, B_body);

    dH_dt = dot(B_body_dot / mu_0, rod_axis);

    % Hard switching based on dH/dt sign
    if dH_dt > 0
        H_shift = -H_c;
    elseif dH_dt < 0
        H_shift = +H_c;
    else
        H_shift = 0;
    end

    % Discontinuous Flatley-Henretty Hysteresis:
    B_hyst = (2 / pi) * B_s * atan(p * (H_parallel + H_shift));
end
%%%%%%%%%%%%%%%%%%%%% MAGNETIC FIELD DERIVATIVE %%%%%%%%%%%%%%%%%%%%%
function B_eci_dot = getGerhardtEarthFieldDerivative(t)
    % Time derivative of Gerhardt Eqns. (17)-(19)
    H_eq = 18.3;              % A/m
    i    = deg2rad(55.0);     % orbit inclination
    mu_0 = 4 * pi * 1e-7;     % H/m
    mu   = 3.986004418e14;    % m^3/s^2

    R  = (6371 + 600) * 1e3;
    n  = sqrt(mu / R^3);

    u0 = 0;
    u  = n * t + u0;

    H1_dot = 6 * H_eq * sin(i) * cos(i) * sin(u) * cos(u) * n;
    H2_dot = -3 * H_eq * sin(i) * (cos(u).^2 - sin(u).^2) * n;
    H3_dot = -6 * H_eq * (sin(i).^2) * sin(u) * cos(u) * n;

    H_eci_dot = [H1_dot; H2_dot; H3_dot];
    B_eci_dot = mu_0 * H_eci_dot;
end
%%%%%%%%%%%%%%%%%%%%%% WALL CLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = wallClockProgress(t, ~, ~, progressTimer, t_end)
    persistent lastUpdate
    if isempty(lastUpdate)
        lastUpdate = 0;
    end
    currentWallTime = toc(progressTimer);
    if (currentWallTime - lastUpdate) > 5
        if ~isempty(t)
            currentTime = t(end);
            percentComplete = (currentTime / t_end) * 100;
            eta = currentWallTime * (100 / percentComplete - 1);
            eta_days = floor(eta / 86400);
            eta_hours = floor(mod(eta, 86400) / 3600);
            eta_minutes = floor(mod(eta, 3600) / 60);
            eta_seconds = floor(mod(eta, 60));
            fprintf('Wall Clock: %.0fs | Progress: %.2f%% | ETA: %dd %02dh %02dm %02ds (Day %.2f of 14)\n', ...
                currentWallTime, percentComplete, eta_days, eta_hours, eta_minutes, eta_seconds, currentTime/86400);
        end
        lastUpdate = currentWallTime;
    end
    status = 0; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%