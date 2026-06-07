%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Plotting Script for AttitudeSim002_Results.mat
%% RAX Scenario 4 SO(3) LGVI Replication

clear; close all; clc;

%% -------------------- LOAD RESULTS --------------------------------------
data_filename = fullfile(pwd, 'AttitudeSim002_Results.mat');
S = load(data_filename);

t = S.t;
orbits = S.orbits;
theta_deg = S.theta_deg;
omega_hist = S.omega_hist;
Pi_hist = S.Pi_hist;
Brod_hist = S.Brod_hist;
Hrod_hist = S.Hrod_hist;
ineq_resid_upper = S.ineq_resid_upper;
ineq_resid_lower = S.ineq_resid_lower;

fprintf('Loaded results from: %s\n', data_filename);

%% -------------------- PLOT SETTINGS -------------------------------------
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultLegendFontName','Times New Roman');

fontAxis  = 12;
fontLabel = 14;
figWidth  = 4.5;
figHeight = 3.5;

c_wx   = [0 0.25 0.85];    % blue
c_wy   = [0.85 0.1 0.1];   % red
c_wz   = [0.93 0.69 0.13]; % yellow
c_beta = [0 0 0];          % black

idx = 1:length(orbits);    % full-resolution plotting

%% -------------------- POINTING ERROR ------------------------------------
fig1 = figure('Color','w', 'Units','inches', ...
              'Position',[1 1 figWidth figHeight], ...
              'Name','Pointing Error');

plot(orbits(idx), theta_deg(idx), '-', 'Color', c_beta, 'LineWidth', 0.8)

grid on
box on

xlabel('Orbit Number', 'FontSize', fontLabel)
ylabel('Pointing Error [deg]', 'FontSize', fontLabel)

xticks(0:1:6)
xlim([0 6])
yticks(0:15:90)
ylim([0 90])

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

exportgraphics(fig1, 'RAX_Scenario-4_Pointing_Error.png', 'Resolution', 600)

%% -------------------- ANGULAR VELOCITY ----------------------------------
fig2 = figure('Color','w', 'Units','inches', ...
              'Position',[6 1 figWidth figHeight], ...
              'Name','Angular Velocity');

plot(orbits(idx), omega_hist(idx,1), '-', 'Color', c_wx, 'LineWidth', 1.2); hold on
plot(orbits(idx), omega_hist(idx,2), '-', 'Color', c_wy, 'LineWidth', 1.2)
plot(orbits(idx), omega_hist(idx,3), '-', 'Color', c_wz, 'LineWidth', 1.2)

grid on
box on

xlabel('Orbit Number', 'FontSize', fontLabel)
ylabel('Angular Velocity [rad/s]', 'FontSize', fontLabel)

xticks(0:1:6)
xlim([0 6])

lgd = legend('\omega_x', '\omega_y', '\omega_z', ...
             'Location', 'southeast', ...
             'FontSize', fontAxis);
set(lgd, 'Interpreter', 'tex', 'Box', 'on');

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

exportgraphics(fig2, 'RAX_Scenario-4_Angular_Velocity.png', 'Resolution', 600)

%% -------------------- HYSTERESIS ROD 1 ----------------------------------
fig3 = figure('Color','w', 'Units','inches', ...
              'Position',[1 1 figWidth figHeight], ...
              'Name','Hysteresis Rod 1');

plot(Hrod_hist(idx,1), Brod_hist(idx,1), '-k', 'LineWidth', 1.2)

grid on
box on

xlabel('Magnetic Field Strength H [A/m]', 'FontSize', fontLabel)
ylabel('Induced Flux Density B [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

exportgraphics(fig3, 'RAX_Scenario-4_Hysteresis_Loop_X-Rod.png', 'Resolution', 600)

%% -------------------- HYSTERESIS ROD 2 ----------------------------------
fig4 = figure('Color','w', 'Units','inches', ...
              'Position',[6 1 figWidth figHeight], ...
              'Name','Hysteresis Rod 2');

plot(Hrod_hist(idx,2), Brod_hist(idx,2), '-k', 'LineWidth', 1.2)

grid on
box on

xlabel('Magnetic Field Strength H [A/m]', 'FontSize', fontLabel)
ylabel('Induced Flux Density B [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, 'LineWidth', 1.0, 'TickDir', 'in')

exportgraphics(fig4, 'RAX_Scenario-4_Hysteresis_Loop_Y-Rod.png', 'Resolution', 600)
