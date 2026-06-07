%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Plotting For SO(3) LGVI GASRATS-3U Attitude Simulation From Cross-Framework Verification Study

clear; close all; clc;

data_filename = fullfile(pwd, 'AttitudeSim004_Results.mat');
S = load(data_filename);

fprintf('Loaded results from: %s\n', data_filename);

t = S.t(:);
omega_hist = S.omega_hist;
theta_deg = S.theta_deg(:);
Brod_hist = S.Brod_hist;
Hrod_hist = S.Hrod_hist;

if size(omega_hist,1) ~= length(t)
    error('length(t) does not match the number of rows in omega_hist.');
end

if length(theta_deg) ~= length(t)
    error('length(t) does not match length(theta_deg).');
end

if size(Brod_hist,1) ~= length(t) || size(Hrod_hist,1) ~= length(t)
    error('length(t) does not match the hysteresis history arrays.');
end

days = t / 86400;
omega_hist_deg = omega_hist * (180/pi);
omega_mag_deg = vecnorm(omega_hist_deg, 2, 2);

maxPlotPts = 200000;
plotStep = max(1, floor(length(t) / maxPlotPts));
idx = 1:plotStep:length(t);

T_orbit_sec = 5556;
is_stable = (omega_mag_deg < 0.5) & (theta_deg < 10);
unstable_indices = find(~is_stable);

if isempty(unstable_indices)
    steady_state_idx = 1;
else
    last_unstable_idx = unstable_indices(end);
    time_stable_at_end = t(end) - t(last_unstable_idx);

    if time_stable_at_end >= T_orbit_sec
        steady_state_idx = last_unstable_idx + 1;
    else
        steady_state_idx = -1;
    end
end

fprintf('\n======================================================\n');
if steady_state_idx > 0
    fprintf('>>> PERMANENT STEADY-STATE CONFIRMED <<<\n');
    fprintf('Time of Final Capture: %.2f hours (Day %.2f)\n', ...
            t(steady_state_idx)/3600, t(steady_state_idx)/86400);
else
    fprintf('>>> PERMANENT STEADY-STATE NOT REACHED <<<\n');
    fprintf('Last recorded spike was at Day %.2f\n', t(unstable_indices(end))/86400);
end
fprintf('======================================================\n\n');

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');

figWidth = 7.2;
figHeight = 3.8;
fontAxis = 21;
fontLabel = 22;

set(groot, 'defaultAxesFontSize', fontAxis);
set(groot, 'defaultTextFontSize', fontAxis);
set(groot, 'defaultLegendFontSize', fontAxis);

c_wx = [0 0.25 0.85];
c_wy = [0.85 0.1 0.1];
c_wz = [0.93 0.69 0.13];
c_wmag = [0 0 0];
c_beta = [0 0 0];
c_ss = [0 0.5 0];

t_max_days = max(days);
x_ticks_days = 0:2:ceil(t_max_days);

fig_point = figure('Color','w', 'Units','inches', ...
                   'Position',[1 1 figWidth figHeight], ...
                   'Name','Pointing Error');

plot(days(idx), theta_deg(idx), '-', 'Color', c_beta, 'LineWidth', 0.8);
hold on

yline(10, '--', 'Color', c_ss, 'LineWidth', 2, 'HandleVisibility', 'off');

text(t_max_days - 0.1, 11.5, '10° pointing requirement', ...
     'Color', c_ss, ...
     'FontSize', fontAxis-2, ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom');

grid on
box on

xlabel('Time [days]', 'FontSize', fontLabel)

yl_p = ylabel('Pointing Error [deg]', 'FontSize', fontLabel);
yl_p.Units = 'normalized';
yl_p.Position(1) = -0.09;

xticks(x_ticks_days)
xlim([0 t_max_days])
ylim([0 90])
yticks(0:15:90)

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in', ...
         'XTickLabelRotation', 0)

pointing_filename = fullfile(pwd, 'GASRATS_LGVI_Pointing_Error.png');
exportgraphics(fig_point, pointing_filename, 'Resolution', 600)

fprintf('Saved pointing error figure to: %s\n', pointing_filename);

fig_vel = figure('Color','w', 'Units','inches', ...
                 'Position',[1 1 figWidth figHeight], ...
                 'Name','Angular Velocity');

plot(days(idx), omega_hist_deg(idx,1), '-', 'Color', c_wx, 'LineWidth', 1.2);
hold on

plot(days(idx), omega_hist_deg(idx,2), '-', 'Color', c_wy, 'LineWidth', 1.2);
plot(days(idx), omega_hist_deg(idx,3), '-', 'Color', c_wz, 'LineWidth', 1.2);
plot(days(idx), omega_mag_deg(idx), '-', 'Color', c_wmag, 'LineWidth', 1.2);

yline(0, '-k', 'HandleVisibility', 'off', 'LineWidth', 1);
yline(0.5, '--', 'Color', c_ss, 'LineWidth', 2, 'HandleVisibility', 'off');

text(t_max_days - 0.1, 0, '0.5°/s settling requirement', ...
     'Color', c_ss, ...
     'FontSize', fontAxis-2, ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'top', ...
     'Clipping', 'off');

grid on
box on

xlabel('Time [days]', 'FontSize', fontLabel)

yl_v = ylabel('Angular Velocity [deg/sec]', 'FontSize', fontLabel);
yl_v.Units = 'normalized';
yl_v.Position(1) = -0.09;
yl_v.Position(2) = 0.47;

xticks(x_ticks_days)
xlim([0 t_max_days])
ylim([-10 15])

lgd = legend('$\omega_x$', '$\omega_y$', '$\omega_z$', '$\|\omega\|$', ...
             'Orientation', 'horizontal', ...
             'Location', 'north', ...
             'FontSize', fontAxis);

set(lgd, 'Interpreter', 'latex', 'Box', 'on');

drawnow

pos = lgd.Position;
pos(1) = 1 - pos(3) - 0.11;
lgd.Position = pos;

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in', ...
         'XTickLabelRotation', 0)

velocity_filename = fullfile(pwd, 'GASRATS_LGVI_Angular_Velocity.png');
exportgraphics(fig_vel, velocity_filename, 'Resolution', 600)

fprintf('Saved angular velocity figure to: %s\n', velocity_filename);

fig_hyst1 = figure('Color','w', 'Units','inches', ...
                   'Position',[1 1 figWidth figHeight], ...
                   'Name','Hysteresis Rod 1');

plot(Hrod_hist(idx,1), Brod_hist(idx,1), '-k', 'LineWidth', 1.2);

grid on
box on

xlabel('H [A/m]', 'FontSize', fontLabel)
ylabel('B_{hyst} [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in')

hyst1_filename = fullfile(pwd, 'GASRATS_LGVI_Hysteresis_Loop_X-Rod.png');
exportgraphics(fig_hyst1, hyst1_filename, 'Resolution', 600)

fprintf('Saved hysteresis rod 1 figure to: %s\n', hyst1_filename);

fig_hyst2 = figure('Color','w', 'Units','inches', ...
                   'Position',[1 1 figWidth figHeight], ...
                   'Name','Hysteresis Rod 2');

plot(Hrod_hist(idx,2), Brod_hist(idx,2), '-k', 'LineWidth', 1.2);

grid on
box on

xlabel('H [A/m]', 'FontSize', fontLabel)
ylabel('B_{hyst} [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in')

hyst2_filename = fullfile(pwd, 'GASRATS_LGVI_Hysteresis_Loop_Y-Rod.png');
exportgraphics(fig_hyst2, hyst2_filename, 'Resolution', 600)

fprintf('Saved hysteresis rod 2 figure to: %s\n', hyst2_filename);