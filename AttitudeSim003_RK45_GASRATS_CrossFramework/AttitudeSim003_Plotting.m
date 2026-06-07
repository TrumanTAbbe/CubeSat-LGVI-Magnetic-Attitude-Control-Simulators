%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Plotting For Quaternion RK45 GASRATS-3U Attitude Simulation From Cross-Framework Verification Study

clear; close all; clc;

data_filename = fullfile(pwd, 'AttitudeSim003_Results.mat');
S = load(data_filename);

fprintf('Loaded results from: %s\n', data_filename);

T = S.T(:);
Y = S.Y;
pointing_error = S.pointing_error(:);
H_body_history = S.H_body_history;
B_hyst_history = S.B_hyst_history;

if size(Y,1) ~= length(T) && size(Y,2) == length(T)
    Y = Y.';
end

if size(Y,1) ~= length(T)
    error('length(T) does not match the number of rows in Y.');
end

if size(Y,2) < 7
    error('Y must contain at least 7 columns.');
end

if length(pointing_error) ~= length(T)
    error('pointing_error must have the same length as T.');
end

if size(H_body_history,2) ~= length(T) || size(B_hyst_history,2) ~= length(T)
    error('H_body_history and B_hyst_history must have length(T) columns.');
end

T_days = T / 86400;
deg_per_rad = 180/pi;

maxPlotPts = 200000;
plotStep = max(1, floor(length(T) / maxPlotPts));
idx = 1:plotStep:length(T);

T_orbit_sec = 5556;
w_total_deg = sqrt(sum(Y(:,5:7).^2, 2)) * deg_per_rad;

is_stable = (w_total_deg < 0.5) & (pointing_error < 10);
unstable_indices = find(~is_stable);

if isempty(unstable_indices)
    steady_state_idx = 1;
else
    last_unstable_idx = unstable_indices(end);
    time_stable_at_end = T(end) - T(last_unstable_idx);

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
            T(steady_state_idx)/3600, T(steady_state_idx)/86400);
else
    fprintf('>>> PERMANENT STEADY-STATE NOT REACHED <<<\n');
    fprintf('Last recorded spike was at Day %.2f\n', T(unstable_indices(end))/86400);
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

t_max_days = max(T_days);
x_ticks_days = 0:2:ceil(t_max_days);

fig_point = figure('Color','w', 'Units','inches', ...
                   'Position',[1 1 figWidth figHeight], ...
                   'Name','Pointing Error');

plot(T_days(idx), pointing_error(idx), '-', 'Color', c_beta, 'LineWidth', 0.8);
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

pointing_filename = fullfile(pwd, 'GASRATS_QRK45_Pointing_Error.png');
exportgraphics(fig_point, pointing_filename, 'Resolution', 600)

fprintf('Saved pointing error figure to: %s\n', pointing_filename);

wx_deg_plot = Y(idx,5) * deg_per_rad;
wy_deg_plot = Y(idx,6) * deg_per_rad;
wz_deg_plot = Y(idx,7) * deg_per_rad;
w_total_deg_plot = sqrt(Y(idx,5).^2 + Y(idx,6).^2 + Y(idx,7).^2) * deg_per_rad;

fig_vel = figure('Color','w', 'Units','inches', ...
                 'Position',[1 1 figWidth figHeight], ...
                 'Name','Angular Velocity');

plot(T_days(idx), wx_deg_plot, '-', 'Color', c_wx, 'LineWidth', 1.2);
hold on

plot(T_days(idx), wy_deg_plot, '-', 'Color', c_wy, 'LineWidth', 1.2);
plot(T_days(idx), wz_deg_plot, '-', 'Color', c_wz, 'LineWidth', 1.2);
plot(T_days(idx), w_total_deg_plot, '-', 'Color', c_wmag, 'LineWidth', 1.2);

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

velocity_filename = fullfile(pwd, 'GASRATS_QRK45_Angular_Velocity.png');
exportgraphics(fig_vel, velocity_filename, 'Resolution', 600)

fprintf('Saved angular velocity figure to: %s\n', velocity_filename);

fig_hyst_x = figure('Color','w', 'Units','inches', ...
                    'Position',[1 1 figWidth figHeight], ...
                    'Name','Hysteresis X Rod');

plot(H_body_history(1,idx), B_hyst_history(1,idx), '-k', 'LineWidth', 1.2);

grid on
box on

xlabel('H [A/m]', 'FontSize', fontLabel)
ylabel('B_{hyst} [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in')

hyst_x_filename = fullfile(pwd, 'GASRATS_QRK45_Hysteresis_X-Rod.png');
exportgraphics(fig_hyst_x, hyst_x_filename, 'Resolution', 600)

fprintf('Saved X-rod hysteresis figure to: %s\n', hyst_x_filename);

fig_hyst_y = figure('Color','w', 'Units','inches', ...
                    'Position',[1 1 figWidth figHeight], ...
                    'Name','Hysteresis Y Rod');

plot(H_body_history(2,idx), B_hyst_history(2,idx), '-k', 'LineWidth', 1.2);

grid on
box on

xlabel('H [A/m]', 'FontSize', fontLabel)
ylabel('B_{hyst} [T]', 'FontSize', fontLabel)

set(gca, 'FontSize', fontAxis, ...
         'LineWidth', 1.0, ...
         'TickDir', 'in')

hyst_y_filename = fullfile(pwd, 'GASRATS_QRK45_Hysteresis_Y-Rod.png');
exportgraphics(fig_hyst_y, hyst_y_filename, 'Resolution', 600)

fprintf('Saved Y-rod hysteresis figure to: %s\n', hyst_y_filename);