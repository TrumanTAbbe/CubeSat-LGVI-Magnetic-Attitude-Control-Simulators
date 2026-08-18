%% Visualize Monte Carlo Initial Conditions
% Using Modified Rodrigues Parameters (MRP) to map quaternions to 3D
% Plots vectors from origin to each point
% Re-evaluates settling criteria from theta_hist (15-deg pointing error ONLY)
% FAILED TRIALS: plotted as black shapes on the 30-day spherical shell
% MRP SHADOW FIX: Ensures q0 >= 0 so MRP maps to unit ball
% COLOR: Jet gradient (Red→Orange→Yellow→Green→Blue→Purple) based on initial velocity
% SHAPE: Triangle=1U, Diamond=1.5U, Circle=2U, Square=3U

clear; clc; close all;

%% LOAD RESULTS
results_file = 'AttitudeSim007_Results.mat';
if ~exist(results_file, 'file')
    error('Results file not found. Please ensure %s is in the current directory.', results_file);
end

% Load data including theta_hist
load(results_file, 'Master_Inputs', 'q_all_matrix', 'results_theta_hist', 't_save');

%% RE-EVALUATE SETTLING CRITERIA FROM THETA_HIST
% Criteria: 15-deg pointing error ONLY
% Must stay below 15 deg for 2 full orbits (10800 seconds)

total_trials = length(Master_Inputs);
results_lock_time = NaN(1, total_trials);

% 2 orbits verification window (10800 seconds)
h = 0.05;
save_cadence_sec = 10;
save_stride = round(save_cadence_sec / h);
required_save_steps = round((10800 / h) / save_stride);

for i = 1:total_trials
    theta_hist = results_theta_hist(i, :);
    violation_idx = find(theta_hist >= 15);
    
    if isempty(violation_idx)
        % Never violated - locked from start
        results_lock_time(i) = t_save(1) / 86400;
    elseif (length(theta_hist) - violation_idx(end)) < required_save_steps
        % Still violating in last 2 orbits - FAILED
        results_lock_time(i) = NaN;
    else
        % Lock occurs after final violation
        lock_idx = violation_idx(end) + 1;
        if lock_idx <= length(t_save)
            lock_time_days = t_save(lock_idx) / 86400;
            % Cap at 30 days (simulation duration)
            if lock_time_days > 30
                lock_time_days = NaN;
            end
            results_lock_time(i) = lock_time_days;
        else
            results_lock_time(i) = NaN;
        end
    end
end

%% PROCESS DATA
q_vectors = zeros(total_trials, 3);
config_labels = cell(total_trials, 1);
success_mask = true(total_trials, 1);
w_mag_init = zeros(total_trials, 1);

% Get initial angular velocity magnitudes
for i = 1:total_trials
    w_mag_init(i) = norm(Master_Inputs(i).w0_deg);
end

% Normalize velocity for color mapping (0 to 35 deg/s)
v_min = 0;
v_max = 35;
v_norm = (w_mag_init - v_min) / (v_max - v_min);
v_norm = max(0, min(1, v_norm));  % Clamp to [0,1]

% Create color map using JET: Red(0) -> Orange -> Yellow -> Green -> Blue -> Purple(35)
jet_colors = jet(256);
color_idx = round(v_norm * 255) + 1;
color_idx = max(1, min(256, color_idx));
color_map = jet_colors(color_idx, :);

for i = 1:total_trials
    q = q_all_matrix(i, :);
    
    % --- MRP SHADOW FIX ---
    % If q0 is negative, negate the quaternion to ensure q0 >= 0
    if q(1) < 0
        q = -q;
    end
    
    q0 = q(1);
    q_vec = q(2:4);
    g = q_vec / (1 + q0);
    
    % Get direction as UNIT vector from MRP
    g_norm = norm(g);
    if g_norm > 0
        direction = g / g_norm;
    else
        direction = [1; 0; 0];
    end
    
    lock_time = results_lock_time(i);
    if isnan(lock_time)
        success_mask(i) = false;
        % Failed trials go on the 30-day shell
        q_vectors(i, :) = direction * 30;
    else
        success_mask(i) = true;
        % Ensure lock time doesn't exceed 30 days
        if lock_time > 30
            lock_time = 30;
        end
        % Position = direction × lock_time (magnitude = lock_time)
        q_vectors(i, :) = direction * lock_time;
    end
    
    config_labels{i} = Master_Inputs(i).config;
end

%% CREATE 3D VECTOR PLOT
figure('Color', 'w', 'Name', 'GASRATS Full Attitude Sphere (MRP)', ...
       'Position', [100, 100, 900, 700]);


set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLegendFontName', 'Times New Roman');


ax = axes('Position', [0.08, 0.1, 0.6, 0.8]);
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
hold(ax, 'on');
grid(ax, 'on');
set(ax, 'GridColor', [0.5, 0.5, 0.5]);
set(ax, 'GridAlpha', 0.5);
set(ax, 'GridLineStyle', '-');
set(ax, 'LineWidth', 0.5);
view(ax, 45, 20);

% --- DRAW SPHERICAL SHELL AT 30 DAYS ---
[sphere_x, sphere_y, sphere_z] = sphere(20);
R_shell = 30;
surf(ax, R_shell * sphere_x, R_shell * sphere_y, R_shell * sphere_z, ...
    'FaceAlpha', 0, 'EdgeColor', [0.5, 0.5, 0.5], ...
    'EdgeAlpha', 0.8, 'LineWidth', 0.5, ...
    'HandleVisibility', 'off');

% Define marker shapes for each configuration
config_markers = containers.Map();
config_markers('1U') = '^';      % Triangle
config_markers('1.5U') = 'd';    % Diamond
config_markers('2U') = 'o';      % Circle
config_markers('3U') = 's';      % Square

% Define marker sizes
marker_size = 30;
fail_marker_size = 30;

% Plot vectors for each configuration
unique_configs = unique(config_labels);
legend_handles = [];
legend_labels = {};

for c = 1:length(unique_configs)
    config_name = unique_configs{c};
    config_idx = strcmp(config_labels, config_name);
    
    success_idx = config_idx & success_mask;
    fail_idx = config_idx & ~success_mask;
    
    marker = config_markers(config_name);
    
    % Successful trials - vectors with dots (colored by velocity)
    if any(success_idx)
        x_data = q_vectors(success_idx, 1);
        y_data = q_vectors(success_idx, 2);
        z_data = q_vectors(success_idx, 3);
        colors = color_map(success_idx, :);
        
        % Draw lines from origin (colored by velocity)
        for j = 1:length(x_data)
            plot3(ax, [0, x_data(j)], [0, y_data(j)], [0, z_data(j)], ...
                'Color', colors(j, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        end
        
        % Endpoints with shapes (colored by velocity)
        for j = 1:length(x_data)
            scatter3(ax, x_data(j), y_data(j), z_data(j), marker_size, ...
                'Marker', marker, 'MarkerFaceColor', colors(j, :), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end
        
        % Add one entry to legend for this config (black shape)
        legend_marker_size = 80;
        h = scatter3(ax, NaN, NaN, NaN, legend_marker_size, ...
            'Marker', marker, 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
            'DisplayName', config_name);
        legend_handles = [legend_handles, h];
        legend_labels{end+1} = config_name;
    end
    
    % Failed trials - Black shapes on the 30-day shell (vector still colored)
    if any(fail_idx)
        x_data = q_vectors(fail_idx, 1);
        y_data = q_vectors(fail_idx, 2);
        z_data = q_vectors(fail_idx, 3);
        colors = color_map(fail_idx, :);
        
        % Draw lines from origin (colored by velocity)
        for j = 1:length(x_data)
            plot3(ax, [0, x_data(j)], [0, y_data(j)], [0, z_data(j)], ...
                'Color', colors(j, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        end
        
        % Endpoints as BLACK shapes on the 30-day shell
        for j = 1:length(x_data)
            scatter3(ax, x_data(j), y_data(j), z_data(j), fail_marker_size, ...
                'Marker', marker, 'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                'HandleVisibility', 'off');
        end
    end
end








% ================================================================
% AXIS LABELS — NATIVE BOUNDING BOX
% ================================================================
% Axis limits and equal aspect
max_val = 35;
xlim(ax, [-max_val, max_val]);
ylim(ax, [-max_val, max_val]);
zlim(ax, [-max_val, max_val]);
axis(ax, 'equal');

% Tick marks and grid lines every 10 days
ticks = -30:10:30;

ax.XTick = ticks;
ax.YTick = ticks;
ax.ZTick = ticks;

% Capture the handles when creating the labels
hx = xlabel(ax, 'Lock Time [days]', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
hy = ylabel(ax, 'Lock Time [days]', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
zlabel(ax, 'Lock Time [days]', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% Move X label up (increase Y position)
hx.Position(1) = hx.Position(1) + 8;
hx.Position(2) = hx.Position(2) + 3;

% Move Y label right (increase X position)
hy.Position(1) = hy.Position(1) - 3;
hy.Position(2) = hy.Position(2) - 4;


% Manually rotate the X and Y labels (in degrees)
hx.Rotation = -19;  
hy.Rotation = 19; 

% Hard-lock the axis position immediately so MATLAB cannot squish the aspect ratio
set(ax, 'Position', [0.08, 0.1, 0.6, 0.8]);











% Title
title(ax, 'Initial Conditions on Attitude Sphere', 'FontSize', 23, 'FontName', 'Times New Roman');

% ================================================================
% LEGEND & COLORBAR POSITIONING
% ================================================================
% Create dedicated dummy objects ONLY for the legend
% Create dedicated dummy objects ONLY for the legend
legend_marker_size = 8;

h_1U = plot3(ax, NaN, NaN, NaN, ...
    '^', ...
    'MarkerSize', legend_marker_size, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);

h_15U = plot3(ax, NaN, NaN, NaN, ...
    'd', ...
    'MarkerSize', legend_marker_size, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);

h_2U = plot3(ax, NaN, NaN, NaN, ...
    'o', ...
    'MarkerSize', legend_marker_size, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);

h_3U = plot3(ax, NaN, NaN, NaN, ...
    's', ...
    'MarkerSize', legend_marker_size, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.5);

lgd = legend(ax, ...
    [h_1U, h_15U, h_2U, h_3U], ...
    {'1U', '1.5U', '2U', '3U', 'Failed (never locked)'}, ...
    'FontSize', 16, ...
    'FontName', 'Times New Roman');

% [left, bottom, width, height]
lgd.Position = [0.715, 0.72, 0.10, 0.15]; 
lgd.ItemTokenSize = [30, 30];

% 2. Colorbar: Positioned directly under the legend
colormap(ax, jet(64));
cb = colorbar(ax);
% [left, bottom, width, height] - top edge sits just below the legend
cb.Position = [0.75, 0.20, 0.025, 0.48]; 
cb.Label.String = 'Initial Angular Velocity (deg/s)';
cb.Label.FontName = 'Times New Roman'; 
cb.Label.FontWeight = 'bold';
cb.FontSize = 14;
cb.FontName = 'Times New Roman';
cb.Label.FontSize = 18;
caxis(ax, [0, 35]);
cb.Ticks = 0:5:35;


%% SAVE THE FIGURE AT 600 DPI
save_dir = pwd;
plot_filename = fullfile(save_dir, 'GASRATS_FullAttitudeSphere_MRP_VelocityGradient.png');
exportgraphics(gcf, plot_filename, 'Resolution', 600);
fprintf('Figure saved to: %s (600 DPI)\n', plot_filename);

% Compute total_success before printing
total_success = sum(success_mask);
fprintf('\n=== Complete ===\n');
fprintf('Total trials: %d\n', total_trials);
fprintf('Success rate (15-deg criteria): %.1f%%\n', total_success/total_trials*100);

