%% Truman Abbe | Utah State University | truman.abbe23@gmail.com
%% Modular Monte Carlo Plotting - Uses .mat file with custom settling criteria

% ==================== SETTINGS ====================
matFile = 'AttitudeSim007_Results.mat';
fprintf('Loading results from: %s\n', matFile);
S = load(matFile);

% ---------- SETTLING CRITERIA (ADJUST THESE) ----------
THETA_THRESHOLD = 15;      % degrees pointing error
W_MAG_THRESHOLD = Inf;    % deg/s angular velocity
VERIFICATION_ORBITS = 2;   % number of orbits to verify lock
ORBIT_PERIOD = 92.5;       % minutes (for 400km orbit)
PLOT_Y_MAX = 21;           % Y-axis upper limit (days)
% -----------------------------------------------------

% Extract data from MAT file first so we can check t_save
Master_Inputs = S.Master_Inputs;
w_all_matrix = S.w_all_matrix;
t_save = S.t_save;
results_theta_hist = S.results_theta_hist;
results_w_mag_hist = S.results_w_mag_hist;
total_trials = length(Master_Inputs);

% Calculate verification time in seconds
verification_time = VERIFICATION_ORBITS * ORBIT_PERIOD * 60;

% FIX: Calculate the downsampled timestep directly from the saved time array
h_save = t_save(2) - t_save(1); 
required_clean_steps = round(verification_time / h_save);

fprintf('Settling Criteria: θ<%.0f°, ω<%.2f°/s, %d orbits\n', ...
    THETA_THRESHOLD, W_MAG_THRESHOLD, VERIFICATION_ORBITS);
fprintf('Verification requires %.0f contiguous seconds (%.0f data steps at %.1fs resolution)\n', ...
    verification_time, required_clean_steps, h_save);

% ============================================================
% Re-evaluate lock times with custom criteria
% ============================================================
fprintf('Re-evaluating lock times with custom criteria...\n');
custom_lock_time = NaN(total_trials, 1);

for i = 1:total_trials
    theta_hist = results_theta_hist(i, :);
    w_mag_hist = results_w_mag_hist(i, :);
    
    violation_idx = find(theta_hist >= THETA_THRESHOLD | w_mag_hist >= W_MAG_THRESHOLD);
    
    if isempty(violation_idx)
        % Never violated - locked from start
        custom_lock_time(i) = t_save(1) / 86400;
    elseif (length(t_save) - violation_idx(end)) < required_clean_steps
        % Violation in the final verification window - NEVER achieved steady state
        custom_lock_time(i) = NaN;
    else
        % Lock occurs at the step immediately following the FINAL violation
        lock_idx = violation_idx(end) + 1;
        if lock_idx <= length(t_save)
            custom_lock_time(i) = t_save(lock_idx) / 86400;
        else
            custom_lock_time(i) = NaN;
        end
    end
end

% ============================================================
% Build configuration maps from Master_Inputs
% ============================================================
config_strings = {Master_Inputs.config}';
sat_config_map = zeros(total_trials, 1);
sat_config_map(strcmp(config_strings, '1U')) = 1;
sat_config_map(strcmp(config_strings, '1.5U')) = 2;
sat_config_map(strcmp(config_strings, '2U')) = 3;
sat_config_map(strcmp(config_strings, '3U')) = 4;

rod_counts = [Master_Inputs.rod_config]' * 2;
rod_config_map = zeros(total_trials, 1);
rod_config_map(rod_counts == 2) = 1;
rod_config_map(rod_counts == 4) = 2;

w_mags = vecnorm(w_all_matrix, 2, 2);
results_lock_time_clean = custom_lock_time;

% ============================================================
% GENERATE PLOT - IDENTICAL TO ORIGINAL
% ============================================================
fprintf('Generating maximized 6.5in plot...\n');

fSizeLabel  = 11;
fSizeTick   = 10;
fSizeLegend = 10;
fSizeAnnot  = 10;
fName       = 'Times New Roman';
figWidth  = 6.5; 
figHeight = 3.5; 

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, figWidth, figHeight], ...
             'Visible', 'on', 'Renderer', 'painters');
ax = axes(fig, 'Units', 'normalized', 'Position', [0.08 0.13 0.90 0.82], ...
          'FontName', fName, 'FontSize', fSizeTick, 'TickDir', 'in'); 
hold(ax, 'on'); grid(ax, 'on');

xline(ax, norm([5,5,5]), '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(ax, 8.2, 21.0, 'Nominal Limit', ...
    'Rotation', 90, ...
    'FontSize', fSizeAnnot, ...
    'FontName', fName, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');

yline(ax, 7, '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(ax, 1.6, 7.45, 'Mission Requirement', ...
    'Rotation', 90, ...
    'FontSize', fSizeAnnot, ...
    'FontName', fName, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');

% 14-Day Line
yline(ax, 14, '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');
text(ax, 4, 14.45, sprintf('Secondary\nThreshold'), ...
    'Rotation', 90, ...
    'FontSize', fSizeAnnot, ...
    'FontName', fName, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');

markers = {'^', 'd', 'o', 's'}; 

for c = 1:4
    for r = 1:2
        c_color = 'r'; if r == 2, c_color = 'b'; end
        c_marker = markers{c};
        idx = find(sat_config_map == c & rod_config_map == r);
        if isempty(idx), continue; end
        
        % Get the actual data for this group
        lock_times = results_lock_time_clean(idx);
        w_mags_group = w_mags(idx);
        
        % ============================================================
        % PLOT TRIALS THAT ACHIEVE STEADY STATE WITHIN PLOT BOUNDS
        % ============================================================
        valid_mask = ~isnan(lock_times);
        within_bounds_mask = valid_mask & (lock_times <= PLOT_Y_MAX);
        
        w_success = w_mags_group(within_bounds_mask);
        t_success = lock_times(within_bounds_mask);
        
        if ~isempty(w_success)
            scatter(ax, w_success, t_success, 20, ...
                'MarkerFaceColor', c_color, 'MarkerEdgeColor', 'k', ...
                'Marker', c_marker, 'LineWidth', 0.3, 'HandleVisibility', 'off');
        end
        
        % ============================================================
        % PLOT TRIALS THAT DON'T ACHIEVE STEADY STATE WITHIN BOUNDS
        % ============================================================
        % Failures: NEVER achieved OR achieved after PLOT_Y_MAX
        fail_mask = isnan(lock_times) | (valid_mask & (lock_times > PLOT_Y_MAX));
        w_fail = w_mags_group(fail_mask);
        
        if ~isempty(w_fail)
            scatter(ax, w_fail, PLOT_Y_MAX * ones(length(w_fail), 1), 25, ...
                'Marker', 'x', 'MarkerEdgeColor', c_color, 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
        end
    end
end

set(ax, 'Box', 'on', 'LineWidth', 0.8); 
xl = xlabel(ax, 'Initial Angular Velocity Magnitude [deg/s]', 'FontSize', fSizeLabel, 'FontName', fName);
xl.Units = 'normalized';
xl.Position(2) = -0.08;

yl = ylabel(ax, 'Time to Steady-State Attitude [days]', 'FontSize', fSizeLabel, 'FontName', fName);
yl.Units = 'normalized';
yl.Position(1) = -0.04; 

xlim(ax, [0 35]); 
ylim(ax, [0 PLOT_Y_MAX + 0.5]); 
yticks(ax, 0:3:PLOT_Y_MAX); 

h = zeros(5, 1);
h(1) = scatter(ax, NaN, NaN, 30, 'k^', 'MarkerFaceColor', 'k', 'DisplayName', '1U');
h(2) = scatter(ax, NaN, NaN, 30, 'kd', 'MarkerFaceColor', 'k', 'DisplayName', '1.5U');
h(3) = scatter(ax, NaN, NaN, 30, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '2U');
h(4) = scatter(ax, NaN, NaN, 30, 'ks', 'MarkerFaceColor', 'k', 'DisplayName', '3U');
%h(5) = scatter(ax, NaN, NaN, 30, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '2 Rods');
h(5) = scatter(ax, NaN, NaN, 30, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', '4 Rods');
lgd = legend(ax, h, 'Location', 'northeast', 'Box', 'on', ...
             'FontName', fName, 'FontSize', fSizeLegend, 'NumColumns', 3);
drawnow;

% Save with criteria info in filename
targetPath = fullfile(pwd, sprintf('MC_Visual_theta%.0f_omega%.2f_orbits%d.png', ...
    THETA_THRESHOLD, W_MAG_THRESHOLD, VERIFICATION_ORBITS));
exportgraphics(fig, targetPath, 'Resolution', 600, 'BackgroundColor', 'none');
fprintf('Plot saved to: %s\n', targetPath);

% ============================================================
% Print statistics
% ============================================================
num_within_bounds = sum(~isnan(custom_lock_time) & custom_lock_time <= PLOT_Y_MAX);
num_beyond_bounds = sum(~isnan(custom_lock_time) & custom_lock_time > PLOT_Y_MAX);
num_never_steady = sum(isnan(custom_lock_time));
num_fail_total = num_beyond_bounds + num_never_steady;

fprintf('\n========== RESULTS ==========\n');
fprintf('Achieved steady state within %.0f days: %d / %d (%.1f%%)\n', ...
    PLOT_Y_MAX, num_within_bounds, total_trials, 100*num_within_bounds/total_trials);
fprintf('Achieved steady state after %.0f days: %d / %d (%.1f%%)\n', ...
    PLOT_Y_MAX, num_beyond_bounds, total_trials, 100*num_beyond_bounds/total_trials);
fprintf('NEVER achieved steady state:     %d / %d (%.1f%%)\n', ...
    num_never_steady, total_trials, 100*num_never_steady/total_trials);
fprintf('Total failures (plotted as X):   %d / %d (%.1f%%)\n', ...
    num_fail_total, total_trials, 100*num_fail_total/total_trials);
fprintf('=============================\n');