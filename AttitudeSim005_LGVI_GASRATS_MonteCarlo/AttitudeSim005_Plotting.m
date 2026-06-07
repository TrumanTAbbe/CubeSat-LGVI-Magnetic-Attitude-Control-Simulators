%% Truman Abbe | Utah State University | truman.abbe23@gmail.com
%% Plotting For 1,001-Trial Monte Carlo Study

matFile = 'AttitudeSim005_Results.mat'; 
reportFile = 'AttitudeSim005_Text_Report.txt';
fprintf('Loading results...\n');
S = load(matFile);
opts = detectImportOptions(reportFile);
opts.VariableNamingRule = 'preserve'; 
T_report = readtable(reportFile, opts); 

config_strings = T_report.Sat_Config;
sat_config_map = zeros(height(T_report), 1);
sat_config_map(strcmp(config_strings, '1U')) = 1;
sat_config_map(strcmp(config_strings, '1.5U')) = 2;
sat_config_map(strcmp(config_strings, '2U')) = 3;
sat_config_map(strcmp(config_strings, '3U')) = 4;

rod_counts = T_report.Total_Rods;
rod_config_map = zeros(height(T_report), 1);
rod_config_map(rod_counts == 2) = 1;
rod_config_map(rod_counts == 4) = 2;

w_mags = T_report.("W_Mag(deg/s)"); 
results_lock_time_clean = T_report.Lock_Day;
failed_mask = strcmp(T_report.Reached_Steady_State, 'No');
results_lock_time_clean(failed_mask) = NaN;

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

xline(ax, 8.5, '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');

text(ax, 8.2, 21.0, 'Nominal Limit', ...
    'Rotation', 90, ...
    'FontSize', fSizeAnnot, ...
    'FontName', fName, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'right');

yline(ax, 7, '--k', 'LineWidth', 1.4, 'HandleVisibility', 'off');

text(ax, 1.6, 7.45, '7-Day Requirement', ...
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
        
        success_idx = idx(~isnan(results_lock_time_clean(idx)));
        if ~isempty(success_idx)
            scatter(ax, w_mags(success_idx), results_lock_time_clean(success_idx), 20, ...
                'MarkerFaceColor', c_color, 'MarkerEdgeColor', 'k', ...
                'Marker', c_marker, 'LineWidth', 0.3, 'HandleVisibility', 'off');
        end
        
        fail_idx = idx(isnan(results_lock_time_clean(idx)));
        if ~isempty(fail_idx)
            scatter(ax, w_mags(fail_idx), 21 * ones(length(fail_idx), 1), 25, ...
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
ylim(ax, [0 21.5]); 
yticks(ax, 0:3:21); 

h = zeros(6, 1);
h(1) = scatter(ax, NaN, NaN, 30, 'k^', 'MarkerFaceColor', 'k', 'DisplayName', '1U');
h(2) = scatter(ax, NaN, NaN, 30, 'kd', 'MarkerFaceColor', 'k', 'DisplayName', '1.5U');
h(3) = scatter(ax, NaN, NaN, 30, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', '2U');
h(4) = scatter(ax, NaN, NaN, 30, 'ks', 'MarkerFaceColor', 'k', 'DisplayName', '3U');
h(5) = scatter(ax, NaN, NaN, 30, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '2 Rods');
h(6) = scatter(ax, NaN, NaN, 30, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', '4 Rods');
lgd = legend(ax, h, 'Location', 'northeast', 'Box', 'on', ...
             'FontName', fName, 'FontSize', fSizeLegend, 'NumColumns', 3);

drawnow;

targetPath = fullfile(pwd, '1,001-Trial_Monte_Carlo_Visual.png');
exportgraphics(fig, targetPath, 'Resolution', 600, 'BackgroundColor', 'none');