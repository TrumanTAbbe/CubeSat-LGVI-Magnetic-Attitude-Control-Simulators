%% Truman Abbe | Utah State University | truman.abbe23@gmail.com
%% Post-Processing: Multiple Omega Threshold Evaluation for 4-Rod Trials
%% Loads modular .mat data and recomputes lock times for various criteria
%% Uses downsampled theta_hist and w_mag_hist

clear; clc; close all;

%% 1. Load the .mat file
% Adjust filename as needed
matFile = 'AttitudeSim007_Results.mat'; % Your modular data file
if ~exist(matFile, 'file')
    error('File "%s" not found. Please provide the correct path.', matFile);
end
fprintf('Loading %s...\n', matFile);
load(matFile); % loads: Master_Inputs, t_save, results_w_mag_hist, results_theta_hist, results_ortho_hist

%% 2. Simulation constants (from original code)
h = 0.05;                 % timestep [s]
t_end = 30 * 86400;       % 30 days [s]
required_clean_time = 10800; % 2 orbits (10800 s) verification window

% Define bin edges based on norm([5,5,5])
baseVec = [5,5,5];
baseNorm = norm(baseVec); % = 8.660254037844386
binEdges = [0, baseNorm, 2*baseNorm, 3*baseNorm, 4*baseNorm];
numBins = length(binEdges) - 1;
binLabels = arrayfun(@(b) sprintf('[%.2f, %.2f)', binEdges(b), binEdges(b+1)), ...
                     1:numBins, 'UniformOutput', false);

% Satellite configurations
configNames = {'1U', '1.5U', '2U', '3U'};
numConfigs = length(configNames);

totalTrials = length(Master_Inputs);
fprintf('Total trials loaded: %d\n', totalTrials);

%% 3. Identify 4-rod trials only
% Master_Inputs(i).rod_config = 1 or 2 (rods per axis)
% Total rods = rod_config * 2
is4Rod = false(totalTrials, 1);
for i = 1:totalTrials
    if Master_Inputs(i).rod_config == 2
        is4Rod(i) = true;
    end
end
fourRodTrials = find(is4Rod);
numFourRod = length(fourRodTrials);
fprintf('4-rod trials: %d out of %d\n', numFourRod, totalTrials);

% Extract data for 4-rod trials only
configIdx = zeros(numFourRod, 1);
wMag_init = zeros(numFourRod, 1);
for j = 1:numFourRod
    i = fourRodTrials(j);
    configStr = Master_Inputs(i).config;
    configIdx(j) = find(strcmp(configStr, configNames));
    wMag_init(j) = norm(Master_Inputs(i).w0_deg);
end

%% 4. Determine velocity bin for each 4-rod trial
% Bins:
%   Bin 1 :   0 <= w <  8.660254038
%   Bin 2 : 8.660254038 <= w < 17.320508076
%   Bin 3 : 17.320508076 <= w < 25.980762114
%   Bin 4 : 25.980762114 <= w <= 35

baseNorm = norm([5 5 5]);

binEdges = [0 ...
            baseNorm ...
            2*baseNorm ...
            3*baseNorm ...
            35];

binLabels = { ...
    sprintf('[%.2f, %.2f)',binEdges(1),binEdges(2)), ...
    sprintf('[%.2f, %.2f)',binEdges(2),binEdges(3)), ...
    sprintf('[%.2f, %.2f)',binEdges(3),binEdges(4)), ...
    sprintf('[%.2f, %.2f]',binEdges(4),binEdges(5))};

numBins = 4;

binIdx = zeros(numFourRod,1);

tol = 1e-10;

for j = 1:numFourRod

    w = wMag_init(j);

    if w < binEdges(2) - tol
        binIdx(j) = 1;

    elseif w < binEdges(3) - tol
        binIdx(j) = 2;

    elseif w < binEdges(4) - tol
        binIdx(j) = 3;

    elseif w <= binEdges(5) + tol
        binIdx(j) = 4;

    else
        error('Initial angular velocity %.12f exceeds 35 deg/s.',w);
    end

end

% Diagnostic printout (should match published paper)
fprintf('\nTrials per configuration/bin:\n');
for c = 1:numConfigs
    fprintf('%s : ',configNames{c});
    for b = 1:numBins
        fprintf('%3d ',sum(configIdx==c & binIdx==b));
    end
    fprintf('\n');
end

%% 5. Define omega thresholds to test
%omega_thresholds = [0.5, 0.75, 1.0, 1.25, 1.5, Inf]; % deg/s
%threshold_labels = {'<0.5', '<0.75', '<1.0', '<1.25', '<1.5', 'theta only'};
%numThresholds = length(omega_thresholds);

%% 5. Define theta thresholds to test (pointing error only)
theta_thresholds = [20, 15, 10, 5]; % deg
threshold_labels = {'<20 deg', '<15 deg', '<10 deg', '<5 deg'};
numThresholds = length(theta_thresholds);

% Initialize storage for lock times and success rates
lock_times_all = NaN(numFourRod, numThresholds);
success_rates = zeros(numBins, numConfigs, numThresholds);
counts = zeros(numBins, numConfigs);

fprintf('Computing lock times for %d different thresholds...\n', numThresholds);

%% 6. Compute lock times for all thresholds
for j = 1:numFourRod
    i = fourRodTrials(j);
    theta_hist = results_theta_hist(i, :);
    w_mag_hist = results_w_mag_hist(i, :);
    t_hist = t_save;
%{    
    for t = 1:numThresholds
        omega_thresh = omega_thresholds(t);
        
        % Find violations based on this threshold
        if isinf(omega_thresh)
            % theta only (no omega requirement)
            violation_idx = find(theta_hist >= 15);
        else
            % theta AND omega threshold
            violation_idx = find(theta_hist >= 15 | w_mag_hist >= omega_thresh);
        end
        
        if isempty(violation_idx)
            lock_times_all(j, t) = 0;
        elseif (t_end - t_hist(violation_idx(end))) < required_clean_time
            lock_times_all(j, t) = NaN;
        else
            if violation_idx(end) < length(t_hist)
                lock_idx = violation_idx(end) + 1;
                lock_times_all(j, t) = t_hist(lock_idx) / 86400;
            else
                lock_times_all(j, t) = NaN;
            end
        end
    end
%}
    for t = 1:numThresholds
        theta_thresh = theta_thresholds(t);
        
        % Find violations: pointing error exceeds threshold
        violation_idx = find(theta_hist >= theta_thresh);
        
        if isempty(violation_idx)
            lock_times_all(j, t) = 0;
        elseif (t_end - t_hist(violation_idx(end))) < required_clean_time
            lock_times_all(j, t) = NaN;
        else
            if violation_idx(end) < length(t_hist)
                lock_idx = violation_idx(end) + 1;
                lock_times_all(j, t) = t_hist(lock_idx) / 86400;
            else
                lock_times_all(j, t) = NaN;
            end
        end
    end
    
end

%% 7. Compute success rates for each bin, config, and threshold
threshold_days = 21;

for b = 1:numBins
    for c = 1:numConfigs
        idx = find(binIdx == b & configIdx == c);
        if isempty(idx)
            continue;
        end
        counts(b,c) = length(idx);
        
        for t = 1:numThresholds
            lock_days = lock_times_all(idx, t);
            success = ~isnan(lock_days) & lock_days < threshold_days;
            success_rates(b,c,t) = sum(success) / length(idx) * 100;
        end
    end
end

%% 8. Display results - One table per threshold
fprintf('\n================================================================================\n');
%fprintf('4-ROD TRIALS ONLY: Success Rates (Lock < %d days) by Omega Threshold\n', threshold_days);
fprintf('4-ROD TRIALS ONLY: Success Rates (Lock < %d days) by Pointing Error Threshold\n', threshold_days);
fprintf('Total 4-rod trials: %d\n', numFourRod);
fprintf('================================================================================\n');

% Determine column width for formatting (MOVED OUTSIDE THE LOOP)
colWidth = 14;

for t = 1:numThresholds
    %fprintf('\n=== THRESHOLD: theta < 15 deg AND omega %s deg/s ===\n', threshold_labels{t});
    fprintf('\n=== THRESHOLD: theta < %d deg (pointing error only) ===\n', theta_thresholds(t));
    %if isinf(omega_thresholds(t))
    %    fprintf('=== (theta only, no omega requirement) ===\n');
    %end
    
    % Compute total successful trials per config for this threshold
    fprintf('Successful trials: ');
    for c = 1:numConfigs
        % Count successes for this config and threshold
        succ_count = 0;
        for b = 1:numBins
            if counts(b,c) > 0
                rate = success_rates(b,c,t);
                succ_count = succ_count + round(rate/100 * counts(b,c));
            end
        end
        fprintf('%s:%d ', configNames{c}, succ_count);
    end
    fprintf('\n');
    
    fprintf('%-15s', 'Velocity Bin');
    for c = 1:numConfigs
        fprintf('%-*s', colWidth+4, configNames{c});
    end
    fprintf('\n');
    
    fprintf('%-15s', '------------');
    for c = 1:numConfigs
        fprintf('%-*s', colWidth+4, '------------');
    end
    fprintf('\n');
    
    for b = 1:numBins
        fprintf('%-15s', binLabels{b});
        for c = 1:numConfigs
            if counts(b,c) > 0
                rate = success_rates(b,c,t);
                n = counts(b,c);
                succ = round(rate/100 * n);
                fprintf('%-*s', colWidth+4, sprintf('%.2f%% (%d/%d)', rate, succ, n));
            else
                fprintf('%-*s', colWidth+4, 'N/A');
            end
        end
        fprintf('\n');
    end
end

%% 9. Overall summary (all bins combined) for each threshold
fprintf('\n================================================================================\n');
%fprintf('OVERALL SUCCESS RATES (All Bins Combined) by Omega Threshold\n');
fprintf('OVERALL SUCCESS RATES (All Bins Combined) by Pointing Error Threshold\n');
fprintf('================================================================================\n');
fprintf('%-10s', 'Config');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, threshold_labels{t});
end
fprintf('%-*s\n', colWidth+2, 'Total');

fprintf('%-10s', '------');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, '----------');
end
fprintf('%-*s\n', colWidth+2, '-----');

for c = 1:numConfigs
    total = 0;
    for b = 1:numBins
        total = total + counts(b,c);
    end
    if total == 0
        continue;
    end
    
    fprintf('%-10s', configNames{c});
    for t = 1:numThresholds
        tot_succ = 0;
        for b = 1:numBins
            if counts(b,c) > 0
                tot_succ = tot_succ + round(success_rates(b,c,t)/100 * counts(b,c));
            end
        end
        fprintf('%-*s', colWidth+2, sprintf('%.2f%%', tot_succ/total*100));
    end
    fprintf('%-*d\n', colWidth+2, total);
end

%% 10. Mean and Median Lock Times for each threshold
fprintf('\n================================================================================\n');
%fprintf('MEAN LOCK TIMES (hrs) by Omega Threshold (Successful Trials Only)\n');
fprintf('MEAN LOCK TIMES (hrs) by Pointing Error Threshold (Successful Trials Only)\n');
fprintf('================================================================================\n');
fprintf('%-10s', 'Config');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, threshold_labels{t});
end
fprintf('\n');

fprintf('%-10s', '------');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, '----------');
end
fprintf('\n');

for c = 1:numConfigs
    idx = find(configIdx == c);
    if isempty(idx)
        continue;
    end
    
    fprintf('%-10s', configNames{c});
    for t = 1:numThresholds
        lock_times = lock_times_all(idx, t);
        success = ~isnan(lock_times);
        if sum(success) > 0
            lock_hrs = lock_times(success) * 24;
            fprintf('%-*s', colWidth+2, sprintf('%.2f', mean(lock_hrs)));
        else
            fprintf('%-*s', colWidth+2, 'N/A');
        end
    end
    fprintf('\n');
end

%% 11. Median Lock Times for each threshold
fprintf('\n================================================================================\n');
%fprintf('MEDIAN LOCK TIMES (hrs) by Omega Threshold (Successful Trials Only)\n');
fprintf('MEDIAN LOCK TIMES (hrs) by Pointing Error Threshold (Successful Trials Only)\n');
fprintf('================================================================================\n');
fprintf('%-10s', 'Config');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, threshold_labels{t});
end
fprintf('\n');

fprintf('%-10s', '------');
for t = 1:numThresholds
    fprintf('%-*s', colWidth+2, '----------');
end
fprintf('\n');

for c = 1:numConfigs
    idx = find(configIdx == c);
    if isempty(idx)
        continue;
    end
    
    fprintf('%-10s', configNames{c});
    for t = 1:numThresholds
        lock_times = lock_times_all(idx, t);
        success = ~isnan(lock_times);
        if sum(success) > 0
            lock_hrs = lock_times(success) * 24;
            fprintf('%-*s', colWidth+2, sprintf('%.2f', median(lock_hrs)));
        else
            fprintf('%-*s', colWidth+2, 'N/A');
        end
    end
    fprintf('\n');
end

%% 12. Save results (optional)
saveResults = false;
if saveResults
    resultsFile = ['LockTimeAnalysis_MultiThreshold_', datestr(now, 'yyyy-mm-dd_HHMM'), '.mat'];
    save(resultsFile, 'lock_times_all', 'success_rates', 'counts', ...
         'omega_thresholds', 'threshold_labels', ...
         'binIdx', 'configIdx', 'wMag_init', 'fourRodTrials');
    fprintf('\nResults saved to: %s\n', resultsFile);
end

fprintf('\nProcessing complete.\n');