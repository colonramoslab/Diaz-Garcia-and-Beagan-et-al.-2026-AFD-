clc
clear
close all
load("C:\Users\dgmdi\DiazGarcia_Beagan_2026_Paper_GithubRepo\Runs_and_polar_plots_Steep_gradient/merged_X.mat")
load("C:\Users\dgmdi\DiazGarcia_Beagan_2026_Paper_GithubRepo\Runs_and_polar_plots_Steep_gradient/merged_Y.mat")
X = merged_X;
Y = merged_Y;
%% =========================
% USER SETTINGS
% =========================
sampling_rate = 2;
sampling_rate_after_downsample = 1;
analysis_minutes = 15;
start_time_point = 1;
end_time_point   = sampling_rate * 60 * analysis_minutes;
pix_range  = [900 1800];
temp_range = [15 25];
temp_bins = [
    20 22
    22 24
];
num_runs_to_display = 40;
trajectory_color = [1 0 0];
ylimits_tracks = [-100 100];
xlimits_tracks = [-60 120];
polar_bin_deg = 10;
polar_radius_lim_mm = 6;
pix_to_mm = 100/900;
%% =========================
% EXTRACT TIME WINDOW
% =========================
x_window = X(start_time_point:end_time_point, :);
y_window = Y(start_time_point:end_time_point, :);
nBins = size(temp_bins, 1);
collected_run_data = cell(nBins, 1);
%% =========================
% PLOT 1: RUN TRAJECTORIES
% =========================
figure;
tiledlayout(1, nBins, 'TileSpacing', 'compact', 'Padding', 'compact');
for b = 1:nBins
    T_low  = temp_bins(b,1);
    T_high = temp_bins(b,2);
    
    x_low  = temp_to_pix(T_low,  pix_range, temp_range);
    x_high = temp_to_pix(T_high, pix_range, temp_range);
    
    x_bin = x_window;
    y_bin = y_window;
    
    mask = x_bin >= x_low & x_bin <= x_high;
    x_bin(~mask) = NaN;
    y_bin(~mask) = NaN;
    
    [~,~,~,~,coordinates_to_plot] = visualize_tracks_runs_and_turns_5_Tc25( ...
        x_bin, y_bin, 'data', 'k', 0);
        
    true_runs = extract_run_segments_from_coordinates_3( ...
        coordinates_to_plot, sampling_rate_after_downsample);
        
    nexttile;
    plot_true_runs_from_origin_topN_2( ...
        true_runs, trajectory_color, ...
        sprintf('Steep: %.1f-%.1f °C', T_low, T_high), ...
        num_runs_to_display, ylimits_tracks, xlimits_tracks);
        
    ax = gca;
    set(ax, 'FontSize', 18)
    
    x_ticks_mm = [-5 0 5 10];
    y_ticks_mm = [-10 -5 0 5 10];
    ax.XTick = x_ticks_mm / pix_to_mm;
    ax.YTick = y_ticks_mm / pix_to_mm;
    ax.XTickLabel = string(x_ticks_mm);
    ax.YTickLabel = string(y_ticks_mm);
    
    xlabel('X Position (mm)')
    ylabel('Y Position (mm)')
end
%% =========================
% PLOT 2: POLAR DISTANCE PLOTS
% =========================
figure;
for b = 1:nBins
    ax_tmp = subplot(1, nBins, b);
    pos = get(ax_tmp, 'Position');
    delete(ax_tmp);
    pax = polaraxes('Position', pos);
    hold(pax, 'on');
    
    T_low  = temp_bins(b,1);
    T_high = temp_bins(b,2);
    
    x_low  = temp_to_pix(T_low,  pix_range, temp_range);
    x_high = temp_to_pix(T_high, pix_range, temp_range);
    
    x_bin = x_window;
    y_bin = y_window;
    
    mask = x_bin >= x_low & x_bin <= x_high;
    x_bin(~mask) = NaN;
    y_bin(~mask) = NaN;
    
    [~,~,~,~,coordinates_to_plot] = visualize_tracks_runs_and_turns_5_Tc25( ...
        x_bin, y_bin, 'data', 'k', 0);
    
    [run_headings, run_distances, ~, run_worm_ids] = extract_runs_from_coordinates_3( ...
        coordinates_to_plot, sampling_rate_after_downsample);
    
    run_distances = run_distances * pix_to_mm;
    
    collected_run_data{b}.headings  = run_headings;
    collected_run_data{b}.distances = run_distances;
    collected_run_data{b}.worm_ids  = run_worm_ids;
    collected_run_data{b}.T_low     = T_low;
    collected_run_data{b}.T_high    = T_high;
    
    if isempty(run_headings) || isempty(run_distances)
        title(pax, sprintf('%.1f-%.1f °C\n(no runs)', T_low, T_high), ...
            'FontWeight', 'normal', 'FontSize', 18);
        continue
    end
    
    axes(pax);
    polar_distance_plot_color_radius_lim_3( ...
        run_headings, ...
        run_distances, ...
        polar_bin_deg, ...
        '', ...
        trajectory_color, ...
        0, polar_radius_lim_mm);
        
    pax.ThetaZeroLocation = 'left';
    pax.ThetaDir          = 'counterclockwise';
    pax.ThetaLim          = [-180 180];
    pax.ThetaTick         = [-180 -90 0 90];
    pax.FontSize          = 18;
    pax.GridColor         = 'k';
    pax.GridAlpha         = 1;
    pax.RTick             = [0 3 6];
    
    title(pax, sprintf('%.1f-%.1f °C', T_low, T_high), ...
        'FontWeight', 'normal', 'FontSize', 18);
end
%% ============================================
% STATISTICAL ANALYSIS (POOLED RUNS)
% ============================================
tbl_reg  = [];
tbl_len  = [];
tbl_type = [];
tbl_worm = [];
tbl_steep = [];

for b = 1:nBins
    if isempty(collected_run_data{b}), continue; end
    
    h_raw = collected_run_data{b}.headings;
    d_raw = collected_run_data{b}.distances;
    w_raw = collected_run_data{b}.worm_ids;
    
    % FILTER: Drop runs oriented to the cooling side
    keep_mask = cosd(h_raw) <= 0;
    
    h = h_raw(keep_mask);
    d = d_raw(keep_mask);
    w = w_raw(keep_mask);
    
    % COMPUTE STEEPNESS
    steepness = abs(sind(h));
    is_steep = steepness > 0.90;
    
    % Populate pooled tables
    for i = 1:length(d)
        if isnan(d(i)), continue; end
        
        tbl_reg  = [tbl_reg; b];
        tbl_len  = [tbl_len; d(i)];
        tbl_worm = [tbl_worm; w(i)];
        tbl_steep = [tbl_steep; steepness(i)];
        
        if is_steep(i)
            tbl_type = [tbl_type; "Steep"];
        else
            tbl_type = [tbl_type; "Shallow"];
        end
    end
end

% Step 5: Create a publication-quality figure
figure('Name', 'Pooled Run Lengths (Warming Side Only)', 'Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

% Combine Region and Type for grouping on the X-axis
x_grp = string(tbl_reg) + " " + tbl_type;
x_grp = categorical(x_grp, ["1 Steep", "1 Shallow", "2 Steep", "2 Shallow"]);

% Use Boxchart for the foundation
bchart = boxchart(x_grp, tbl_len, 'BoxFaceColor', [0.7 0.7 0.7], 'MarkerStyle', 'none');

% Overlay each individual run as a scatter point
scatter(x_grp, tbl_len, 15, 'k', 'filled', 'XJitter', 'density', 'XJitterWidth', 0.4);

ylabel('Run Length (mm)');
title({'Pooled Run Lengths: Steep vs Shallow by Region','[Cooling side runs excluded]'}, 'FontWeight', 'normal');
set(gca, 'FontSize', 18);
box off;

% Step 6 & 7: Statistical testing (Unpaired) & Command Window Report
fprintf('\n\n%% ============================================\n');
fprintf('%% COMMAND WINDOW STATISTICAL REPORT (Pooled Runs)\n');
fprintf('%% ============================================\n\n');

% Extract arrays for easy statistical comparison
reg1_steep   = tbl_len(tbl_reg == 1 & tbl_type == "Steep");
reg1_shallow = tbl_len(tbl_reg == 1 & tbl_type == "Shallow");
reg2_steep   = tbl_len(tbl_reg == 2 & tbl_type == "Steep");
reg2_shallow = tbl_len(tbl_reg == 2 & tbl_type == "Shallow");

% 1. Region 1: Steep vs Shallow
do_unpaired_stats('Region 1: Steep vs Shallow', reg1_steep, reg1_shallow);
% 2. Region 2: Steep vs Shallow
do_unpaired_stats('Region 2: Steep vs Shallow', reg2_steep, reg2_shallow);
% 3. Steep: Region 1 vs Region 2
do_unpaired_stats('Steep: Region 1 vs Region 2', reg1_steep, reg2_steep);
% 4. Shallow: Region 1 vs Region 2
do_unpaired_stats('Shallow: Region 1 vs Region 2', reg1_shallow, reg2_shallow);

%% ============================================
% MIXED-EFFECTS MODEL ON RUN-LEVEL DATA
% ============================================
if ~isempty(tbl_worm)

    tbl = table( ...
        categorical(tbl_worm), ...
        categorical(tbl_reg), ...
        tbl_len, ...
        tbl_steep, ...
        'VariableNames', ...
        {'AnimalID','Region','RunLength','Steepness'});

    %% ============================================
    % DIAGNOSTICS
    %% ============================================

    fprintf('\n==============================\n');
    fprintf('ANIMAL DIAGNOSTICS\n');
    fprintf('==============================\n');
    
    fprintf('Total runs: %d\n',height(tbl));
    fprintf('Unique animals: %d\n',numel(unique(tbl.AnimalID)));
    
    % ---- Runs per animal ----
    [Gid,~,idx] = unique(tbl.AnimalID);
    
    nRuns = accumarray(idx,1);
    meanRun = accumarray(idx,tbl.RunLength,[],@mean);
    
    disp(table(Gid,nRuns));
    
    fprintf('\nVariance of animal mean run lengths = %.4f\n',var(meanRun));
    
    % ---- Mean run length per animal ----
    figure('Name','Mean Run Length per Animal');
    scatter(1:numel(meanRun),meanRun,50,'filled');
    xlabel('Animal');
    ylabel('Mean Run Length');
    
    % ---- Histogram of animal means ----
    figure('Name','Distribution of Animal Means');
    histogram(meanRun);
    xlabel('Mean Run Length');
    ylabel('Number of Animals');
    
    % ---- Boxplot ----
    figure('Name','Run Length by Animal');
    boxplot(tbl.RunLength,idx);
    xlabel('Animal');
    ylabel('Run Length');

    %% ============================================
    % FIT MIXED MODEL
    %% ============================================
    lme = fitlme(tbl, 'RunLength ~ Region*Steepness + (1|AnimalID)');
    
    fprintf('\n==============================\n');
    fprintf('MIXED EFFECTS MODEL RESULTS\n');
    fprintf('==============================\n');
    disp(lme);
    
    fprintf('\n--- ANOVA Table ---\n');
    disp(anova(lme));

    fprintf('\n==============================\n');
    fprintf('COMPARISON TO LINEAR MODEL\n');
    fprintf('==============================\n');
    
    lm = fitlm(tbl,'RunLength ~ Region*Steepness');
    
    disp(lm)
    
    fprintf('\nLME coefficients:\n');
    disp(lme.Coefficients(:,{'Estimate','pValue'}))
    
    fprintf('\nLM coefficients:\n');
    disp(lm.Coefficients(:,{'Estimate','pValue'}))

    fprintf('==============================\n\n');
end

%% =========================
% HELPER FUNCTIONS
% =========================
function x_pix = temp_to_pix(T, pix_range, temp_range)
    x_pix = pix_range(1) + (T - temp_range(1)) * ...
        (pix_range(2) - pix_range(1)) / (temp_range(2) - temp_range(1));
end

function [all_heading_angles, all_distances, all_mean_headings_per_worm, all_mean_distances_per_worm, coordinates_to_plot] = ...
    visualize_tracks_runs_and_turns_5_Tc25(x_coord, y_coord, plot_title, ttx_color, plot_figures)
    windowWidth = 7;           
    downsampling_order = 2;
    polynomial_order = 3;      
    reversal_angle_threshold = 120;
    min_reversal_ruration = 2;
    window_size = 5;
    if mod(windowWidth, 2) == 0, windowWidth = windowWidth + 1; end
    if windowWidth < polynomial_order + 1
        windowWidth = polynomial_order + 1;
        if mod(windowWidth, 2) == 0, windowWidth = windowWidth + 1; end
    end
    all_heading_angles = []; all_distances = [];
    [T, M] = size(x_coord);
    all_mean_headings_per_worm = nan(M, 1);
    all_mean_distances_per_worm = nan(M, 1);
    emptyEntry = struct('worm_id', [], 'x_smooth', {{}}, 'y_smooth', {{}}, ...
        'rev_index', {{}}, 'heading_angle', {{}}, 'seg_start_idx', [], 'seg_end_idx', [], 'x0_worm', []);
    coordinates_to_plot = repmat(emptyEntry, 1, M);
    for col = 1:M
        xcol = x_coord(:, col); ycol = y_coord(:, col);
        good = ~isnan(xcol) & ~isnan(ycol);
        if ~any(good)
            coordinates_to_plot(col).worm_id = col; continue
        end
        first_valid_idx = find(good, 1, 'first');
        x0_worm = xcol(first_valid_idx);
        d = diff([false; good; false]);
        starts = find(d == 1); ends_  = find(d == -1) - 1;
        worm_heading_concat = []; worm_distance_concat = []; seg_counter = 0;
        for seg = 1:numel(starts)
            xs = xcol(starts(seg):ends_(seg)); ys = ycol(starts(seg):ends_(seg));
            if numel(xs) < max(windowWidth, polynomial_order + 1), continue; end
            smoothedX = downsample(sgolayfilt(xs, polynomial_order, windowWidth), downsampling_order);
            smoothedY = downsample(sgolayfilt(ys, polynomial_order, windowWidth), downsampling_order);
            reversal_indices = findWormReversals_2(smoothedX, smoothedY, ...
                'WindowSize', window_size, 'AngleThreshold', reversal_angle_threshold, ...
                'MinReversalDuration', min_reversal_ruration, 'PlotResults', false);
            if plot_figures == 1
                hold on; plot(smoothedX - x0_worm, smoothedY, 'LineStyle', '-', 'Color', ttx_color, 'LineWidth', 1)
                scatter(smoothedX - x0_worm, smoothedY, 4, ttx_color, 'filled')
                if ~isempty(reversal_indices)
                    plot(smoothedX(reversal_indices) - x0_worm, smoothedY(reversal_indices), 'ro', 'MarkerSize', 3, 'LineWidth', 2, 'MarkerFaceColor', 'r');
                end
                xlim([-1000 400]); xlabel('X-coordinates'); yticks([]); set(gca, 'fontsize', 16); set(gca, 'YColor', 'none');
            end
            [heading_worm, distance_worm] = calculateRelativeHeadingandDistance(smoothedX, smoothedY, [-1 0]);
            all_heading_angles = [all_heading_angles; heading_worm]; all_distances = [all_distances; distance_worm]; 
            worm_heading_concat  = [worm_heading_concat;  heading_worm]; worm_distance_concat = [worm_distance_concat; distance_worm]; 
            seg_counter = seg_counter + 1;
            if seg_counter == 1
                coordinates_to_plot(col).worm_id = col; coordinates_to_plot(col).x0_worm = x0_worm;
                coordinates_to_plot(col).seg_start_idx = starts(seg); coordinates_to_plot(col).seg_end_idx = ends_(seg);
            else
                coordinates_to_plot(col).seg_start_idx(end + 1) = starts(seg); coordinates_to_plot(col).seg_end_idx(end + 1) = ends_(seg);  
            end
            coordinates_to_plot(col).x_smooth{seg_counter, 1} = smoothedX; coordinates_to_plot(col).y_smooth{seg_counter, 1} = smoothedY;
            coordinates_to_plot(col).rev_index{seg_counter, 1} = reversal_indices; coordinates_to_plot(col).heading_angle{seg_counter} = heading_worm;
        end
        if ~isempty(worm_heading_concat), all_mean_headings_per_worm(col) = mean(worm_heading_concat, 'omitnan'); end
        if ~isempty(worm_distance_concat), all_mean_distances_per_worm(col) = mean(worm_distance_concat, 'omitnan'); end
    end
    if plot_figures == 1, set(gca, 'fontsize', 16); title(plot_title, 'fontweight', 'normal', 'fontsize', 18); end
end

function run_segments = extract_run_segments_from_coordinates_3(coordinates, sampling_rate)
    run_segments = {};
    min_run_seconds = 5;
    min_run_frames  = max(2, round(min_run_seconds * sampling_rate));
    for i = 1:numel(coordinates)
        if ~isfield(coordinates, 'x_smooth') || isempty(coordinates(i).x_smooth), continue; end
        xcell = coordinates(i).x_smooth; ycell = coordinates(i).y_smooth; revcell = coordinates(i).rev_index;
        if ~iscell(xcell), xcell = {xcell}; end; if ~iscell(ycell), ycell = {ycell}; end; if ~iscell(revcell), revcell = {revcell}; end
        nSeg = numel(xcell);
        for s = 1:nSeg
            x = xcell{s}; y = ycell{s};
            if s <= numel(revcell), rev_idx = revcell{s}; else, rev_idx = []; end
            if isempty(x) || isempty(y), continue; end
            x = x(:); y = y(:); N = numel(x);
            if ~isempty(rev_idx), rev_idx = rev_idx(:); rev_idx = rev_idx(~isnan(rev_idx) & rev_idx >= 1 & rev_idx <= N); end
            rev_idx = unique([1; rev_idx; N]);
            for r = 1:numel(rev_idx)-1
                run_start = rev_idx(r); run_end = rev_idx(r+1);
                if run_end <= run_start || (run_end - run_start + 1) < min_run_frames, continue; end
                x_run = x(run_start:run_end); y_run = y(run_start:run_end);
                mask = ~isnan(x_run) & ~isnan(y_run); x_run = x_run(mask); y_run = y_run(mask);
                if numel(x_run) < 2, continue; end
                run_segments{end+1,1} = [x_run y_run]; 
            end
        end
    end
end

function plot_true_runs_from_origin_topN_2(run_segments, run_color, panel_title, n_runs_to_plot, ylimits, xlimits)
    if nargin < 4 || isempty(n_runs_to_plot), n_runs_to_plot = 20; end
    if nargin < 5, ylimits = []; end
    if nargin < 6, xlimits = []; end
    cla reset; hold on
    n_total = numel(run_segments);
    if n_total == 0
        title(panel_title, 'FontWeight', 'normal'); xlabel('X Position (pixels)'); ylabel('Y Position (pixels)'); box off; return;
    end
    run_lengths = nan(n_total,1);
    for k = 1:n_total
        run_xy = run_segments{k};
        run_lengths(k) = sum(sqrt(diff(run_xy(:,1)).^2 + diff(run_xy(:,2)).^2), 'omitnan');
    end
    [~, sort_idx] = sort(run_lengths, 'descend');
    top_idx = sort_idx(1:min(n_runs_to_plot, n_total));
    all_x = []; all_y = [];
    for ii = 1:numel(top_idx)
        k = top_idx(ii); run_xy = run_segments{k};
        x = run_xy(:,1); y = run_xy(:,2); x = x(:); y = y(:);
        valid = ~(isnan(x) | isnan(y)); x = x(valid); y = y(valid);
        if numel(x) < 2, continue; end
        x_shift = x - x(1); y_shift = y - y(1);
        plot(x_shift, y_shift, '-', 'Color', run_color, 'LineWidth', 2);
        all_x = [all_x; x_shift]; all_y = [all_y; y_shift]; 
    end
    xlabel('X Position (pixels)'); ylabel('Y Position (pixels)'); title(panel_title, 'FontWeight', 'normal'); box off;
    if ~isempty(ylimits), ylim(ylimits); elseif ~isempty(all_y), ylim([min(all_y) max(all_y)]); end
    if ~isempty(xlimits), xlim(xlimits); elseif ~isempty(all_x), xlim([min(all_x) max(all_x)]); end
    xl = xlim; yl = ylim; plot([xl(1) xl(2)], [0 0], 'k-', 'LineWidth', 1); plot([0 0], [yl(1) yl(2)], 'k-', 'LineWidth', 1);
    xlim(xl); ylim(yl);
end

function [all_run_headings, all_run_distances, all_run_durations, all_run_worm_ids] = ...
    extract_runs_from_coordinates_3(coordinates, sampling_rate)
    all_run_headings  = []; all_run_distances = []; all_run_durations = []; all_run_worm_ids  = [];
    min_run_seconds = 5; min_run_frames  = max(2, round(min_run_seconds * sampling_rate));
    for i = 1:numel(coordinates)
        if ~isfield(coordinates, 'x_smooth') || isempty(coordinates(i).x_smooth), continue; end
        worm_id = coordinates(i).worm_id; if isempty(worm_id), worm_id = NaN; end
        xcell = coordinates(i).x_smooth; ycell = coordinates(i).y_smooth; revcell = coordinates(i).rev_index;
        if ~iscell(xcell), xcell = {xcell}; end; if ~iscell(ycell), ycell = {ycell}; end; if ~iscell(revcell), revcell = {revcell}; end
        nSeg = numel(xcell);
        for s = 1:nSeg
            x = xcell{s}; y = ycell{s};
            if s <= numel(revcell), rev_idx = revcell{s}; else, rev_idx = []; end
            if isempty(x) || isempty(y), continue; end
            N = numel(x);
            if ~isempty(rev_idx), rev_idx = rev_idx(:); rev_idx = rev_idx(~isnan(rev_idx) & rev_idx >= 1 & rev_idx <= N); end
            rev_idx = unique([1; rev_idx; N]);
            for r = 1:numel(rev_idx)-1
                run_start = rev_idx(r); run_end = rev_idx(r+1);
                if run_end <= run_start || (run_end - run_start + 1) < min_run_frames, continue; end
                x_run = x(run_start:run_end); y_run = y(run_start:run_end);
                mask = ~isnan(x_run) & ~isnan(y_run); x_run = x_run(mask); y_run = y_run(mask);
                if numel(x_run) < 2, continue; end
                run_duration = (run_end - run_start) / sampling_rate;
                run_distance = sum(sqrt(diff(x_run).^2 + diff(y_run).^2), 'omitnan');
                [heading_run, ~] = calculateRelativeHeadingandDistance(x_run, y_run, [-1 0]);
                mean_heading_run = circ_mean_deg(heading_run);
                all_run_durations(end+1,1) = run_duration;      
                all_run_distances(end+1,1) = run_distance;      
                all_run_headings(end+1,1)  = mean_heading_run;  
                all_run_worm_ids(end+1,1)  = worm_id;   
            end
        end
    end
end

function mu = circ_mean_deg(theta_deg)
    th = deg2rad(theta_deg(:));
    mu = rad2deg(atan2(mean(sin(th),'omitnan'), mean(cos(th),'omitnan')));
end

function [diameter_of_polar_plot,directionality,distance_at_zero,distance_at_one_eighty,normalized_directionality,isotherm_diameter]=polar_distance_plot_color_radius_lim_3(heading_deg, distance_bout, bin_width,plot_title,line_color,rlim_min,rlim_max)
    if nargin < 3, bin_width = 10; end
    if mod(bin_width, 2) ~= 0, error('InputError:OddValue', 'The entered value of bin width must be an even number.'); end
    edges = -180:bin_width:180; centers = edges(1:end-1) + bin_width/2;
    binned_distance = zeros(size(centers));
    for i = 1:length(centers)
        in_bin = heading_deg >= edges(i) & heading_deg < edges(i+1);
        binned_distance(i) = mean(distance_bout(in_bin), 'omitnan'); 
    end
    theta_rad = deg2rad(centers); theta_rad = [theta_rad, theta_rad(1)]; binned_distance = [binned_distance, binned_distance(1)];
    polarplot(theta_rad, binned_distance, 'LineWidth', 3,'Color',line_color); rlim([rlim_min rlim_max]);
    title(plot_title, 'fontsize',14,'FontWeight', 'normal');
    pax=gca; pax.ThetaZeroLocation="left"; pax.GridColor = 'k'; pax.GridAlpha = 0.5; pax.MinorGridAlpha = 0.5; pax.Color = 'none';              
    index_for_zero = round((length(binned_distance) + 1) / 2); distance_at_zero=binned_distance(index_for_zero);
    index_for_one_eighty = length(binned_distance); distance_at_one_eighty=binned_distance(index_for_one_eighty);
    index_for_ninty= round(index_for_zero+((index_for_one_eighty-index_for_zero)/2)); distance_at_ninty=binned_distance(index_for_ninty);
    [~, index_for_two_seventy] = min(abs(centers - (-90))); distance_at_two_seventy = binned_distance(index_for_two_seventy);
    diameter_of_polar_plot=distance_at_zero+distance_at_one_eighty; directionality=distance_at_zero-distance_at_one_eighty;
    isotherm_diameter=distance_at_ninty+distance_at_two_seventy; normalized_directionality=directionality./distance_at_ninty;
end

% NEW HELPER: Runs unpaired stats for pooled runs
function do_unpaired_stats(comp_name, data1, data2)
    v1 = data1(~isnan(data1));
    v2 = data2(~isnan(data2));
    n1 = length(v1);
    n2 = length(v2);
    
    fprintf('==============================\n');
    fprintf('%s\n', comp_name);
    fprintf('Total Runs (Group 1): %d\n', n1);
    fprintf('Total Runs (Group 2): %d\n', n2);
    
    if n1 < 3 || n2 < 3
        fprintf('Not enough runs for a valid statistical test.\n==============================\n\n');
        return;
    end
    
    % Lillietest for normality
    p1 = NaN; p2 = NaN;
    if var(v1) > 0 && n1 >= 4, [~, p1] = lillietest(v1); end
    if var(v2) > 0 && n2 >= 4, [~, p2] = lillietest(v2); end
    
    fprintf('Normality p values: %.4f, %.4f\n', p1, p2);
    
    % Choose test (Rank-Sum for independent groups)
    if (isnan(p1) || p1 > 0.05) && (isnan(p2) || p2 > 0.05)
        test_name = 'Two-sample t-test (Unpaired)';
        [~, p, ~, stats] = ttest2(v1, v2);
        stat_val = stats.tstat;
    else
        test_name = 'Wilcoxon rank-sum test (Mann-Whitney U)';
        [p, ~, stats] = ranksum(v1, v2);
        stat_val = stats.ranksum;
    end
    
    m1 = mean(v1); med1 = median(v1); sem1 = std(v1)/sqrt(n1);
    m2 = mean(v2); med2 = median(v2); sem2 = std(v2)/sqrt(n2);
    
    fprintf('Test selected: %s\n', test_name);
    fprintf('Test statistic: %.4f\n', stat_val);
    fprintf('p value: %.4e\n', p);
    fprintf('Group 1 - Mean: %.2f, Median: %.2f, SEM: %.2f\n', m1, med1, sem1);
    fprintf('Group 2 - Mean: %.2f, Median: %.2f, SEM: %.2f\n', m2, med2, sem2);
    fprintf('==============================\n\n');
end
