clc
clear
close all

load("path/merged_X.mat")
load("path/merged_Y.mat")

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

    [run_headings, run_distances, ~] = extract_runs_from_coordinates_3( ...
        coordinates_to_plot, sampling_rate_after_downsample);

    run_distances = run_distances * pix_to_mm;

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

%% =========================
% HELPER FUNCTION
% =========================

function x_pix = temp_to_pix(T, pix_range, temp_range)
    x_pix = pix_range(1) + (T - temp_range(1)) * ...
        (pix_range(2) - pix_range(1)) / (temp_range(2) - temp_range(1));
end

function [all_heading_angles, all_distances, all_mean_headings_per_worm, all_mean_distances_per_worm, coordinates_to_plot] = ...
    visualize_tracks_runs_and_turns_5_Tc25(x_coord, y_coord, plot_title, ttx_color, plot_figures)
% x_coord, y_coord are T x M with possible interior NaNs.
% Segment aware, but outputs are aligned to worms (M outputs).
% Per worm means are computed across all its valid segments.
% Plotting recenters by the worm's first valid x, not per segment.

    % Parameters
    windowWidth = 7;           % Savitzky Golay window, must be odd
    downsampling_order = 2;
    polynomial_order = 3;      % sgolay poly order
    reversal_angle_threshold = 120;
    min_reversal_ruration = 2;
    window_size = 5;

    % Ensure valid sgolay settings
    if mod(windowWidth, 2) == 0
        windowWidth = windowWidth + 1;
    end
    if windowWidth < polynomial_order + 1
        windowWidth = polynomial_order + 1;
        if mod(windowWidth, 2) == 0
            windowWidth = windowWidth + 1;
        end
    end

    % Outputs
    all_heading_angles = [];
    all_distances = [];

    % Preallocate per worm arrays to fixed length M
    [T, M] = size(x_coord);
    all_mean_headings_per_worm = nan(M, 1);
    all_mean_distances_per_worm = nan(M, 1);

    % Prepare coordinates_to_plot as 1xM struct with segment cells
    emptyEntry = struct( ...
        'worm_id', [], ...
        'x_smooth', {{}}, ...
        'y_smooth', {{}}, ...
        'rev_index', {{}}, ...
        'heading_angle', {{}}, ...
        'seg_start_idx', [], ...
        'seg_end_idx', [], ...
        'x0_worm', [] ...
    );
    coordinates_to_plot = repmat(emptyEntry, 1, M);

    for col = 1:M
        xcol = x_coord(:, col);
        ycol = y_coord(:, col);

        good = ~isnan(xcol) & ~isnan(ycol);
        if ~any(good)
            % Leave NaNs in mean arrays and an empty struct entry
            coordinates_to_plot(col).worm_id = col;
            continue
        end

        % Per worm recentering offset: first valid x of the worm
        first_valid_idx = find(good, 1, 'first');
        x0_worm = xcol(first_valid_idx);

        % Find contiguous segments of valid data
        d = diff([false; good; false]);
        starts = find(d == 1);
        ends_  = find(d == -1) - 1;

        % Accumulators for this worm
        worm_heading_concat = [];
        worm_distance_concat = [];

        seg_counter = 0;

        for seg = 1:numel(starts)
            xs = xcol(starts(seg):ends_(seg));
            ys = ycol(starts(seg):ends_(seg));

            % Skip too short segments for the chosen filter parameters
            if numel(xs) < max(windowWidth, polynomial_order + 1)
                continue
            end

            % Smooth
            smoothedX = sgolayfilt(xs, polynomial_order, windowWidth);
            smoothedY = sgolayfilt(ys, polynomial_order, windowWidth);

            % Downsample
            smoothedX = downsample(smoothedX, downsampling_order);
            smoothedY = downsample(smoothedY, downsampling_order);

            % Reversal detection on downsampled trajectories
            reversal_indices = findWormReversals_2( ...
                smoothedX, smoothedY, ...
                'WindowSize', window_size, ...
                'AngleThreshold', reversal_angle_threshold, ...
                'MinReversalDuration', min_reversal_ruration, ...
                'PlotResults', false);

            % Plot per segment but recentered by the worm's first valid x
            if plot_figures == 1
                hold on
                plot(smoothedX - x0_worm, smoothedY, 'LineStyle', '-', 'Color', ttx_color, 'LineWidth', 1)
                scatter(smoothedX - x0_worm, smoothedY, 4, ttx_color, 'filled')
                if ~isempty(reversal_indices)
                    plot(smoothedX(reversal_indices) - x0_worm, smoothedY(reversal_indices), 'ro', ...
                        'MarkerSize', 3, 'LineWidth', 2, 'MarkerFaceColor', 'r');
                end
                xlim([-1000 400])
                xlabel('X-coordinates'); yticks([]); set(gca, 'fontsize', 16); set(gca, 'YColor', 'none');
            end

            % Kinematics
            [heading_worm, distance_worm] = calculateRelativeHeadingandDistance(smoothedX, smoothedY, [1 0]);

            % Accumulate global pooled distributions
            all_heading_angles = [all_heading_angles; heading_worm]; %#ok<AGROW>
            all_distances      = [all_distances;      distance_worm]; %#ok<AGROW>

            % Accumulate per worm for later averaging
            worm_heading_concat  = [worm_heading_concat;  heading_worm];  %#ok<AGROW>
            worm_distance_concat = [worm_distance_concat; distance_worm]; %#ok<AGROW>

            % Store per segment into the worm entry
            seg_counter = seg_counter + 1;
            if seg_counter == 1
                coordinates_to_plot(col).worm_id = col;
                coordinates_to_plot(col).x0_worm = x0_worm;
                coordinates_to_plot(col).seg_start_idx = starts(seg);
                coordinates_to_plot(col).seg_end_idx   = ends_(seg);
            else
                coordinates_to_plot(col).seg_start_idx(end + 1) = starts(seg); %#ok<AGROW>
                coordinates_to_plot(col).seg_end_idx(end + 1)   = ends_(seg);  %#ok<AGROW>
            end
            coordinates_to_plot(col).x_smooth{seg_counter, 1}   = smoothedX;
            coordinates_to_plot(col).y_smooth{seg_counter, 1}   = smoothedY;
            coordinates_to_plot(col).rev_index{seg_counter, 1}  = reversal_indices;
            coordinates_to_plot(col).heading_angle{seg_counter} = heading_worm;
        end

        % Final per worm means across all its segments
        if ~isempty(worm_heading_concat)
            all_mean_headings_per_worm(col) = mean(worm_heading_concat, 'omitnan');
        end
        if ~isempty(worm_distance_concat)
            all_mean_distances_per_worm(col) = mean(worm_distance_concat, 'omitnan');
        end
    end

    if plot_figures == 1
        set(gca, 'fontsize', 16)
        title(plot_title, 'fontweight', 'normal', 'fontsize', 18)
    end
end

function run_segments = extract_run_segments_from_coordinates_3(coordinates, sampling_rate)
% Returns true run segments, where each run is the path between two reversals.
% Output:
%   run_segments = cell array, each cell contains [x_run y_run]

run_segments = {};

min_run_seconds = 5;
min_run_frames  = max(2, round(min_run_seconds * sampling_rate));

for i = 1:numel(coordinates)
    if ~isfield(coordinates, 'x_smooth') || isempty(coordinates(i).x_smooth)
        continue
    end

    xcell   = coordinates(i).x_smooth;
    ycell   = coordinates(i).y_smooth;
    revcell = coordinates(i).rev_index;

    if ~iscell(xcell),   xcell   = {xcell};   end
    if ~iscell(ycell),   ycell   = {ycell};   end
    if ~iscell(revcell), revcell = {revcell}; end

    nSeg = numel(xcell);

    for s = 1:nSeg
        x = xcell{s};
        y = ycell{s};

        if s <= numel(revcell)
            rev_idx = revcell{s};
        else
            rev_idx = [];
        end

        if isempty(x) || isempty(y)
            continue
        end

        x = x(:);
        y = y(:);
        N = numel(x);

        if isempty(rev_idx)
            rev_idx = [];
        else
            rev_idx = rev_idx(:);
            rev_idx = rev_idx(~isnan(rev_idx) & rev_idx >= 1 & rev_idx <= N);
        end

        % define runs as intervals between consecutive reversals
        rev_idx = unique([1; rev_idx; N]);

        for r = 1:numel(rev_idx)-1
            run_start = rev_idx(r);
            run_end   = rev_idx(r+1);

            if run_end <= run_start
                continue
            end

            if (run_end - run_start + 1) < min_run_frames
                continue
            end

            x_run = x(run_start:run_end);
            y_run = y(run_start:run_end);

            mask = ~isnan(x_run) & ~isnan(y_run);
            x_run = x_run(mask);
            y_run = y_run(mask);

            if numel(x_run) < 2
                continue
            end

            run_segments{end+1,1} = [x_run y_run]; %#ok<AGROW>
        end
    end
end
end


function plot_true_runs_from_origin_topN_2(run_segments, run_color, panel_title, n_runs_to_plot, ylimits, xlimits)

if nargin < 4 || isempty(n_runs_to_plot)
    n_runs_to_plot = 20;
end
if nargin < 5
    ylimits = [];
end
if nargin < 6
    xlimits = [];
end

cla reset
hold on

n_total = numel(run_segments);
if n_total == 0
    title(panel_title, 'FontWeight', 'normal');
    xlabel('X Position (pixels)');
    ylabel('Y Position (pixels)');
    box off;
    return;
end

run_lengths = nan(n_total,1);
for k = 1:n_total
    run_xy = run_segments{k};
    dx = diff(run_xy(:,1));
    dy = diff(run_xy(:,2));
    run_lengths(k) = sum(sqrt(dx.^2 + dy.^2), 'omitnan');
end

[~, sort_idx] = sort(run_lengths, 'descend');
top_idx = sort_idx(1:min(n_runs_to_plot, n_total));

all_x = [];
all_y = [];

for ii = 1:numel(top_idx)
    k = top_idx(ii);
    run_xy = run_segments{k};

    x = run_xy(:,1);
    y = run_xy(:,2);

    x = x(:);
    y = y(:);

    valid = ~(isnan(x) | isnan(y));
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        continue
    end

    x_shift = x - x(1);
    y_shift = y - y(1);

    plot(x_shift, y_shift, '-', 'Color', run_color, 'LineWidth', 2);

    all_x = [all_x; x_shift]; %#ok<AGROW>
    all_y = [all_y; y_shift]; %#ok<AGROW>
end

xlabel('X Position (pixels)');
ylabel('Y Position (pixels)');
title(panel_title, 'FontWeight', 'normal');
box off;

% Set limits after plotting everything
if ~isempty(ylimits)
    ylim(ylimits);
elseif ~isempty(all_y)
    ylim([min(all_y) max(all_y)]);
end

if ~isempty(xlimits)
    xlim(xlimits);
elseif ~isempty(all_x)
    xlim([min(all_x) max(all_x)]);
end

% Draw axis lines AFTER limits are finalized
xl = xlim;
yl = ylim;
plot([xl(1) xl(2)], [0 0], 'k-', 'LineWidth', 1);
plot([0 0], [yl(1) yl(2)], 'k-', 'LineWidth', 1);

xlim(xl);
ylim(yl);
end

function [all_run_headings, all_run_distances, all_run_durations] = ...
    extract_runs_from_coordinates_3(coordinates, sampling_rate)
% As extract_runs_from_coordinates_2, but:
%  - uses circular mean for run headings
%  - skips very short runs (>= 5 s)

all_run_headings  = [];
all_run_distances = [];
all_run_durations = [];

% ---- policy knob: minimum run duration ----
min_run_seconds = 5;
min_run_frames  = max(2, round(min_run_seconds * sampling_rate));

for i = 1:numel(coordinates)
    if ~isfield(coordinates, 'x_smooth') || isempty(coordinates(i).x_smooth)
        continue
    end

    xcell   = coordinates(i).x_smooth;
    ycell   = coordinates(i).y_smooth;
    revcell = coordinates(i).rev_index;

    if ~iscell(xcell),   xcell   = {xcell};   end
    if ~iscell(ycell),   ycell   = {ycell};   end
    if ~iscell(revcell), revcell = {revcell}; end

    nSeg = numel(xcell);
    for s = 1:nSeg
        x = xcell{s};
        y = ycell{s};
        if s <= numel(revcell)
            rev_idx = revcell{s};
        else
            rev_idx = [];
        end
        if isempty(x) || isempty(y), continue; end

        N = numel(x);

        % clean reversal indices and bound to [1, N]
        if isempty(rev_idx)
            rev_idx = [];
        else
            rev_idx = rev_idx(:);
            rev_idx = rev_idx(~isnan(rev_idx) & rev_idx >= 1 & rev_idx <= N);
        end

        % add artificial start and end, ensure sorted unique
        rev_idx = unique([1; rev_idx; N]);

        for r = 1:numel(rev_idx)-1
            run_start = rev_idx(r);
            run_end   = rev_idx(r+1);

            if run_end <= run_start
                continue
            end

            % ---- skip runs shorter than 5 s ----
            if (run_end - run_start + 1) < min_run_frames
                continue
            end

            % extract run segment
            x_run = x(run_start:run_end);
            y_run = y(run_start:run_end);

            % drop interior NaNs pairwise
            mask = ~isnan(x_run) & ~isnan(y_run);
            x_run = x_run(mask);
            y_run = y_run(mask);

            if numel(x_run) < 2
                continue
            end

            % duration in seconds
            run_duration = (run_end - run_start) / sampling_rate;

            % path length
            dx = diff(x_run);
            dy = diff(y_run);
            run_distance = sum(sqrt(dx.^2 + dy.^2), 'omitnan');

            % mean heading over the run (circular mean)
            [heading_run, ~] = calculateRelativeHeadingandDistance(x_run, y_run, [-1 0]);
            mean_heading_run = circ_mean_deg(heading_run);

            % collect
            all_run_durations(end+1,1) = run_duration;      %#ok<AGROW>
            all_run_distances(end+1,1) = run_distance;      %#ok<AGROW>
            all_run_headings(end+1,1)  = mean_heading_run;  %#ok<AGROW>
        end
    end
end
end

function mu = circ_mean_deg(theta_deg)
% Circular mean of angles in degrees; robust to wraparound.
th = deg2rad(theta_deg(:));
mu = rad2deg(atan2(mean(sin(th),'omitnan'), mean(cos(th),'omitnan')));
end

function [diameter_of_polar_plot,directionality,distance_at_zero,distance_at_one_eighty,normalized_directionality,isotherm_diameter]=polar_distance_plot_color_radius_lim_3(heading_deg, distance_bout, bin_width,plot_title,line_color,rlim_min,rlim_max)
    % Inputs:
    % - heading_deg: array of heading angles in degrees (-180 to 180)
    % - distance_bout: array of distances corresponding to each heading
    % - bin_width: angular bin width in degrees (e.g., 10)

    if nargin < 3
        bin_width = 10; % default bin width
    end
    
     if mod(bin_width, 2) ~= 0
        error('InputError:OddValue', 'The entered value of bin width must be an even number.');
     end

    % Bin edges from -180 to 180 degrees
    edges = -180:bin_width:180;
    centers = edges(1:end-1) + bin_width/2;

    % Initialize binned distance
    binned_distance = zeros(size(centers));

    % Bin the data: sum distances in each heading bin
    for i = 1:length(centers)
        in_bin = heading_deg >= edges(i) & heading_deg < edges(i+1);
        binned_distance(i) = mean(distance_bout(in_bin), 'omitnan'); % Use mean or sum as desired
    end

    % Convert to polar coordinates (r, theta)
    theta_rad = deg2rad(centers);  % convert to radians for polar plot

    % Repeat first bin at the end to close the polar plot
    theta_rad = [theta_rad, theta_rad(1)];
    binned_distance = [binned_distance, binned_distance(1)];

    % Create polar plot
    polarplot(theta_rad, binned_distance, 'LineWidth', 3,'Color',line_color);
    % rlim([0 max(binned_distance) * 1.1]);
    rlim([rlim_min rlim_max]);
    title(plot_title, 'fontsize',14,'FontWeight', 'normal');
    pax=gca;
    pax.ThetaZeroLocation="left";
    pax.GridColor = 'k'; % Set grid color to black
    pax.GridAlpha = 0.5;               % Force grid to be fully opaque
    pax.MinorGridAlpha = 0.5;          % Also apply to minor grid lines
    pax.Color = 'none';              % Avoid background blending
    
    index_for_zero = round((length(binned_distance) + 1) / 2);
    distance_at_zero=binned_distance(index_for_zero);

    index_for_one_eighty = length(binned_distance);
    distance_at_one_eighty=binned_distance(index_for_one_eighty);

    index_for_ninty= round(index_for_zero+((index_for_one_eighty-index_for_zero)/2));
    distance_at_ninty=binned_distance(index_for_ninty);

    [~, index_for_two_seventy] = min(abs(centers - (-90)));
    distance_at_two_seventy = binned_distance(index_for_two_seventy);

    diameter_of_polar_plot=distance_at_zero+distance_at_one_eighty;
    directionality=distance_at_zero-distance_at_one_eighty;
    isotherm_diameter=distance_at_ninty+distance_at_two_seventy;
    normalized_directionality=directionality./distance_at_ninty;

end
