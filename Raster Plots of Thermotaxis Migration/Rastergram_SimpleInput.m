% generate_rastergram.m

% clear; clc;

%% === USER INPUT ===
% load('your_data.mat', 'x_mat', 'y_mat');

% number of columns to select
n_select = 50;

% total number of columns
n_cols = size(shallow_x_wt, 2);

% set random seed for reproducibility
rng(1);
% randomly select 50 unique column indices
col_idx = randperm(n_cols, n_select);

% extract the selected columns
shallow_x_subset = shallow_x_wt(:, col_idx);

x_mat = shallow_x_subset;

%% === OPTIONS ===
useTempConversion = true;
num_time_bins = 60;   % <<< number of time bins


%% === COLORS ===
hex_colors = {"#0061a0", "#b7babf", "#ef4041"};
rgb_colors = hex2rgb(hex_colors);

color_low  = rgb_colors(1,:);  % blue (15°C)
color_mid  = rgb_colors(2,:);  % gray (center)
color_high = rgb_colors(3,:);  % red (25°C)
%% === TEMP CALIBRATION (optional) ===
% Get extreme values from data
%Shallow gradient boundary
x_min = 251.2803;
x_max = 2230;
% %Steep gradient boundary
% x_min = 1005;
% x_max = 1726;

% %Use automatic boundary
% x_min = min(x_mat(:), [], 'omitnan');
% x_max = max(x_mat(:), [], 'omitnan');


% Target temperature range
% T_min = 15;
% T_max = 25;
T_min = 18;
T_max = 22;

%% === PREP ===
[num_rows, num_worms] = size(x_mat);

% Time in minutes
time_min = (1:num_rows) / 120;

%% === TIME BINNING (row-based, correct) ===
t_min = min(time_min);
t_max = max(time_min);

edges = linspace(t_min, t_max, num_time_bins + 1);

% Assign each row to a bin
[~,~,time_bin_per_row] = histcounts(time_min, edges);

% Fix last edge issue
time_bin_per_row(time_bin_per_row == 0) = num_time_bins;

%% === DATA-DRIVEN TEMPERATURE CONVERSION ===
if useTempConversion

    % Linear mapping: x → temperature
    slope = (T_max - T_min) / (x_max - x_min);
    intercept = T_min - slope * x_min;

    x_plot = slope * x_mat + intercept;

    % Ticks directly in temperature space
    xticktemp = T_min:T_max;
    xtickloc = xticktemp;

    xlabel_str = ['Temperature (' char(176) 'C)'];
else
    x_plot = x_mat;

    xtickloc = linspace(min(x_plot(:)), max(x_plot(:)), 10);
    xticktemp = round(xtickloc,2);

    xlabel_str = 'X position (pixel)';
end

%% === BUILD RASTER ===
total_spikes = numel(x_plot);

X = NaN(3, total_spikes);
Y = NaN(3, total_spikes);

sp = 0;

for i = 1:num_rows
    bin_id = time_bin_per_row(i);

    for j = 1:num_worms
        x_val = x_plot(i,j);

        if isnan(x_val)
            continue;
        end

        sp = sp + 1;

        X(:,sp) = [x_val; x_val; NaN];
        Y(:,sp) = [bin_id - 0.5; bin_id + 0.5; NaN];
    end
end

X = X(:,1:sp);
Y = Y(:,1:sp);

Xv = X(:);
Yv = Y(:);

%% === PLOT ===
figure;
plot(Xv, Yv, 'k');
% plot(Xv, Yv, 'Color', color_mid);
hold off;

set(gca, 'YDir', 'reverse');

xlabel(xlabel_str);
ylabel('Time (min)');

%% === TICKS ===
% xticks(xtickloc);
% xticklabels(xticktemp);
xticks(15:25);
xticklabels(15:25);

% Map bins → real minutes
yticks(linspace(1, num_time_bins, 7));
yticklabels(round(linspace(t_min, t_max, 7),1));

%% === LIMITS ===
% xoffset = range(Xv(~isnan(Xv))) * 0.05;
% yoffset = num_time_bins * 0.05;
% 
% xlim([min(Xv(~isnan(Xv))) - xoffset, max(Xv(~isnan(Xv))) + xoffset]);
% ylim([1 - yoffset, num_time_bins + yoffset]);

yoffset = num_time_bins * 0.05;

% amount of padding (in temperature units)
xpad = (T_max - T_min) * 0.05;   % 5% padding (adjust if needed)
% xpad = 3.5;

% Force x-axis to match temperature range
if useTempConversion
    xlim([T_min - xpad, T_max + xpad]);
else
    xoffset = range(Xv(~isnan(Xv))) * 0.05;
    xlim([min(Xv(~isnan(Xv))) - xoffset, max(Xv(~isnan(Xv))) + xoffset]);
end

ylim([1 - yoffset, num_time_bins + yoffset]);

%% === OPTIONAL TEMP LINES ===
if useTempConversion
    hold on;
    xline(T_min, '--', 'Color', color_low, ...
        'Label', [num2str(T_min) char(176) 'C'], ...
        'LineWidth', 2, ...
        'LabelHorizontalAlignment','center', ...
        'LabelVerticalAlignment','middle');

    xline(T_max, '--', 'Color', color_high, ...
        'Label', [num2str(T_max) char(176) 'C'], ...
        'LineWidth', 2, ...
        'LabelHorizontalAlignment','center', ...
        'LabelVerticalAlignment','middle');

    xline((T_min+T_max)/2, '--', 'Color', color_mid, ...
        'Label', [num2str((T_min+T_max)/2) char(176) 'C'], ...
        'LineWidth', 2, ...
        'LabelHorizontalAlignment','center', ...
        'LabelVerticalAlignment','middle');
    hold off;
end

%% === STYLE ===
set(gca, 'FontName', 'Arial', 'FontSize', 14);
set(gca, 'Color', 'none');

title(['Rastergram (' num2str(num_time_bins) ' bins)']);