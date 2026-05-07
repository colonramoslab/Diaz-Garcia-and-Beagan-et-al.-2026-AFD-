%% plot_n2_runs_and_polar_summary_by_temperature_bin.m
% Plot N2 thermotaxis run trajectories and polar heading-distance summaries
% across defined temperature bins.
%
% This script generates two figures:
%   1) Run trajectories aligned to the origin
%   2) Polar heading-distance summaries
%
% Required input files:
%   - merged_X.mat containing variable merged_X
%   - merged_Y.mat containing variable merged_Y

clc;
clear;
close all;

%% =========================
% USER SETTINGS
% =========================

xMatPath = "/Users/sk3526/Documents/Yale/modelling/figure_5_AFD_paper/datasets/merged_data_first_15_mins/merged_X.mat";
yMatPath = "/Users/sk3526/Documents/Yale/modelling/figure_5_AFD_paper/datasets/merged_data_first_15_mins/merged_Y.mat";

saveFigures = false;
saveDir = "/Users/sk3526/Desktop/github_figures";

if saveFigures && ~exist(saveDir, "dir")
    mkdir(saveDir);
end

%% =========================
% LOAD DATA
% =========================

xData = load(xMatPath, "merged_X");
yData = load(yMatPath, "merged_Y");

xCoordinates = xData.merged_X;
yCoordinates = yData.merged_Y;

%% =========================
% ANALYSIS PARAMETERS
% =========================

rawSamplingRateHz = 2;              % frames/sec in original coordinate matrices
segmentationSamplingRateHz = 1;     % frames/sec after downsampling inside segmentation

analysisStartMin = 0;
analysisEndMin   = 15;

startFrame = rawSamplingRateHz * 60 * analysisStartMin + 1;
endFrame   = rawSamplingRateHz * 60 * analysisEndMin;

% Pixel-to-temperature calibration.
pixelRange = [900 1800];
temperatureRange = [15 25];

% Spatial bins in pixels.
% For pixelRange = [900 1800] and temperatureRange = [15 25],
% xBinEdges = 1350:180:1800 corresponds approximately to:
%   20-22 °C and 22-24 °C.
xBinEdges = 1350:180:1800;

% Pixel-to-mm conversion.
% 900 pixels = 100 mm.
pixelsToMm = 100 / 900;

% Plot settings.
numRunsToDisplay = 40;
runColor = [1 0 0];
fontSize = 18;

trackYLimitsPixels = [-100 100];
trackXLimitsPixels = [-60 120];

polarBinWidthDeg = 10;
polarRadiusLimitsMm = [0 6];
polarRadialTicksMm = [0 3 6];

%% =========================
% FIGURE 1: RUN TRAJECTORIES FROM ORIGIN
% =========================

figure;
tiledlayout(1, numel(xBinEdges) - 1, ...
    "TileSpacing", "compact", ...
    "Padding", "compact");

for binIdx = 1:(numel(xBinEdges) - 1)

    xBinLow  = xBinEdges(binIdx);
    xBinHigh = xBinEdges(binIdx + 1);

    tempLow  = pixelToTemperature(xBinLow,  pixelRange, temperatureRange);
    tempHigh = pixelToTemperature(xBinHigh, pixelRange, temperatureRange);

    [xBin, yBin] = extractCoordinatesInXBin( ...
        xCoordinates, yCoordinates, ...
        startFrame, endFrame, ...
        xBinLow, xBinHigh);

    trackSegments = segmentSmoothedTracksAndDetectReversals( ...
        xBin, yBin, ...
        "N2", "k", false);

    runSegments = extractRunSegmentsBetweenReversals( ...
        trackSegments, segmentationSamplingRateHz);

    nexttile;

    plotRunSegmentsFromOrigin( ...
        runSegments, ...
        runColor, ...
        sprintf("Steep: %.1f-%.1f °C", tempLow, tempHigh), ...
        numRunsToDisplay, ...
        trackYLimitsPixels, ...
        trackXLimitsPixels);

    formatTrajectoryAxesInMillimeters(gca, pixelsToMm, fontSize);
end

if saveFigures
    saveas(gcf, fullfile(saveDir, "n2_run_tracks_by_temperature_bin.svg"));
    saveas(gcf, fullfile(saveDir, "n2_run_tracks_by_temperature_bin.pdf"));
end

%% =========================
% FIGURE 2: POLAR HEADING-DISTANCE SUMMARIES
% =========================

figure;

numBins = numel(xBinEdges) - 1;

for binIdx = 1:numBins

    polarAx = createPolarAxesAtSubplotPosition(1, numBins, binIdx);

    xBinLow  = xBinEdges(binIdx);
    xBinHigh = xBinEdges(binIdx + 1);

    tempLow  = pixelToTemperature(xBinLow,  pixelRange, temperatureRange);
    tempHigh = pixelToTemperature(xBinHigh, pixelRange, temperatureRange);

    [xBin, yBin] = extractCoordinatesInXBin( ...
        xCoordinates, yCoordinates, ...
        startFrame, endFrame, ...
        xBinLow, xBinHigh);

    trackSegments = segmentSmoothedTracksAndDetectReversals( ...
        xBin, yBin, ...
        "N2", "k", false);

    [runHeadingsDeg, runDistancesPixels] = computeRunHeadingDistanceDuration( ...
        trackSegments, segmentationSamplingRateHz);

    runDistancesMm = runDistancesPixels * pixelsToMm;

    fprintf("Bin %d: %.1f-%.1f °C, %d headings, %d distances\n", ...
        binIdx, tempLow, tempHigh, numel(runHeadingsDeg), numel(runDistancesMm));

    if isempty(runHeadingsDeg) || isempty(runDistancesMm)
        title(polarAx, sprintf("%.1f-%.1f °C\n(no runs)", tempLow, tempHigh), ...
            "FontWeight", "normal", ...
            "FontSize", fontSize);
        continue;
    end

    axes(polarAx);

    plotPolarHeadingDistanceSummary( ...
        runHeadingsDeg, ...
        runDistancesMm, ...
        polarBinWidthDeg, ...
        runColor, ...
        polarRadiusLimitsMm);

    formatPolarHeadingAxes( ...
        polarAx, ...
        fontSize, ...
        polarRadialTicksMm, ...
        sprintf("%.1f-%.1f °C", tempLow, tempHigh));
end

if saveFigures
    saveas(gcf, fullfile(saveDir, "n2_polar_heading_distance_by_temperature_bin.svg"));
    saveas(gcf, fullfile(saveDir, "n2_polar_heading_distance_by_temperature_bin.pdf"));
end

%% ========================================================================
% LOCAL HELPER FUNCTIONS
% ========================================================================

function [xBin, yBin] = extractCoordinatesInXBin( ...
    xCoordinates, yCoordinates, startFrame, endFrame, xBinLow, xBinHigh)
% Extract coordinates within a time window and retain only positions inside
% the specified x-coordinate bin.

    xBin = xCoordinates(startFrame:endFrame, :);
    yBin = yCoordinates(startFrame:endFrame, :);

    inBin = xBin >= xBinLow & xBin <= xBinHigh;

    xBin(~inBin) = NaN;
    yBin(~inBin) = NaN;
end

function temperature = pixelToTemperature(xPixel, pixelRange, temperatureRange)
% Convert x-position in pixels to temperature using linear interpolation.

    temperature = temperatureRange(1) + ...
        (xPixel - pixelRange(1)) * ...
        (temperatureRange(2) - temperatureRange(1)) / ...
        (pixelRange(2) - pixelRange(1));
end

function formatTrajectoryAxesInMillimeters(ax, pixelsToMm, fontSize)
% Display trajectory axis tick labels in millimeters while keeping plotted
% coordinates in pixels.

    xTicksMm = [-5 0 5 10];
    yTicksMm = [-10 -5 0 5 10];

    ax.FontSize = fontSize;

    ax.XTick = xTicksMm / pixelsToMm;
    ax.YTick = yTicksMm / pixelsToMm;

    ax.XTickLabel = string(xTicksMm);
    ax.YTickLabel = string(yTicksMm);

    xlabel(ax, "X Position (mm)");
    ylabel(ax, "Y Position (mm)");
end

function polarAx = createPolarAxesAtSubplotPosition(numRows, numCols, tileIdx)
% Create polar axes at the position of a standard subplot tile.

    tempAx = subplot(numRows, numCols, tileIdx);
    axPosition = get(tempAx, "Position");
    delete(tempAx);

    polarAx = polaraxes("Position", axPosition);
    hold(polarAx, "on");
end

function formatPolarHeadingAxes(polarAx, fontSize, radialTicksMm, plotTitle)
% Apply consistent formatting to polar heading-distance plots.

    polarAx.ThetaZeroLocation = "left";
    polarAx.ThetaDir = "counterclockwise";
    polarAx.ThetaLim = [-180 180];
    polarAx.ThetaTick = [-180 -90 0 90];

    polarAx.RTick = radialTicksMm;
    polarAx.FontSize = fontSize;
    polarAx.GridColor = "k";
    polarAx.GridAlpha = 1;

    title(polarAx, plotTitle, ...
        "FontWeight", "normal", ...
        "FontSize", fontSize);
end

function trackSegments = segmentSmoothedTracksAndDetectReversals( ...
    xCoordinates, yCoordinates, plotTitle, plotColor, showDiagnosticPlot)
% Smooth worm tracks, split them into contiguous valid segments, detect
% reversals in each segment, and return a structured representation of each
% worm's segmented trajectory.

    smoothingWindowFrames = 7;
    downsampleFactor = 2;
    smoothingPolynomialOrder = 3;

    reversalAngleThresholdDeg = 120;
    minReversalDurationFrames = 2;
    reversalWindowFrames = 5;

    if mod(smoothingWindowFrames, 2) == 0
        smoothingWindowFrames = smoothingWindowFrames + 1;
    end

    if smoothingWindowFrames < smoothingPolynomialOrder + 1
        smoothingWindowFrames = smoothingPolynomialOrder + 1;

        if mod(smoothingWindowFrames, 2) == 0
            smoothingWindowFrames = smoothingWindowFrames + 1;
        end
    end

    [~, numWorms] = size(xCoordinates);

    emptyTrack = struct( ...
        "wormId", [], ...
        "xSmooth", {{}}, ...
        "ySmooth", {{}}, ...
        "reversalIndices", {{}}, ...
        "headingAngles", {{}}, ...
        "segmentStartFrames", [], ...
        "segmentEndFrames", [], ...
        "firstValidX", []);

    trackSegments = repmat(emptyTrack, 1, numWorms);

    for wormIdx = 1:numWorms

        xWorm = xCoordinates(:, wormIdx);
        yWorm = yCoordinates(:, wormIdx);

        validPoints = ~isnan(xWorm) & ~isnan(yWorm);

        trackSegments(wormIdx).wormId = wormIdx;

        if ~any(validPoints)
            continue;
        end

        firstValidFrame = find(validPoints, 1, "first");
        firstValidX = xWorm(firstValidFrame);

        [segmentStartFrames, segmentEndFrames] = findContiguousValidSegments(validPoints);

        segmentCounter = 0;

        for segIdx = 1:numel(segmentStartFrames)

            xSegment = xWorm(segmentStartFrames(segIdx):segmentEndFrames(segIdx));
            ySegment = yWorm(segmentStartFrames(segIdx):segmentEndFrames(segIdx));

            if numel(xSegment) < max(smoothingWindowFrames, smoothingPolynomialOrder + 1)
                continue;
            end

            xSmooth = sgolayfilt(xSegment, smoothingPolynomialOrder, smoothingWindowFrames);
            ySmooth = sgolayfilt(ySegment, smoothingPolynomialOrder, smoothingWindowFrames);

            xSmooth = downsample(xSmooth, downsampleFactor);
            ySmooth = downsample(ySmooth, downsampleFactor);

            reversalIndices = detectWormReversals( ...
                xSmooth, ySmooth, ...
                "WindowSize", reversalWindowFrames, ...
                "AngleThreshold", reversalAngleThresholdDeg, ...
                "MinReversalDuration", minReversalDurationFrames, ...
                "PlotResults", false);

            [headingAngles, ~] = calculateRelativeHeadingAndDistance( ...
                xSmooth, ySmooth, [1 0]);

            segmentCounter = segmentCounter + 1;

            if segmentCounter == 1
                trackSegments(wormIdx).firstValidX = firstValidX;
            end

            trackSegments(wormIdx).segmentStartFrames(segmentCounter, 1) = segmentStartFrames(segIdx);
            trackSegments(wormIdx).segmentEndFrames(segmentCounter, 1) = segmentEndFrames(segIdx);
            trackSegments(wormIdx).xSmooth{segmentCounter, 1} = xSmooth;
            trackSegments(wormIdx).ySmooth{segmentCounter, 1} = ySmooth;
            trackSegments(wormIdx).reversalIndices{segmentCounter, 1} = reversalIndices;
            trackSegments(wormIdx).headingAngles{segmentCounter, 1} = headingAngles;

            if showDiagnosticPlot
                hold on;
                plot(xSmooth - firstValidX, ySmooth, ...
                    "LineStyle", "-", ...
                    "Color", plotColor, ...
                    "LineWidth", 1);

                scatter(xSmooth - firstValidX, ySmooth, 4, plotColor, "filled");

                if ~isempty(reversalIndices)
                    plot(xSmooth(reversalIndices) - firstValidX, ySmooth(reversalIndices), "ro", ...
                        "MarkerSize", 3, ...
                        "LineWidth", 2, ...
                        "MarkerFaceColor", "r");
                end
            end
        end
    end

    if showDiagnosticPlot
        xlim([-1000 400]);
        xlabel("X Position");
        yticks([]);
        set(gca, "FontSize", 16, "YColor", "none");
        title(plotTitle, "FontWeight", "normal", "FontSize", 18);
    end
end

function [segmentStartFrames, segmentEndFrames] = findContiguousValidSegments(validPoints)
% Find start and end indices of contiguous true regions in a logical vector.

    segmentTransitions = diff([false; validPoints(:); false]);

    segmentStartFrames = find(segmentTransitions == 1);
    segmentEndFrames = find(segmentTransitions == -1) - 1;
end

function reversalIndices = detectWormReversals(x, y, varargin)
% Detect reversal events in a worm trajectory.
%
% Optional name-value parameters:
%   WindowSize          Number of frames before and after the candidate point.
%   AngleThreshold      Minimum angle between past and future vectors.
%   SpeedThreshold      Minimum local speed required for a candidate reversal.
%   MinReversalDuration Minimum number of consecutive reversal frames.
%   PlotResults         Whether to show a diagnostic plot.

    parser = inputParser;

    addRequired(parser, "x", @isnumeric);
    addRequired(parser, "y", @isnumeric);

    addParameter(parser, "WindowSize", 10, ...
        @(v) isnumeric(v) && isscalar(v) && v > 0);

    addParameter(parser, "AngleThreshold", 150, ...
        @(v) isnumeric(v) && isscalar(v) && v > 0 && v <= 180);

    addParameter(parser, "SpeedThreshold", 0.01, ...
        @(v) isnumeric(v) && isscalar(v) && v >= 0);

    addParameter(parser, "MinReversalDuration", 3, ...
        @(v) isnumeric(v) && isscalar(v) && v > 0);

    addParameter(parser, "PlotResults", false, ...
        @(v) islogical(v) || isnumeric(v) && ismember(v, [0 1]));

    parse(parser, x, y, varargin{:});

    windowSize = round(parser.Results.WindowSize);
    angleThresholdDeg = parser.Results.AngleThreshold;
    speedThreshold = parser.Results.SpeedThreshold;
    minReversalDurationFrames = round(parser.Results.MinReversalDuration);
    plotResults = logical(parser.Results.PlotResults);

    x = x(:);
    y = y(:);

    if numel(x) ~= numel(y)
        error("Input vectors x and y must have the same number of elements.");
    end

    coordinates = [x, y];
    numPoints = size(coordinates, 1);
    isReversalPoint = false(numPoints, 1);

    for pointIdx = (windowSize + 1):(numPoints - windowSize)

        localCoordinates = coordinates(pointIdx - windowSize:pointIdx + windowSize, :);

        if any(isnan(localCoordinates(:)))
            continue;
        end

        pastVector = coordinates(pointIdx, :) - coordinates(pointIdx - windowSize, :);
        futureVector = coordinates(pointIdx + windowSize, :) - coordinates(pointIdx, :);

        pastNorm = norm(pastVector);
        futureNorm = norm(futureVector);

        localSpeed = norm(coordinates(pointIdx + 1, :) - coordinates(pointIdx - 1, :)) / 2;

        if pastNorm == 0 || futureNorm == 0 || localSpeed < speedThreshold
            continue;
        end

        cosAngle = dot(pastVector, futureVector) / (pastNorm * futureNorm);
        cosAngle = max(-1, min(1, cosAngle));

        turnAngleDeg = acosd(cosAngle);

        if turnAngleDeg > angleThresholdDeg
            isReversalPoint(pointIdx) = true;
        end
    end

    reversalStarts = find(diff([false; isReversalPoint]) == 1);
    reversalEnds = find(diff([false; isReversalPoint]) == -1) - 1;

    if numel(reversalEnds) < numel(reversalStarts)
        reversalEnds(end + 1) = numPoints;
    end

    reversalDurations = reversalEnds - reversalStarts + 1;
    validReversals = reversalDurations >= minReversalDurationFrames;

    reversalIndices = reversalStarts(validReversals);

    if plotResults
        figure;
        hold on;

        plot(x, y, "k-", "LineWidth", 1.5);

        if ~isempty(reversalIndices)
            plot(x(reversalIndices), y(reversalIndices), "ro", ...
                "MarkerSize", 5, ...
                "LineWidth", 2, ...
                "MarkerFaceColor", "r");
            legend("Trajectory", "Reversal Start");
        else
            legend("Trajectory");
        end

        title("Worm Trajectory with Detected Reversals");
        xlabel("X Position");
        ylabel("Y Position");
        box off;
        hold off;
    end
end

function [headingDeg, stepDistance] = calculateRelativeHeadingAndDistance(x, y, referenceVector)
% Calculate heading relative to a reference vector and stepwise distance.

    if nargin < 3
        referenceVector = [-1 0];
    end

    dx = diff(x);
    dy = diff(y);

    absoluteHeadingDeg = atan2d(dy, dx);
    referenceAngleDeg = atan2d(referenceVector(2), referenceVector(1));

    relativeHeadingDeg = absoluteHeadingDeg - referenceAngleDeg;

    headingDeg = mod(relativeHeadingDeg + 180, 360) - 180;
    stepDistance = sqrt(dx.^2 + dy.^2);
end

function runSegments = extractRunSegmentsBetweenReversals(trackSegments, samplingRateHz)
% Return run trajectory segments, where each run is the path between two
% consecutive reversal indices.

    runSegments = {};

    minRunDurationSec = 5;
    minRunFrames = max(2, round(minRunDurationSec * samplingRateHz));

    for wormIdx = 1:numel(trackSegments)

        if ~isfield(trackSegments, "xSmooth") || isempty(trackSegments(wormIdx).xSmooth)
            continue;
        end

        xCells = trackSegments(wormIdx).xSmooth;
        yCells = trackSegments(wormIdx).ySmooth;
        reversalCells = trackSegments(wormIdx).reversalIndices;

        numSegments = numel(xCells);

        for segIdx = 1:numSegments

            x = xCells{segIdx};
            y = yCells{segIdx};

            if segIdx <= numel(reversalCells)
                reversalIdx = reversalCells{segIdx};
            else
                reversalIdx = [];
            end

            if isempty(x) || isempty(y)
                continue;
            end

            x = x(:);
            y = y(:);
            numPoints = numel(x);

            reversalIdx = cleanReversalIndices(reversalIdx, numPoints);
            reversalIdx = unique([1; reversalIdx; numPoints]);

            for runIdx = 1:(numel(reversalIdx) - 1)

                runStart = reversalIdx(runIdx);
                runEnd = reversalIdx(runIdx + 1);

                if runEnd <= runStart
                    continue;
                end

                if (runEnd - runStart + 1) < minRunFrames
                    continue;
                end

                xRun = x(runStart:runEnd);
                yRun = y(runStart:runEnd);

                validPoints = ~isnan(xRun) & ~isnan(yRun);

                xRun = xRun(validPoints);
                yRun = yRun(validPoints);

                if numel(xRun) < 2
                    continue;
                end

                runSegments{end + 1, 1} = [xRun yRun]; %#ok<AGROW>
            end
        end
    end
end

function [runHeadingsDeg, runDistances, runDurationsSec] = computeRunHeadingDistanceDuration( ...
    trackSegments, samplingRateHz)
% Compute one mean heading, distance, and duration for each run between
% consecutive reversals.

    runHeadingsDeg = [];
    runDistances = [];
    runDurationsSec = [];

    minRunDurationSec = 5;
    minRunFrames = max(2, round(minRunDurationSec * samplingRateHz));

    for wormIdx = 1:numel(trackSegments)

        if ~isfield(trackSegments, "xSmooth") || isempty(trackSegments(wormIdx).xSmooth)
            continue;
        end

        xCells = trackSegments(wormIdx).xSmooth;
        yCells = trackSegments(wormIdx).ySmooth;
        reversalCells = trackSegments(wormIdx).reversalIndices;

        numSegments = numel(xCells);

        for segIdx = 1:numSegments

            x = xCells{segIdx};
            y = yCells{segIdx};

            if segIdx <= numel(reversalCells)
                reversalIdx = reversalCells{segIdx};
            else
                reversalIdx = [];
            end

            if isempty(x) || isempty(y)
                continue;
            end

            x = x(:);
            y = y(:);
            numPoints = numel(x);

            reversalIdx = cleanReversalIndices(reversalIdx, numPoints);
            reversalIdx = unique([1; reversalIdx; numPoints]);

            for runIdx = 1:(numel(reversalIdx) - 1)

                runStart = reversalIdx(runIdx);
                runEnd = reversalIdx(runIdx + 1);

                if runEnd <= runStart
                    continue;
                end

                if (runEnd - runStart + 1) < minRunFrames
                    continue;
                end

                xRun = x(runStart:runEnd);
                yRun = y(runStart:runEnd);

                validPoints = ~isnan(xRun) & ~isnan(yRun);

                xRun = xRun(validPoints);
                yRun = yRun(validPoints);

                if numel(xRun) < 2
                    continue;
                end

                runDurationsSec(end + 1, 1) = (runEnd - runStart) / samplingRateHz; %#ok<AGROW>

                dx = diff(xRun);
                dy = diff(yRun);

                runDistances(end + 1, 1) = sum(sqrt(dx.^2 + dy.^2), "omitnan"); %#ok<AGROW>

                [stepHeadingsDeg, ~] = calculateRelativeHeadingAndDistance(xRun, yRun, [-1 0]);
                runHeadingsDeg(end + 1, 1) = circularMeanDeg(stepHeadingsDeg); %#ok<AGROW>
            end
        end
    end
end

function reversalIdx = cleanReversalIndices(reversalIdx, numPoints)
% Remove invalid reversal indices and constrain them to the trajectory length.

    if isempty(reversalIdx)
        reversalIdx = [];
        return;
    end

    reversalIdx = reversalIdx(:);
    reversalIdx = reversalIdx(~isnan(reversalIdx) & reversalIdx >= 1 & reversalIdx <= numPoints);
end

function plotRunSegmentsFromOrigin( ...
    runSegments, runColor, panelTitle, numRunsToPlot, yLimitsPixels, xLimitsPixels)
% Plot the longest run segments after shifting each run to start at (0, 0).

    if nargin < 4 || isempty(numRunsToPlot)
        numRunsToPlot = 20;
    end

    if nargin < 5
        yLimitsPixels = [];
    end

    if nargin < 6
        xLimitsPixels = [];
    end

    cla reset;
    hold on;

    numRunsTotal = numel(runSegments);

    if numRunsTotal == 0
        title(panelTitle, "FontWeight", "normal");
        xlabel("X Position (pixels)");
        ylabel("Y Position (pixels)");
        box off;
        return;
    end

    runLengths = nan(numRunsTotal, 1);

    for runIdx = 1:numRunsTotal
        runXY = runSegments{runIdx};

        dx = diff(runXY(:, 1));
        dy = diff(runXY(:, 2));

        runLengths(runIdx) = sum(sqrt(dx.^2 + dy.^2), "omitnan");
    end

    [~, sortIdx] = sort(runLengths, "descend");
    topRunIdx = sortIdx(1:min(numRunsToPlot, numRunsTotal));

    allX = [];
    allY = [];

    for idx = 1:numel(topRunIdx)

        runXY = runSegments{topRunIdx(idx)};

        x = runXY(:, 1);
        y = runXY(:, 2);

        validPoints = ~isnan(x) & ~isnan(y);

        x = x(validPoints);
        y = y(validPoints);

        if numel(x) < 2
            continue;
        end

        xShifted = x - x(1);
        yShifted = y - y(1);

        plot(xShifted, yShifted, "-", ...
            "Color", runColor, ...
            "LineWidth", 2);

        allX = [allX; xShifted]; %#ok<AGROW>
        allY = [allY; yShifted]; %#ok<AGROW>
    end

    xlabel("X Position (pixels)");
    ylabel("Y Position (pixels)");
    title(panelTitle, "FontWeight", "normal");
    box off;

    if ~isempty(yLimitsPixels)
        ylim(yLimitsPixels);
    elseif ~isempty(allY)
        ylim([min(allY) max(allY)]);
    end

    if ~isempty(xLimitsPixels)
        xlim(xLimitsPixels);
    elseif ~isempty(allX)
        xlim([min(allX) max(allX)]);
    end

    xLimits = xlim;
    yLimits = ylim;

    plot([xLimits(1) xLimits(2)], [0 0], "k-", "LineWidth", 1);
    plot([0 0], [yLimits(1) yLimits(2)], "k-", "LineWidth", 1);

    xlim(xLimits);
    ylim(yLimits);
end

function plotPolarHeadingDistanceSummary( ...
    headingDeg, distance, binWidthDeg, lineColor, radiusLimits)
% Plot mean run distance as a function of run heading in polar coordinates.

    if nargin < 3 || isempty(binWidthDeg)
        binWidthDeg = 10;
    end

    if nargin < 4 || isempty(lineColor)
        lineColor = [1 0 0];
    end

    if nargin < 5 || isempty(radiusLimits)
        radiusLimits = [0 max(distance, [], "omitnan")];
    end

    if mod(binWidthDeg, 2) ~= 0
        error("InputError:OddValue", "The bin width must be an even number.");
    end

    binEdgesDeg = -180:binWidthDeg:180;
    binCentersDeg = binEdgesDeg(1:end - 1) + binWidthDeg / 2;

    meanDistanceByHeading = nan(size(binCentersDeg));

    for binIdx = 1:numel(binCentersDeg)
        inBin = headingDeg >= binEdgesDeg(binIdx) & headingDeg < binEdgesDeg(binIdx + 1);
        meanDistanceByHeading(binIdx) = mean(distance(inBin), "omitnan");
    end

    thetaRad = deg2rad(binCentersDeg);

    thetaRad = [thetaRad thetaRad(1)];
    meanDistanceByHeading = [meanDistanceByHeading meanDistanceByHeading(1)];

    polarplot(thetaRad, meanDistanceByHeading, ...
        "LineWidth", 3, ...
        "Color", lineColor);

    rlim(radiusLimits);

    polarAx = gca;
    polarAx.ThetaZeroLocation = "left";
    polarAx.GridColor = "k";
    polarAx.GridAlpha = 0.5;
    polarAx.MinorGridAlpha = 0.5;
    polarAx.Color = "none";
end

function meanAngleDeg = circularMeanDeg(thetaDeg)
% Circular mean of angles in degrees.

    thetaRad = deg2rad(thetaDeg(:));

    meanAngleDeg = rad2deg(atan2( ...
        mean(sin(thetaRad), "omitnan"), ...
        mean(cos(thetaRad), "omitnan")));
end