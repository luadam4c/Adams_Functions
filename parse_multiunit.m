function varargout = parse_multiunit (vVecs, siMs, varargin)
%% Parses multiunit recordings: detect spikes
% Usage: [parsedParams, parsedData, figs] = parse_multiunit (vVecs, siMs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       vVecs       - original voltage vector(s) in mV
%                   must be a numeric array or a cell array of numeric arrays
%       siMs        - sampling interval in ms
%                   must be a positive vector
%       varargin    - 'PlotFlag': whether to plot traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileBase': base of filename (without extension)
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%                   - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == detect from pulse vector
%                   - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   
% Requires:
%       cd/argfun.m TODO
%       cd/count_samples.m TODO
%       cd/count_vectors.m TODO
%       cd/iscellnumeric.m TODO
%       cd/find_stim_start.m TODO
%       cd/plot_raster.m TODO
%       cd/plot_horizontal_line.m TODO
%       cd/create_logical_array.m
%       cd/create_time_vectors.m
%       cd/extract_elements.m
%       cd/force_column_cell.m
%       cd/compute_axis_limits.m
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-02-19 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
plotFlagDefault = false;
outFolderDefault = pwd;
fileBaseDefault = {};           % set later
stimStartMsDefault = [];        % set later
pulseVectorsDefault = [];       % don't use pulse vectors by default
tVecsDefault = [];              % set later

% TODO
baseWindows = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PulseVectors', pulseVectorsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['PulseVectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vVecs, siMs, varargin{:});
plotFlag = iP.Results.PlotFlag;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
stimStartMs = iP.Results.StimStartMs;
pulseVectors = iP.Results.PulseVectors;
tVecs = iP.Results.tVecs;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vVecs);

% Count the number of samples for each vector
nSamples = count_samples(vVecs);

% Match time vector(s) with sampling interval(s) and number(s) of samples
[tVecs, siMs, nSamples] = match_time_info(tVecs, siMs, nSamples);

% Initialize figures array
figs = gobjects(nVectors + 1, 1);

% Create figure path base
figPathBase = fullfile(outFolder, [fileBase, '_spike_detection']);

% Create a figure title base
figTitleBase = replace(fileBase, '_', '\_');

%% Do the job
% Detect stimulation start time if not provided
%   Otherwise find the corresponding index in the time vector
if isempty(stimStartMs)
    % TODO: Make this a function find_stim_start.m
    if ~isempty(pulseVectors)
        % Parse the pulse vectors
        [pulseParams, pulseData] = ...
            parse_pulse(pulseVectors, 'SamplingIntervalMs', siMs);

        % Use the indices after pulse starts for stimulation start
        idxStimStart = pulseParams{:, 'idxAfterStart'};

        % Use the time vectors 
        stimStartMs = extract_elements(tVecs, 'specific', ...
                                        'Index', idxStimStart);
    else
        error('One of stimStartMs and pulseVectors must be provided!');
    end
else
    % Find the indices of stimulation start
    if ~isempty(tVecs)
        % Use the indices of tVecs with values closest to stimStartMs
        % TODO: find_closest.m
        idxStimStart = find_closest(tVecs, stimStartMs);
    else
        % Assume tVecs start from 0 and use siMs
        idxStimStart = round(stimStartMs ./ siMs);
    end
end

% Construct default baseline windows
if isempty(baseWindows)
    % Get the starting time(s)
    timeStartMs = extract_elements(tVecs, 'first');

    % Use timeStartMs to stimStartMs by default
    baseWindows = transpose([timeStartMs, stimStartMs]);
end

% Force as a cell array of vectors
[vVecs, tVecs, baseWindows] = ...
    argfun(@force_column_cell, vVecs, tVecs, baseWindows);

% Parse all of them in a parfor loop
parsedParamsCell = cell(nVectors, 1);
parsedDataCell = cell(nVectors, 1);
parfor iVec = 1:nVectors
%for iVec = 1:1
    [parsedParamsCell{iVec}, parsedDataCell{iVec}, figs(iVec)] = ...
        parse_multiunit_helper(iVec, vVecs{iVec}, tVecs{iVec}, siMs(iVec), ...
                                idxStimStart(iVec), stimStartMs(iVec), ...
                                baseWindows{iVec}, ...
                                plotFlag, figPathBase, figTitleBase);
end

% Convert to a struct array
%   Note: This removes all entries that are empty
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@(x) [x{:}], parsedParamsCell, parsedDataCell);

% Convert to a table
[parsedParams, parsedData] = ...
    argfun(@(x) struct2table(x, 'AsArray', true), ...
            parsedParamsStruct, parsedDataStruct);

%% Plot raster plot
if plotFlag
    % Modify the figure base
    figPathBaseThis = [figPathBase, '_raster'];

    % Extract the spike times
    spikeTimesMs = parsedData.spikeTimesMs;
    stimStartMs = parsedParams.stimStartMs;
    detectStartMs = parsedParams.detectStartMs;
    firstSpikeMs = parsedParams.firstSpikeMs;

    % Create figure and plot
    figs(nVectors + 1) = figure(2);
    clf
    [hLines, eventTimes, yEnds, yTicksTable] = ...
        plot_raster(spikeTimesMs, 'LineWidth', 0.5);
    vLines = plot_vertical_line(mean(stimStartMs), 'Color', 'g', 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Trace #');
    title(['Spike times for ', figTitleBase]);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(nVectors + 1), figPathBaseThis, ...
                    mean(stimStartMs), mean(detectStartMs), mean(firstSpikeMs));
end

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;
varargout{3} = figs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData, fig] = ...
                parse_multiunit_helper(iVec, vVec, tVec, siMs, ...
                                idxStimStart, stimStartMs, baseWindow, ...
                                plotFlag, figPathBase, figTitleBase)
% Parse a single multiunit recording

% Hard-coded constants
MS_PER_S = 1000;

% Hard-coded parameters
signal2Noise = 4; %2.5
minDelaySamples = 2000;
binWidthMs = 10;
filterWidthMs = 100;

% Modify the figure base
figPathBaseThis = [figPathBase, '_trace', num2str(iVec)];
figTitleBaseThis = [figTitleBase, '\_trace', num2str(iVec)];

% Find the starting index for detecting a spike
idxDetectStart = idxStimStart + minDelaySamples;

% Find the corresponding time
detectStartMs = tVec(idxDetectStart);

% Compute the number of samples
nSamples = numel(vVec);

% Compute all instantaneous slopes in V/s
slopes = diff(vVec) / siMs;

% Compute a baseline slope noise in V/s
baseSlopeNoise = compute_baseline_noise(slopes, tVec(1:(end-1)), baseWindow);

% Compute a slope threshold in V/s
slopeThreshold = baseSlopeNoise * signal2Noise;

% Determine whether each slope is a local maximum
[~, indPeakSlopes] = findpeaks(slopes);
isPeakSlope = create_logical_array(indPeakSlopes, [nSamples - 1, 1]);

% Create all indices minus 1
allIndices = transpose(1:nSamples);

% Detect spikes after idxStimStart + minDelaySamples
isSpike = [false; slopes > slopeThreshold] & [false; isPeakSlope] & ...
            allIndices > idxDetectStart;
idxSpikes = find(isSpike);

% Compute the overall spike count
spikeCountTotal = numel(idxSpikes);

% Index of first spike
idxFirstSpike = idxSpikes(1);

% Store spike times
spikeTimesMs = tVec(idxSpikes);
firstSpikeMs = spikeTimesMs(1);

% Query the maximum and range of vVec after detectStartMs
vVecTrunc = vVec(idxDetectStart:idxDetectStart + 1e5);
vMin = min(vVecTrunc);
vMax = max(vVecTrunc);
vRange = vMax - vMin;

% Query the maximum and range of slope after detectStartMs
slopesTrunc = slopes(idxDetectStart:idxDetectStart + 1e5);
slopeMin = min(slopesTrunc);
slopeMax = max(slopesTrunc);
slopeRange = slopeMax - slopeMin;

% Compute a spike histogram
[spikeCounts, edgesMs] = compute_bins(spikeTimesMs, 'BinWidth', binWidthMs);

% Compute the number of bins
nBins = numel(spikeCounts);

% Record the starting time of the histogram
histLeftMs = edgesMs(1);

% Compute the bin width in seconds
binWidthSec = binWidthMs / MS_PER_S;

% Compute an unnormalized autocorrelogram in Hz^2
autoCorr = xcorr(spikeCounts, 'unbiased') / binWidthSec ^ 2;

% Take just the positive side
acf = autoCorr(nBins:end);

% Compute a normalized autocorrelation function
% autocorr(spikeCounts, nBins - 1);
% acf = autocorr(spikeCounts, nBins - 1);

% Smooth the autocorrelogram with a moving-average filter
acfFiltered = movingaveragefilter(acf, filterWidthMs, binWidthMs);

% Record the amplitude of the primary peak
ampPeak1 = acfFiltered(1);

% Find the index and amplitude of the secondary peak
[peakAmp, peakInd] = findpeaks(acfFiltered, 'MinPeakProminence', 0.1 * ampPeak1);
idxPeak2 = peakInd(1);
ampPeak2 = peakAmp(1);

% Find the amplitude of the first trough
[troughNegAmp, troughInd] = findpeaks(-acfFiltered);
idxTrough1 = troughInd(1);
ampTrough1 = -troughNegAmp(1);

% Compute the average amplitude of first two peaks
ampPeak12 = mean([ampPeak1, ampPeak2]);

% Compute the oscillatory index
oscIndex = (ampPeak12 - ampTrough1) / ampPeak12;

% Compute the oscillation duration
% TODO

% Compute the average spikes per oscillation
% TODO

% Compute the oscillation period
oscPeriodMs = idxPeak2 * binWidthMs;


% Store in outputs
parsedParams.signal2Noise = signal2Noise;
parsedParams.minDelaySamples = minDelaySamples;
parsedParams.binWidthMs = binWidthMs;
parsedParams.filterWidthMs = filterWidthMs;
parsedParams.siMs = siMs;
parsedParams.idxStimStart = idxStimStart;
parsedParams.stimStartMs = stimStartMs;
parsedParams.baseWindow = baseWindow;
parsedParams.baseSlopeNoise = baseSlopeNoise;
parsedParams.slopeThreshold = slopeThreshold;
parsedParams.idxDetectStart = idxDetectStart;
parsedParams.detectStartMs = detectStartMs;
parsedParams.spikeCountTotal = spikeCountTotal;
parsedParams.idxFirstSpike = idxFirstSpike;
parsedParams.firstSpikeMs = firstSpikeMs;
parsedParams.vMin = vMin;
parsedParams.vMax = vMax;
parsedParams.vRange = vRange;
parsedParams.slopeMin = slopeMin;
parsedParams.slopeMax = slopeMax;
parsedParams.slopeRange = slopeRange;
parsedParams.nBins = nBins;
parsedParams.histLeftMs = histLeftMs;
parsedParams.binWidthSec = binWidthSec;
parsedParams.ampPeak1 = ampPeak1;
parsedParams.idxPeak2 = idxPeak2;
parsedParams.ampPeak2 = ampPeak2;
parsedParams.idxTrough1 = idxTrough1;
parsedParams.ampTrough1 = ampTrough1;
parsedParams.ampPeak12 = ampPeak12;
parsedParams.oscIndex = oscIndex;
parsedParams.oscPeriodMs = oscPeriodMs;

parsedData.tVec = tVec;
parsedData.vVec = vVec;
parsedData.slopes = slopes;
parsedData.idxSpikes = idxSpikes;
parsedData.spikeTimesMs = spikeTimesMs;
parsedData.spikeCounts = spikeCounts;
parsedData.edgesMs = edgesMs;
parsedData.autoCorr = autoCorr;
parsedData.acf = acf;
parsedData.acfFiltered = acfFiltered;

%% Plots
% if plotFlag
if plotFlag && iVec == 1
    %% Plot spike detection
    [fig, ax, lines, markers, raster] = ...
        plot_spike_detection(tVec, vVec, slopes, idxSpikes, ...
                            baseSlopeNoise, slopeThreshold, ...
                            vMin, vMax, vRange, slopeMin, slopeMax, ...
                            detectStartMs, figTitleBaseThis);
        
    % Save the figure zoomed to several x limits
    save_all_zooms(fig, figPathBaseThis, stimStartMs, detectStartMs, histLeftMs);

    %% Plot the spike histogram
    edgesSeconds = edgesMs / MS_PER_S;
    xLimitsHist = edgesSeconds(1) + [0, 7];
    histFig = figure;
    [histBars, histFig] = ...
        plot_histogram([], 'Counts', spikeCounts, 'Edges', edgesSeconds, ...
                        'XLimits', xLimitsHist, 'XLabel', 'Time (seconds)', ...
                        'YLabel', 'Spike Count per 10 ms', ...
                        'FigTitle', ['Spike histogram for ', figTitleBaseThis], ...
                        'FigHandle', histFig);
    saveas(histFig, [figPathBaseThis, '_spike_histogram'], 'png');

    %% Plot the autocorrelogram
    % Create time values 
    tAcfTemp = create_time_vectors(nBins - 1, 'SamplingIntervalMs', binWidthMs, ...
                                'TimeUnits', 's');
    tAcf = [0; tAcfTemp];
    tAutoCorr = [-flipud(tAcfTemp); 0; tAcfTemp];
    timePeak1Sec = 0;
    timeTrough1Sec = idxTrough1 * binWidthSec;
    timePeak2Sec = idxPeak2 * binWidthSec;
    yLimits = compute_axis_limits(acf(1:floor(7/binWidthSec)), 'y', 'Coverage', 90);

    % Plot the autocorrelogram
    acfFig = figure;
    acfLine = plot(tAutoCorr, autoCorr);
    xlim([-7, 7]);
    ylim(yLimits);
    xlabel('Lag (s)');
    ylabel('Spike rate squared (Hz^2)');
    title(['Autocorrelation for ', figTitleBaseThis]);
    saveas(acfFig, [figPathBaseThis, '_autocorrelogram'], 'png');

    %% Plot the filtered autocorrelogram
    acfFilteredFig = figure;
    hold on;
    acfLine = plot(tAcf, acf, 'k');
    acfFilteredLine = plot(tAcf, acfFiltered, 'g', 'LineWidth', 1);
    plot(timePeak1Sec, ampPeak1, 'ro', 'LineWidth', 2);
    plot(timeTrough1Sec, ampTrough1, 'bx', 'LineWidth', 2);
    plot(timePeak2Sec, ampPeak2, 'ro', 'LineWidth', 2);
    text(0.5, 0.95, sprintf('Oscillatory Index = %g', oscIndex), ...
        'Units', 'normalized');
    text(0.5, 0.9, sprintf('Oscillation Period = %g ms', oscPeriodMs), ...
        'Units', 'normalized');
    xlim([0, 7]);
    ylim(yLimits);
    xlabel('Lag (s)');
    ylabel('Spike rate squared (Hz^2)');
    title(['Smoothed autocorrelation for ', figTitleBaseThis]);
    saveas(acfFilteredFig, [figPathBaseThis, '_smoothed_autocorrelogram'], 'png');
else
    fig = gobjects(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig, ax, lines, markers, raster] = ...
                plot_spike_detection(tVec, vVec, slopes, idxSpikes, ...
                                    baseSlopeNoise, slopeThreshold, ...
                                    vMin, vMax, vRange, slopeMin, slopeMax, ...
                                    detectStartMs, figTitle)
%% Plots the spike detection

% Hard-coded constants
barWidth2Range = 1/10;

% Compute the midpoint and bar width for the raster
barWidth = vRange * barWidth2Range;
yMid = vMax + barWidth;

% Compute y axis limits
yLimits1 = compute_axis_limits([slopeMin, slopeMax], 'y', 'Coverage', 100);
yLimits2 = compute_axis_limits([vMin, vMax], 'y', 'Coverage', 100);
yLimits3 = compute_axis_limits([vMin, yMid], 'y', 'Coverage', 100);

% Initialize graphics object handles
ax = gobjects(3, 1);
lines = gobjects(5, 1);
markers = gobjects(2, 1);

% Make a figure for spike detection
fig = figure(1); 
clf; 

% Plot the slope trace
ax(1) = subplot(3, 1, 1);
cla; hold on
lines(1) = plot(tVec(1:(end-1)), slopes, 'k');
lines(4) = plot_horizontal_line(baseSlopeNoise, 'Color', 'b', 'LineStyle', '--');
lines(5) = plot_horizontal_line(slopeThreshold, 'Color', 'g', 'LineStyle', '--');
markers(1) = plot(tVec(idxSpikes - 1), slopes(idxSpikes - 1), 'rx', 'LineWidth', 2);
ylim(yLimits1);
ylabel('Slope (V/s)');
title('Detection of peaks in the slope vector');

% Plot the original trace
ax(2) = subplot(3, 1, 2);
cla; hold on
lines(2) = plot(tVec, vVec, 'k');
markers(2) = plot(tVec(idxSpikes), vVec(idxSpikes), 'rx', 'LineWidth', 2);
ylim(yLimits2);
ylabel('Voltage (mV)');
title('Corresponding positions in the voltage vector');

% Plot the original trace
ax(3) = subplot(3, 1, 3);
cla; hold on
lines(3) = plot(tVec, vVec, 'k');
raster = plot_raster(tVec(idxSpikes), 'YMid', yMid, 'BarWidth', barWidth, ...
                    'LineWidth', 0.5, 'Colors', {'Red'}, ...
                    'YLimits', 'suppress', 'YTickLocs', 'suppress', ...
                    'YTickLabels', 'suppress');
ylim(yLimits3);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Original voltage vector with spikes');

% Create an overarching title
suptitle(figTitle);

% Link the x axes
linkaxes(ax, 'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_all_zooms(fig, figPathBase, stimStartMs, detectStartMs, firstSpikeMs)
%% Save the figure as .fig and 4 zooms as .png
% TODO: Make this more general

% Get the figure
figure(fig)

% Save the full figure
%save_all_figtypes(fig, [figPathBase, '_full'], {'png', 'fig'});
saveas(fig, [figPathBase, '_full'], 'png');

% Zoom #1
xlim([stimStartMs, stimStartMs + 1e4]);
saveas(fig, [figPathBase, '_zoom1'], 'png');

% Zoom #2
xlim([detectStartMs, detectStartMs + 2e3]);
saveas(fig, [figPathBase, '_zoom2'], 'png');

% Zoom #3
xlim([firstSpikeMs, firstSpikeMs + 60]);
saveas(fig, [figPathBase, '_zoom3'], 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute baseline rms noise from window
baseNoises = compute_baseline_noise(vVecs, tVec, baseWindow);

% Compute a baseline slope noise in V/s
baseSlopeNoise = baseNoise / siMs;

parsedParams.baseNoise = baseNoise;

idxDetectStart = find(tVec > detectStartMs, 1);

xlim([detectStartMs, detectStartMs + 1e4]);
xlim([detectStartMs, detectStartMs + 2e3]);
xlim([3410, 3470]);

% Query the maximum and range of vVec after detectStartMs
vVecTrunc = vVec(idxDetectStart:end);
vMean = mean(vVecTrunc);
vStd = std(vVecTrunc);
vMin = vMean - 10 * vStd;
vMax = vMean + 10 * vStd;
vRange = vMax - vMin;

% Query the maximum and range of slope after detectStartMs
slopesTrunc = slopes(idxDetectStart:end);
slopesMean = mean(slopesTrunc);
slopesStd = std(slopesTrunc);
slopeMin = slopesMean - 10 * slopesStd;
slopeMax = slopesMean + 10 * slopesStd;
slopeRange = slopeMax - slopeMin;

filterCutoffHz = 3;
acfFiltered = freqfilter(acf, filterCutoffHz, binWidthMs / 1000);
parsedParams.filterCutoffHz = filterCutoffHz;


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%