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
%       cd/extract_elements.m
%       cd/force_column_cell.m
%       cd/compute_axis_limits.m
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/match_time_info.m
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
%parfor iVec = 1:nVectors
for iVec = 1:1
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
    spikeTimes = parsedData.spikeTimes;
    stimStartMs = parsedParams.stimStartMs;
    detectStartMs = parsedParams.detectStartMs;

    % Create figure and plot
    figs(nVectors + 1) = figure(2);
    clf
    [hLines, eventTimes, yEnds, yTicksTable] = ...
        plot_raster(spikeTimes, 'LineWidth', 0.5);
    vLines = plot_vertical_line(mean(stimStartMs), 'Color', 'g', 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Trace #');
    title(['Spike times for ', figTitleBase]);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(nVectors + 1), figPathBaseThis, ...
                    mean(stimStartMs), mean(detectStartMs));
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

% Hard-coded parameters
signal2Noise = 4; %3
minDelaySamples = 2000;

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

% Store spike times
spikeTimes = tVec(idxSpikes);

% Query the maximum and range of vVec after detectStartMs
vMin = min(vVec(idxDetectStart:end));
vMax = max(vVec(idxDetectStart:end));
vRange = vMax - vMin;

% Query the maximum and range of slope after detectStartMs
slopeMin = min(slope(idxDetectStart:end));
slopeMax = max(slope(idxDetectStart:end));
slopeRange = slopeMax - slopeMin;

% Store in outputs
parsedParams.signal2Noise = signal2Noise;
parsedParams.minDelaySamples = minDelaySamples;
parsedParams.siMs = siMs;
parsedParams.idxStimStart = idxStimStart;
parsedParams.stimStartMs = stimStartMs;
parsedParams.baseWindow = baseWindow;
parsedParams.baseSlopeNoise = baseSlopeNoise;
parsedParams.slopeThreshold = slopeThreshold;
parsedParams.idxDetectStart = idxDetectStart;
parsedParams.detectStartMs = detectStartMs;
parsedParams.vMin = vMin;
parsedParams.vMax = vMax;
parsedParams.vRange = vRange;
parsedParams.slopeMin = slopeMin;
parsedParams.slopeMax = slopeMax;
parsedParams.slopeRange = slopeRange;
parsedData.tVec = tVec;
parsedData.vVec = vVec;
parsedData.slopes = slopes;
parsedData.idxSpikes = idxSpikes;
parsedData.spikeTimes = spikeTimes;

% Plots the spike detection
if plotFlag && iVec == 1
    % Plot spike detection
    [fig, ax, lines, markers, raster] = ...
        plot_spike_detection(tVec, vVec, slopes, idxSpikes, ...
                            baseSlopeNoise, slopeThreshold, ...
                            vMin, vMax, vRange, slopeMin, slopeMax, ...
                            detectStartMs, figTitleBaseThis);
        
    % Save the figure zoomed to several x limits
    save_all_zooms(fig, figPathBaseThis, stimStartMs, detectStartMs);
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

function save_all_zooms(fig, figPathBase, stimStartMs, detectStartMs)
%% Save the figure as .fig and 4 zooms as .png
% TODO: Make this more general

% Get the figure
figure(fig)

% Save the full figure
save_all_figtypes(fig, [figPathBase, '_full'], {'png', 'fig'});

% Zoom #1
xlim([stimStartMs, stimStartMs + 1e4]);
saveas(fig, [figPathBase, '_zoom1'], 'png');

% Zoom #2
xlim([detectStartMs, detectStartMs + 2e3]);
saveas(fig, [figPathBase, '_zoom2'], 'png');

% Zoom #3
xlim([3410, 3470]);
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%