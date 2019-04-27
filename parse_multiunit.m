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
%                   - 'PhaseBoundaries': vector of phase boundaries
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   
% Requires:
%       TODO cd/find_stim_start.m
%       cd/argfun.m
%       cd/check_dir.m
%       cd/check_subdir.m
%       cd/compute_axis_limits.m
%       cd/compute_baseline_noise.m
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_logical_array.m
%       cd/create_time_vectors.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_nearest_multiple.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%       cd/plot_horizontal_line.m
%       cd/plot_raster.m
%       cd/plot_table.m
%       cd/transform_vectors.m
%
% Used by:
%       cd/parse_all_multiunit.m

% File History:
% 2019-02-19 Created by Adam Lu
% 2019-02-24 Added computation of oscillation index and period
% 2019-02-25 Added computation of oscillation duration
% 2019-02-26 Updated computation of oscillation duration
% 2019-02-26 Updated computation of oscillation index
% 2019-03-14 Nows places figures in subdirectories
% 2019-03-14 Nows computes appropriate x and y limits for all traces
% 2019-03-14 Fixed plotting of oscillation duration in histograms
% 2019-03-14 Redefined the oscillation period so that it is between the primary
%               peak and the next largest-amplitude peak
% 2019-03-14 Redefined the oscillatory index so that it is the reciprocal of 
%               the coefficient of variation of the lag differences 
%               between consecutive peaks
% 2019-03-15 Redefined the oscillatory index so that it is 
%               1 minus the average of all distances 
%               (normalized by half the oscillation period)
%               to the closest multiple of the period over all peaks
% 2019-03-17 Added nSpikesPerBurstInOsc, nSpikesInOsc, nBurstsInOsc, etc ...
% 2019-03-19 Added nSpikesPerBurstIn10s, nSpikesIn10s, nBurstsIn10s, etc ...
% 2019-03-24 Fixed bugs in prepare_for_plot_horizontal_line.m
% 2019-03-24 Renamed setNumber -> phaseNumber, setName -> phaseName

% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
rawDir = 'raw';
rasterDir = 'rasters';
autoCorrDir = 'autocorrelograms';
acfDir = 'autocorrelation_functions';
spikeHistDir = 'spike_histograms';
spikeDetectionDir = 'spike_detections';
measuresDir = 'measures';
measuresToPlot = {'oscIndex1', 'oscIndex2', 'oscIndex3', 'oscIndex4', ...
                    'oscPeriod1Ms', 'oscPeriod2Ms', ...
                    'oscDurationSec', ...
                    'nSpikesTotal', 'nSpikesIn10s', 'nSpikesInOsc', ...
                    'nBurstsTotal', 'nBurstsIn10s', 'nBurstsInOsc', ...
                    'nSpikesPerBurst', 'nSpikesPerBurstIn10s', ...
                    'nSpikesPerBurstInOsc'};

%% Default values for optional arguments
plotFlagDefault = false;
outFolderDefault = pwd;
fileBaseDefault = {};           % set later
stimStartMsDefault = [];        % set later
pulseVectorsDefault = [];       % don't use pulse vectors by default
setBoundariesDefault = [];      % no set boundaries by default
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
addParameter(iP, 'PhaseBoundaries', setBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
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
phaseBoundaries = iP.Results.PhaseBoundaries;
tVecs = iP.Results.tVecs;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vVecs);

% Count the number of samples for each vector
nSamples = count_samples(vVecs);

% Match time vector(s) with sampling interval(s) and number(s) of samples
[tVecs, siMs, nSamples] = match_time_info(tVecs, siMs, nSamples);

% Count the number of measures to plot
nMeasures = numel(measuresToPlot);

% Initialize figures array
figs = gobjects(nMeasures + 2, 1);

% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

%% Do the job
% Detect stimulation start time if not provided
%   Otherwise find the corresponding index in the time vector
fprintf('Detecting stimulation start for %s ...\n', fileBase);
if isempty(stimStartMs)
    % TODO: Make this a function find_stim_start.m
    if ~isempty(pulseVectors)
        % Parse the pulse vectors
        [pulseParams, ~] = ...
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
fprintf('Constructing baseline window for %s ...\n', fileBase);
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
fprintf('Parsing recording for %s ...\n', fileBase);
parsedParamsCell = cell(nVectors, 1);
parsedDataCell = cell(nVectors, 1);
parfor iVec = 1:nVectors
%for iVec = 1:nVectors
%for iVec = 1:1
    [parsedParamsCell{iVec}, parsedDataCell{iVec}] = ...
        parse_multiunit_helper(iVec, vVecs{iVec}, tVecs{iVec}, siMs(iVec), ...
                                idxStimStart(iVec), stimStartMs(iVec), ...
                                baseWindows{iVec}, ...
                                fileBase, titleBase, phaseBoundaries);
end

% Convert to a struct array
%   Note: This removes all entries that are empty
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@(x) [x{:}], parsedParamsCell, parsedDataCell);

% Convert to a table
[parsedParams, parsedData] = ...
    argfun(@(x) struct2table(x, 'AsArray', true), ...
            parsedParamsStruct, parsedDataStruct);

% Save the parameters table
writetable(parsedParams, fullfile(outFolder, [fileBase, '_params.csv']));

%% Plot spike detection
if plotFlag
    fprintf('Plotting spike detection for %s ...\n', fileBase);

    % Retrieve data for plotting
    tVec = parsedData.tVec;
    vVec = parsedData.vVec;
    slopes = parsedData.slopes;
    idxSpikes = parsedData.idxSpikes;

    stimStartMs = parsedParams.stimStartMs;
    detectStartMs = parsedParams.detectStartMs;
    firstSpikeMs = parsedParams.firstSpikeMs;
    baseSlopeNoise = parsedParams.baseSlopeNoise;
    slopeThreshold = parsedParams.slopeThreshold;
    vMin = parsedParams.vMin;
    vMax = parsedParams.vMax;
    vRange = parsedParams.vRange;
    slopeMin = parsedParams.slopeMin;
    slopeMax = parsedParams.slopeMax;
    figTitleBase = parsedParams.figTitleBase;
    figPathBase = parsedParams.figPathBase;

    % Create output directory
    outFolderSpikeDetection = fullfile(outFolder, spikeDetectionDir);
    check_dir(outFolderSpikeDetection);

    parfor iVec = 1:nVectors
        % Plot spike detection
        [fig, ax, lines, markers, raster] = ...
            plot_spike_detection(tVec{iVec}, vVec{iVec}, ...
                                slopes{iVec}, idxSpikes{iVec}, ...
                                baseSlopeNoise(iVec), slopeThreshold(iVec), ...
                                vMin(iVec), vMax(iVec), vRange(iVec), ...
                                slopeMin(iVec), slopeMax(iVec), ...
                                [], figTitleBase{iVec});

        % Set zoom windows
        zoomWin1 = stimStartMs(iVec) + [0, 1e4];
        zoomWin2 = detectStartMs(iVec) + [0, 2e3];
        if ~isnan(firstSpikeMs(iVec))
            zoomWin3 = firstSpikeMs(iVec) + [0, 60];
        else
            zoomWin3 = [0, 60];
        end            

        % Save the figure zoomed to several x limits
        save_all_zooms(fig, outFolderSpikeDetection, ...
                        figPathBase{iVec}, zoomWin1, zoomWin2, zoomWin3);

        close all force hidden
    end
end

%% Plot spike histograms
if plotFlag
    fprintf('Plotting spike histograms for %s ...\n', fileBase);

    % Retrieve data for plotting
    spikeCounts = parsedData.spikeCounts;
    edgesSec = parsedData.edgesSec;
    timeBurstStartsSec = parsedData.timeBurstStartsSec;
    timeBurstEndsSec = parsedData.timeBurstEndsSec;

    binWidthSec = parsedParams.binWidthSec;
    histLeftSec = parsedParams.histLeftSec;
    timeOscEndSec = parsedParams.timeOscEndSec;
    maxInterBurstIntervalSec = parsedParams.maxInterBurstIntervalSec;
    oscDurationSec = parsedParams.oscDurationSec;
    nSpikesInOsc = parsedParams.nSpikesInOsc;
    figTitleBase = parsedParams.figTitleBase;
    figPathBase = parsedParams.figPathBase;

    % Find appropriate x limits
    histLeft = min(histLeftSec);
    % histRight = nanmean(timeOscEndSec) + 1.96 * stderr(timeOscEndSec) + ...
    %                 1.5 * max(maxInterBurstIntervalSec);
    histRight = 10;
    xLimitsHist = [histLeft, histRight];

    % Compute x limits for durations
%    durationWindows = force_column_cell(transpose([histLeft, timeOscEndSec]));
    durationWindows = cellfun(@(x, y) prepare_for_plot_horizontal_line(x, y), ...
                            timeBurstStartsSec, timeBurstEndsSec, ...
                            'UniformOutput', false);

    % Find the last bin to show for all traces
    lastBinToShow = floor((histRight - histLeft) ./ binWidthSec) + 1;
    
    % Find appropriate y limits
    spikeCountsOfInterest = extract_subvectors(spikeCounts, ...
                            'IndexEnd', lastBinToShow);
    largestSpikeCount = apply_iteratively(@max, spikeCountsOfInterest);
    yLimitsHist = [0, largestSpikeCount * 1.1];

    % Create output directory
    outFolderHist = fullfile(outFolder, spikeHistDir);
    check_dir(outFolderHist);

    % Plot histograms
    parfor iVec = 1:nVectors
        [histBars, histFig] = ...
            plot_spike_histogram(spikeCounts{iVec}, edgesSec{iVec}, ...
                                durationWindows{iVec}, oscDurationSec(iVec), ...
                                nSpikesInOsc(iVec), ...
                                xLimitsHist, yLimitsHist, figTitleBase{iVec});

        saveas(histFig, fullfile(outFolderHist, ...
                        [figPathBase{iVec}, '_spike_histogram']), 'png');
        close all force hidden
    end
end

%% Plot autocorrelograms
if plotFlag
    fprintf('Plotting autocorrelograms for %s ...\n', fileBase);

    % Retrieve data for plotting
    autoCorr = parsedData.autoCorr;
    acf = parsedData.acf;
    acfFiltered = parsedData.acfFiltered;
    indPeaks = parsedData.indPeaks;
    indTroughs = parsedData.indTroughs;
    ampPeaks = parsedData.ampPeaks;
    ampTroughs = parsedData.ampTroughs;

    binWidthSec = parsedParams.binWidthSec;
    nBins = parsedParams.nBins;
    halfNBins = parsedParams.halfNBins;
    oscIndex1 = parsedParams.oscIndex1;
    oscIndex2 = parsedParams.oscIndex2;
    oscIndex3 = parsedParams.oscIndex3;
    oscIndex4 = parsedParams.oscIndex4;
    oscPeriod2Ms = parsedParams.oscPeriod2Ms;
    oscPeriod1Ms = parsedParams.oscPeriod1Ms;
    oscDurationSec = parsedParams.oscDurationSec;
    nSpikesInOsc = parsedParams.nSpikesInOsc;
    figTitleBase = parsedParams.figTitleBase;
    figPathBase = parsedParams.figPathBase;

    % Compute appropriate x limits
    allLastPeaksBins = extract_elements(indPeaks, 'last');
    allLastPeaksSec = allLastPeaksBins .* binWidthSec;
    allOscDur = oscDurationSec;
    bestRightForAll = max([allOscDur, allLastPeaksSec], [], 2) + 1;
    acfFilteredRight = compute_stats(bestRightForAll, 'upper95', ...
                                    'RemoveOutliers', true);
    % xLimitsAutoCorr = [-acfFilteredRight, acfFilteredRight];
    xLimitsAutoCorr = [-10, 10];
    % xLimitsAcfFiltered = [0, acfFilteredRight];
    xLimitsAcfFiltered = [0, 10];

    % Find the last index to show
    lastIndexToShow = floor(acfFilteredRight ./ binWidthSec) + 1;
    
    % Compute appropriate y limits
    acfOfInterest = extract_subvectors(acf, 'IndexEnd', lastIndexToShow);
    largestAcfValues = extract_elements(acfOfInterest, 'max');
    bestUpperLimit = compute_stats(largestAcfValues, 'upper95', ...
                                    'RemoveOutliers', true);
    yLimitsAutoCorr = compute_axis_limits([0, bestUpperLimit], ...
                                            'y', 'Coverage', 95);
    yLimitsAcfFiltered = compute_axis_limits([0, bestUpperLimit], ...
                                            'y', 'Coverage', 90);
    yOscDur = -(bestUpperLimit * 0.025);

    % Create output directories
    outFolderAutoCorr = fullfile(outFolder, autoCorrDir);
    check_dir(outFolderAutoCorr);
    outFolderAcf = fullfile(outFolder, acfDir);
    check_dir(outFolderAcf);

    parfor iVec = 1:nVectors
        [autoCorrFig, acfFig] = ...
            plot_autocorrelogram(autoCorr{iVec}, acf{iVec}, acfFiltered{iVec}, ...
                indPeaks{iVec}, indTroughs{iVec}, ...
                ampPeaks{iVec}, ampTroughs{iVec}, ...
                binWidthSec(iVec), nBins(iVec), halfNBins(iVec), ...
                oscIndex1(iVec), oscIndex2(iVec), ...
                oscIndex3(iVec), oscIndex4(iVec), ...
                oscPeriod1Ms(iVec), oscPeriod2Ms(iVec), ...
                oscDurationSec(iVec), nSpikesInOsc(iVec), ...
                xLimitsAutoCorr, yLimitsAutoCorr, ...
                xLimitsAcfFiltered, yLimitsAcfFiltered, ...
                yOscDur, figTitleBase{iVec});

        saveas(autoCorrFig, fullfile(outFolderAutoCorr, ...
                [figPathBase{iVec}, '_autocorrelogram']), 'png');
        saveas(acfFig, fullfile(outFolderAcf, ...
                [figPathBase{iVec}, '_autocorrelation_function']), 'png');

        close all force hidden
    end
end

%% Plot raw traces
%if plotFlag
    fprintf('Plotting raw traces for %s ...\n', fileBase);

    % Modify the figure base
    figBaseRaw = [fileBase, '_raw'];

    % Create output directory
    outFolderRaw = fullfile(outFolder, rawDir);
    check_dir(outFolderRaw);

    % Extract parameters
    stimStartSec = parsedParams.stimStartSec;
    detectStartSec = parsedParams.detectStartSec;
    firstSpikeSec = parsedParams.firstSpikeSec;
    timeOscEndSec = parsedParams.timeOscEndSec;

    % Convert time vector to seconds
    tVecsSec = transform_vectors(tVecs, MS_PER_S, 'divide');

    % Prepare for the plot
    xLabel = 'Time (s)';
    figTitle = ['Raw traces for ', titleBase];

    % Compute the original y limits from data
    yLimitsOrig = compute_axis_limits(vVecs, 'y', 'AutoZoom', autoZoom);

    % Compute the amount of y to stagger
    yAmountToStagger = range(yLimitsOrig);

    % Create figure and plot
    figs(1) = figure('Visible', 'off');
    clf
    plot_traces(tVecsSec, vVecs, 'Verbose', false, ...
                'PlotMode', 'staggered', 'SubplotOrder', 'list', ...
                'YLimits', yLimitsOrig, 'YAmountToStagger', yAmountToStagger, ...
                'XLabel', xLabel, 'LinkAxesOption', 'y', ...
                'YLabel', 'suppress', 'TraceLabels', 'suppress', ...
                'FigTitle', figTitle, 'FigHandle', figs(1), ...
                'Color', 'k');

    vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                    'LineStyle', '--');
    if ~isempty(phaseBoundaries)
        yBoundaries = nVectors - phaseBoundaries + 1;
        horzLine = plot_horizontal_line(yBoundaries, 'Color', 'g', ...
                                        'LineStyle', '--', 'LineWidth', 2);
    end

    % Save the figure zoomed to several x limits
    zoomWin1 = mean(stimStartSec) + [0, 10];
    zoomWin2 = mean(detectStartSec) + [0, 2];
    meanFirstSpike = nanmean(firstSpikeSec);
    if ~isnan(meanFirstSpike)
        zoomWin3 = meanFirstSpike + [0, 0.06];
    else
        zoomWin3 = [0, 0.06];
    end            
    save_all_zooms(figs(1), outFolderRaw, ...
                    figBaseRaw, zoomWin1, zoomWin2, zoomWin3);
%end

%% Plot raster plot
% TODO: Plot burst duration
% TODO: Plot oscillatory index
if plotFlag
    fprintf('Plotting raster plot for %s ...\n', fileBase);

    % Modify the figure base
    figBaseRaster = [fileBase, '_raster'];

    % Create output directory
    outFolderRaster = fullfile(outFolder, rasterDir);
    check_dir(outFolderRaster);

    % Extract the spike times
    spikeTimesSec = parsedData.spikeTimesSec;
    timeBurstStartsSec = parsedData.timeBurstStartsSec;
    timeBurstEndsSec = parsedData.timeBurstEndsSec;

    stimStartSec = parsedParams.stimStartSec;
    detectStartSec = parsedParams.detectStartSec;
    firstSpikeSec = parsedParams.firstSpikeSec;
    timeOscEndSec = parsedParams.timeOscEndSec;

    % Convert oscillatory index to a window
    % TODO

    % Oscillation window
    oscWindow = transpose([stimStartSec, timeOscEndSec]);

    % Burst windows
    burstWindows = cellfun(@(x, y) prepare_for_plot_horizontal_line(x, y), ...
                            timeBurstStartsSec, timeBurstEndsSec, ...
                            'UniformOutput', false);

    % Create colors
    nSweeps = numel(spikeTimesSec);
    colorsRaster = repmat({'Black'}, nSweeps, 1);

    % Create figure and plot
    figs(2) = figure('Visible', 'off');
    clf
    [hLines, eventTimes, yEnds, yTicksTable] = ...
        plot_raster(spikeTimesSec, 'DurationWindow', burstWindows, ...
                    'LineWidth', 0.5, 'Colors', colorsRaster);
    % [hLines, eventTimes, yEnds, yTicksTable] = ...
    %     plot_raster(spikeTimesSec, 'DurationWindow', oscWindow, ...
    %                 'LineWidth', 0.5, 'Colors', colorsRaster);
    vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                    'LineStyle', '--');
    if ~isempty(phaseBoundaries)
        yBoundaries = nVectors - phaseBoundaries + 1;
        horzLine = plot_horizontal_line(yBoundaries, 'Color', 'g', ...
                                        'LineStyle', '--', 'LineWidth', 2);
    end
    xlabel('Time (s)');
    ylabel('Trace #');
    title(['Spike times for ', titleBase]);

    % Save the figure zoomed to several x limits
    zoomWin1 = mean(stimStartSec) + [0, 10];
    zoomWin2 = mean(detectStartSec) + [0, 2];
    meanFirstSpike = nanmean(firstSpikeSec);
    if ~isnan(meanFirstSpike)
        zoomWin3 = meanFirstSpike + [0, 0.06];
    else
        zoomWin3 = [0, 0.06];
    end            
    save_all_zooms(figs(2), outFolderRaster, ...
                    figBaseRaster, zoomWin1, zoomWin2, zoomWin3);
end

%% Plot time series of measures
if plotFlag
    fprintf('Plotting time series of measures for %s ...\n', fileBase);    

    % Create output directory and subdirectories for each measure
    outFolderMeasures = fullfile(outFolder, measuresDir);
    check_dir(outFolderMeasures);
    check_subdir(outFolderMeasures, measuresToPlot);

    % Create full figure paths
    figPathsMeasures = fullfile(outFolderMeasures, measuresToPlot, ...
                                strcat(fileBase, '_', measuresToPlot));

    % Create custom figure titles
    figTitlesMeasures = [measuresToPlot, ' for ', titleBase];

    % Plot table
    figs(3:(nMeasures + 2)) = ...
        plot_table(parsedParams, 'VariableNames', measuresToPlot, ...
                    'XLabel', 'Time (min)', 'FigNames', figPathsMeasures, ...
                    'FigTitles', figTitlesMeasures, ...
                    'XBoundaries', phaseBoundaries, ...
                    'RemoveOutliers', true, 'PlotSeparately', true);
end

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;
varargout{3} = figs;

fprintf('%s analyzed! ...\n\n', fileBase);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
                parse_multiunit_helper(iVec, vVec, tVec, siMs, ...
                                idxStimStart, stimStartMs, baseWindow, ...
                                fileBase, figTitleBase, phaseBoundaries)

% Parse a single multiunit recording

% Hard-coded constants
MS_PER_S = 1000;

% Hard-coded parameters
signal2Noise = 3; %4
minDelayMs = 25;
binWidthMs = 10;
filterWidthMs = 100;
minRelProm = 0.02;
minSpikeRateInBurstHz = 100;
minBurstLengthMs = 20;
maxInterBurstIntervalMs = 1000; %2000;

%% Preparation
% Compute the minimum delay in samples
minDelaySamples = ceil(minDelayMs ./ siMs);

% Compute the bin width in seconds
binWidthSec = binWidthMs ./ MS_PER_S;

% Compute the number of phase boundaries
nBoundaries = numel(phaseBoundaries);

% Determine which phase number this sweep belongs to
if nBoundaries > 0
    % For the first n - 1 sets, use find
    phaseNumber = find(phaseBoundaries > iVec, 1, 'first');

    % For the last phase, use numel(phaseBoundaries) + 1
    if isempty(phaseNumber)
        phaseNumber = numel(phaseBoundaries) + 1;
    end
else
    phaseNumber = NaN;
end

% Create phase names
if nBoundaries == 2
    if phaseNumber == 1
        phaseName = 'baseline';
    elseif phaseNumber == 2
        phaseName = 'washon';
    elseif phaseNumber == 3
        phaseName = 'washoff';
    end
else
    phaseName = '';
end

%% Detect spikes
% Find the starting index for detecting a spike
idxDetectStart = idxStimStart + minDelaySamples;

% Find the corresponding time
detectStartMs = tVec(idxDetectStart);

% Compute the number of samples
nSamples = numel(vVec);

% Compute all instantaneous slopes in V/s
slopes = diff(vVec) ./ siMs;

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
nSpikesTotal = numel(idxSpikes);

% Index of first spike
if nSpikesTotal == 0
    idxFirstSpike = NaN;
else
    idxFirstSpike = idxSpikes(1);
end

% Store spike times
if nSpikesTotal == 0
    spikeTimesMs = [];
    firstSpikeMs = NaN;
else
    spikeTimesMs = tVec(idxSpikes);
    firstSpikeMs = spikeTimesMs(1);    
end

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

%% Compute the spike histogram, spikes per burst & oscillation duration
if nSpikesTotal == 0
    spikeCounts = [];
    edgesMs = [];
    nBins = 0;
    halfNBins = 0;
    histLeftMs = NaN;
    nBurstsTotal = 0;
    nBurstsIn10s = 0;
    nBurstsInOsc = 0;
    iBinBurstStarts = [];
    iBinBurstEnds = [];
    iBinBurstIn10sStarts = [];
    iBinBurstIn10sEnds = [];
    iBinBurstInOscStarts = [];
    iBinBurstInOscEnds = [];
    iBinLastOfLastBurst = NaN;
    iBinLastOfLastBurstIn10s = NaN;
    iBinLastOfLastBurstInOsc = NaN;
    spikeCountsEachBurst = [];
    spikeCountsEachBurstIn10s = [];
    spikeCountsEachBurstInOsc = [];
    nSpikesPerBurst = NaN;
    nSpikesPerBurstIn10s = NaN;
    nSpikesPerBurstInOsc = NaN;
    nSpikesIn10s = 0;
    nSpikesInOsc = 0;
    timeBurstStartsMs = [];
    timeBurstEndsMs = [];
    timeBurstIn10sStartsMs = [];
    timeBurstIn10sEndsMs = [];
    timeBurstInOscStartsMs = [];
    timeBurstInOscEndsMs = [];
    timeOscEndMs = stimStartMs;
    oscDurationMs = 0;    
else
    % Compute a spike histogram
    [spikeCounts, edgesMs] = compute_bins(spikeTimesMs, 'BinWidth', binWidthMs);

    % Compute the number of bins
    nBins = numel(spikeCounts);

    % Compute the half number of bins
    halfNBins = floor(nBins/2);

    % Record the starting time of the histogram
    histLeftMs = edgesMs(1);

    % Compute the minimum number of bins in a burst
    minBinsInBurst = ceil(minBurstLengthMs ./ binWidthMs);

    % Compute the sliding window length in seconds
    slidingWinSec = minBinsInBurst * binWidthSec;

    % Compute the minimum spikes per sliding window if in a burst
    minSpikesPerWindowInBurst = ceil(minSpikeRateInBurstHz * slidingWinSec);

    % Compute the maximum number of bins between consecutive bursts
    maxIbiBins = floor(maxInterBurstIntervalMs ./ binWidthMs);

    % Count spikes for each sliding window (ending at each bin)
    spikeCountsWin = spikeCounts;
    spikeCountsPrev = spikeCounts;
    for i = 1:(minBinsInBurst - 1)
        % Compute spike counts from the previous ith bin
        spikeCountsPrev = [0; spikeCountsPrev(1:(end-1))];

        % Add the spike counts in the previous i bins
        spikeCountsWin = spikeCountsWin + spikeCountsPrev;
    end

    % Determine whether each sliding window passes 
    %   the number of spikes criterion
    isInBurst = spikeCountsWin >= minSpikesPerWindowInBurst;

    % Find the last bins of each burst
    %   Note: this sliding window is in burst but the next sliding window 
    %           is not in burst
    iWinLastInBurst = find(isInBurst & [~isInBurst(2:end); true]);
    iBinBurstEnds = iWinLastInBurst;

    % Find the first bins of each burst
    %   Note: this sliding window is in burst but the previous sliding window 
    %           is not in burst, then minus the number of bins in a window
    iWinFirstInBurst = find(isInBurst & [true; ~isInBurst(1:(end-1))]);
    iBinBurstStarts = iWinFirstInBurst - minBinsInBurst + 1;
    iBinBurstStarts(iBinBurstStarts < 1) = 1;
    
    % Find the total number of bursts
    if numel(iBinBurstEnds) ~= numel(iBinBurstStarts)
        error('Code logic error!');
    else
        nBurstsTotal = numel(iBinBurstEnds);
    end

    % Count bursts within the first 10 seconds after stimulation start
    if isempty(iBinBurstEnds)
        nSpikesIn10s = 0;
        nBurstsIn10s = 0;
        iBinBurstIn10sStarts = [];
        iBinBurstIn10sEnds = [];
    else
        timeAfterHistLeftFor10sMs = 10000 - (histLeftMs - stimStartMs);

        % Find the last bin 10 secs after stimulation start
        %   Note: cannot exceed nBins
        iBin10s = min(nBins, ceil(timeAfterHistLeftFor10sMs ./ binWidthMs) + 1);

        % Count the number of spikes within 10 secs after stimulation start
        nSpikesIn10s = sum(spikeCounts(1:iBin10s));

        % Count the number of bursts within 10 secs after stimulation start
        nBurstsIn10s = find(iBinBurstEnds <= iBin10s, 1, 'last');

        % Restrict the starting and ending bins vectors to ten seconds
        %   after stimulation start
        if isempty(nBurstsIn10s)
            nBurstsIn10s = 0;
            iBinBurstIn10sStarts = [];
            iBinBurstIn10sEnds = [];
        else
            iBinBurstIn10sStarts = iBinBurstStarts(1:nBurstsIn10s);
            iBinBurstIn10sEnds = iBinBurstEnds(1:nBurstsIn10s);
        end
    end

    % Detect the oscillation end and count bursts within oscillation
    %   Note: Find the last bin of the last burst, using maxIbiBins
    if isempty(iBinBurstEnds)
        nBurstsInOsc = 0;
        iBinBurstInOscStarts = [];
        iBinBurstInOscEnds = [];
    else
        % Compute the inter-burst intervals in bins
        ibiBins = diff(iBinBurstEnds);

        % Find the first inter-burst interval greater than maxIbiBins
        firstIbiTooLarge = find(ibiBins > maxIbiBins, 1, 'first');

        % Determine the number of bursts in oscillation
        if isempty(firstIbiTooLarge)
            % All bursts are close enough together
            nBurstsInOsc = nBurstsTotal;
        else
            % Actually (firstIbiTooLarge - 1) + 1
            nBurstsInOsc = firstIbiTooLarge;
        end

        % Restrict the starting and ending bins vectors
        iBinBurstInOscStarts = iBinBurstStarts(1:nBurstsInOsc);
        iBinBurstInOscEnds = iBinBurstEnds(1:nBurstsInOsc);
    end

    % Compute burst statistics
    [iBinLastOfLastBurst, spikeCountsEachBurst, nSpikesPerBurst] = ...
            compute_burst_statistics(spikeCounts, ...
                                        iBinBurstStarts, iBinBurstEnds);
    [iBinLastOfLastBurstIn10s, spikeCountsEachBurstIn10s, ...
        nSpikesPerBurstIn10s] = ...
            compute_burst_statistics(spikeCounts, ...
                                iBinBurstIn10sStarts, iBinBurstIn10sEnds);
    [iBinLastOfLastBurstInOsc, spikeCountsEachBurstInOsc, ...
        nSpikesPerBurstInOsc] = ...
            compute_burst_statistics(spikeCounts, ...
                                iBinBurstInOscStarts, iBinBurstInOscEnds);

    % Compute the total number of spikes in the oscillation
    if ~isnan(iBinLastOfLastBurstInOsc)
        nSpikesInOsc = sum(spikeCounts(1:iBinLastOfLastBurstInOsc));
    else
        nSpikesInOsc = 0;
    end

    % Compute the burst start and ends in ms
    timeBurstStartsMs = histLeftMs + (iBinBurstStarts - 1) * binWidthMs;
    timeBurstEndsMs = histLeftMs + iBinBurstEnds * binWidthMs;
    timeBurstIn10sStartsMs = histLeftMs + (iBinBurstIn10sStarts - 1) * binWidthMs;
    timeBurstIn10sEndsMs = histLeftMs + iBinBurstIn10sEnds * binWidthMs;
    timeBurstInOscStartsMs = histLeftMs + (iBinBurstInOscStarts - 1) * binWidthMs;
    timeBurstInOscEndsMs = histLeftMs + iBinBurstInOscEnds * binWidthMs;

    % Find the time of oscillation end in ms
    if isnan(iBinLastOfLastBurstInOsc)
        timeOscEndMs = stimStartMs;
    else
        % Compute the time of oscillation end in ms
        timeOscEndMs = histLeftMs + iBinLastOfLastBurstInOsc * binWidthMs;
    end

    % Compute the oscillation duration in ms
    oscDurationMs = timeOscEndMs - stimStartMs;
end

%% Compute the autocorrelogram, oscillation period & oscillatory index
if nSpikesTotal == 0
    oscIndex1 = 0;
    oscIndex2 = 0;
    oscIndex3 = NaN;
    oscIndex4 = 0;
    oscPeriod1Ms = 0;
    oscPeriod2Ms = 0;
    minOscPeriod2Bins = 0;
    maxOscPeriod2Bins = 0;
    autoCorr = [];
    acf = [];
    acfFiltered = [];
    acfFilteredOfInterest = [];
    indPeaks = [];
    indTroughs = [];
    ampPeaks = [];
    ampTroughs = [];
    halfPeriodsToMultiple = [];
else
    % Compute an unnormalized autocorrelogram in Hz^2
    autoCorr = xcorr(spikeCounts, 'unbiased') ./ binWidthSec ^ 2;

    % Take just half of the positive side to get the autocorrelation function
    acf = autoCorr(nBins:(nBins + halfNBins));

    % Compute a normalized autocorrelation function
    % autocorr(spikeCounts, nBins - 1);
    % acf = autocorr(spikeCounts, nBins - 1);

    % Smooth the autocorrelation function with a moving-average filter
    acfFiltered = movingaveragefilter(acf, filterWidthMs, binWidthMs);

    % Record the amplitude of the primary peak
    ampPeak1 = acfFiltered(1);

    % Compute the oscillation duration in bins
    oscDurationBins = floor(oscDurationMs ./ binWidthMs);

    % Find the maximum bin of interest
    maxBinOfInterest = min(1 + oscDurationBins, numel(acfFiltered));
    
    % Restrict the autocorrelation function to oscillation duration
    acfFilteredOfInterest = acfFiltered(1:maxBinOfInterest);

    % Find the index and amplitude of peaks within oscillation duration
    if numel(acfFilteredOfInterest) > 3
        [peakAmp, peakInd] = ...
            findpeaks(acfFilteredOfInterest, ...
                        'MinPeakProminence', minRelProm * ampPeak1);

        % Record all peak indices and amplitudes
        indPeaks = [1; peakInd];
        ampPeaks = [ampPeak1; peakAmp];
    else
        indPeaks = 1;
        ampPeaks = ampPeak1;
    end

    % Compute the number of peaks
    nPeaks = numel(indPeaks);

    % Compute the lags of peaks in bins
    if nPeaks <= 1
        lagsPeaksBins = [];
    else
        lagsPeaksBins = indPeaks(2:end) - indPeaks(1);
    end

    % Compute the oscillation period
    %   Note: TODO
    % TODO: Try both:
    %   1. Use fminsearch on the distance to multiples function
    %   2. Use the largest peak in the frequency spectrum
    if nPeaks <= 1
        oscPeriod1Bins = 0;
        oscPeriod2Bins = 0;
        minOscPeriod2Bins = 0;
        maxOscPeriod2Bins = 0;
    else
        % Compute the oscillation period version 1
        %   Note: The lag between the primary peak and the second peak
        oscPeriod1Bins = indPeaks(2) - indPeaks(1);

        % Compute the oscillation period version 2
        %   Note: TODO
        if nPeaks <= 2
            % Just use the lag between first two peaks
            oscPeriod2Bins = oscPeriod1Bins;
            minOscPeriod2Bins = oscPeriod1Bins;
            maxOscPeriod2Bins = oscPeriod1Bins;
        else
            % Create a function for computing the average distance
            %   of each peak lag to a multiple of x
            average_distance_for_a_period = ...
                @(x) compute_average_distance_to_a_multiple(lagsPeaksBins, x);

            % Define the range the actual oscillation period can lie in
            minOscPeriod2Bins = oscPeriod1Bins * 2 / 3;
            maxOscPeriod2Bins = oscPeriod1Bins * 3 / 2;

            % Find the oscillation period in bins by looking for the 
            %   value of x that minimizes myFun, using the range
            %   oscPeriod1Bins / 3 to oscPeriod1Bins * 3
            oscPeriod2Bins = ...
                fminbnd(average_distance_for_a_period, ...
                        minOscPeriod2Bins, maxOscPeriod2Bins);
        end
    end

    % Convert to ms
    oscPeriod1Ms = oscPeriod1Bins .* binWidthMs;
    oscPeriod2Ms = oscPeriod2Bins .* binWidthMs;

    % Find the troughs
    if numel(indPeaks) <= 1
        indTroughs = [];
        ampTroughs = [];
    else
        % Find the indices and amplitudes of the troughs in between 
        %   each pair of peaks
        [ampTroughs, indTroughs] = ...
            find_troughs_from_peaks(acfFilteredOfInterest, indPeaks);
    end

    % Compute the oscillatory index version 1
    %   Note: This is defined in Sohal's paper 
    %           as the ratio of the difference between 
    %           the average of first two peaks and the first trough
    %           and the average of first two peaks
    if numel(indPeaks) <= 1
        oscIndex1 = 0;
    else
        % Compute the average amplitude of the first two peaks
        ampAvgFirstTwoPeaks = mean(ampPeaks(1:2));

        % Compute the oscillatory index
        oscIndex1 = (ampAvgFirstTwoPeaks - ampTroughs(1)) / ampAvgFirstTwoPeaks;
    end

    % Compute the oscillatory index version 2
    %   Note: This is the average of all oscillatory indices as defined
    %           by Sohal's paper between adjacent peaks
    if numel(indPeaks) <= 1
        oscIndex2 = 0;
    else
        % Compute the average amplitudes between adjacent peaks
        ampAdjPeaks = mean([ampPeaks(1:(end-1)), ampPeaks(2:end)], 2);

        % Take the average of oscillatory indices as defined by Sohal's paper
        oscIndex2 = mean((ampAdjPeaks - ampTroughs) ./ ampAdjPeaks);
    end

    % Compute the oscillatory index version 3
    %   Note: This is 1 minus the average of all distances 
    %           (normalized by half the oscillation period)
    %           to the closest multiple of the period over all peaks from the
    %           2nd and beyond. However, if there are less than two peaks,
    %           consider it non-oscillatory
    if numel(indPeaks) <= 2
        % If there are no more than two peaks, don't compute
        halfPeriodsToMultiple = [];
        oscIndex3 = NaN;
    else
        % Compute the distance to the closest multiple of the oscillation period 
        %   for each peak from the 2nd and beyond,
        %   normalize by half the oscillation period
        [averageDistance, halfPeriodsToMultiple] = ...
            compute_average_distance_to_a_multiple(lagsPeaksBins, ...
                                                    oscPeriod2Bins);

        % Compute the oscillatory index 
        oscIndex3 = 1 - averageDistance;
    end

    % Compute the oscillatory index version 4
    %   Note: This is 
    if numel(indPeaks) <= 1
        oscIndex4 = 0;
    else
        % Get the amplitude of the primary peak
        primaryPeakAmp = ampPeaks(1);

        % Find the peak with maximum amplitude other than the primary peak
        %   call this the "secondary peak"
        [~, iSecPeak] = max(ampPeaks(2:end));

        % Get the amplitude of the "secondary peak"
        secPeakAmp = ampPeaks(iSecPeak + 1);

        % Get the minimum trough amplitude between the primary peak
        %   and the secondary peak
        troughAmp = min(ampTroughs(1:iSecPeak));

        % The oscillatory index is the ratio between 
        %   the difference of maximum peak to trough in between
        %   and the different of primary peak to trough in between
        oscIndex4 = (secPeakAmp - troughAmp) / (primaryPeakAmp - troughAmp);
    end
end

%% For plotting later
% Modify the figure base
figPathBase = [fileBase, '_trace', num2str(iVec)];
figTitleBase = [figTitleBase, '\_trace', num2str(iVec)];

% Convert to seconds
[stimStartSec, detectStartSec, firstSpikeSec, ...
    histLeftSec, timeOscEndSec, oscDurationSec, ...
    maxInterBurstIntervalSec, spikeTimesSec, edgesSec, ...
    timeBurstStartsSec, timeBurstEndsSec] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            stimStartMs, detectStartMs, firstSpikeMs, ...
            histLeftMs, timeOscEndMs, oscDurationMs, ...
            maxInterBurstIntervalMs, spikeTimesMs, edgesMs, ...
            timeBurstInOscStartsMs, timeBurstInOscEndsMs);

%% Store in outputs
parsedParams.phaseNumber = phaseNumber;
parsedParams.phaseName = phaseName;
parsedParams.signal2Noise = signal2Noise;
parsedParams.minDelayMs = minDelayMs;
parsedParams.minDelaySamples = minDelaySamples;
parsedParams.binWidthMs = binWidthMs;
parsedParams.binWidthSec = binWidthSec;
parsedParams.filterWidthMs = filterWidthMs;
parsedParams.minRelProm = minRelProm;
parsedParams.minSpikeRateInBurstHz = minSpikeRateInBurstHz;
parsedParams.minBurstLengthMs = minBurstLengthMs;
parsedParams.maxInterBurstIntervalMs = maxInterBurstIntervalMs;
parsedParams.maxInterBurstIntervalSec = maxInterBurstIntervalSec;
parsedParams.siMs = siMs;
parsedParams.idxStimStart = idxStimStart;
parsedParams.stimStartMs = stimStartMs;
parsedParams.stimStartSec = stimStartSec;
parsedParams.baseWindow = baseWindow;
parsedParams.baseSlopeNoise = baseSlopeNoise;
parsedParams.slopeThreshold = slopeThreshold;
parsedParams.idxDetectStart = idxDetectStart;
parsedParams.detectStartMs = detectStartMs;
parsedParams.detectStartSec = detectStartSec;
parsedParams.nSpikesTotal = nSpikesTotal;
parsedParams.idxFirstSpike = idxFirstSpike;
parsedParams.firstSpikeMs = firstSpikeMs;
parsedParams.firstSpikeSec = firstSpikeSec;
parsedParams.vMin = vMin;
parsedParams.vMax = vMax;
parsedParams.vRange = vRange;
parsedParams.slopeMin = slopeMin;
parsedParams.slopeMax = slopeMax;
parsedParams.slopeRange = slopeRange;
parsedParams.nBins = nBins;
parsedParams.halfNBins = halfNBins;
parsedParams.histLeftMs = histLeftMs;
parsedParams.histLeftSec = histLeftSec;
parsedParams.nSpikesPerBurst = nSpikesPerBurst;
parsedParams.nSpikesPerBurstIn10s = nSpikesPerBurstIn10s;
parsedParams.nSpikesPerBurstInOsc = nSpikesPerBurstInOsc;
parsedParams.nSpikesIn10s = nSpikesIn10s;
parsedParams.nSpikesInOsc = nSpikesInOsc;
parsedParams.nBurstsTotal = nBurstsTotal;
parsedParams.nBurstsIn10s = nBurstsIn10s;
parsedParams.nBurstsInOsc = nBurstsInOsc;
parsedParams.iBinLastOfLastBurst = iBinLastOfLastBurst;
parsedParams.iBinLastOfLastBurstIn10s = iBinLastOfLastBurstIn10s;
parsedParams.iBinLastOfLastBurstInOsc = iBinLastOfLastBurstInOsc;
parsedParams.timeOscEndMs = timeOscEndMs;
parsedParams.timeOscEndSec = timeOscEndSec;
parsedParams.oscDurationMs = oscDurationMs;
parsedParams.oscDurationSec = oscDurationSec;
parsedParams.oscIndex1 = oscIndex1;
parsedParams.oscIndex2 = oscIndex2;
parsedParams.oscIndex3 = oscIndex3;
parsedParams.oscIndex4 = oscIndex4;
parsedParams.oscPeriod1Ms = oscPeriod1Ms;
parsedParams.oscPeriod2Ms = oscPeriod2Ms;
parsedParams.minOscPeriod2Bins = minOscPeriod2Bins;
parsedParams.maxOscPeriod2Bins = maxOscPeriod2Bins;
parsedParams.figPathBase = figPathBase;
parsedParams.figTitleBase = figTitleBase;

parsedData.tVec = tVec;
parsedData.vVec = vVec;
parsedData.slopes = slopes;
parsedData.idxSpikes = idxSpikes;
parsedData.spikeTimesMs = spikeTimesMs;
parsedData.spikeTimesSec = spikeTimesSec;
parsedData.spikeCounts = spikeCounts;
parsedData.edgesMs = edgesMs;
parsedData.edgesSec = edgesSec;
parsedData.spikeCountsEachBurst = spikeCountsEachBurst;
parsedData.spikeCountsEachBurstIn10s = spikeCountsEachBurstIn10s;
parsedData.spikeCountsEachBurstInOsc = spikeCountsEachBurstInOsc;
parsedData.iBinBurstStarts = iBinBurstStarts;
parsedData.iBinBurstEnds = iBinBurstEnds;
parsedData.iBinBurstIn10sStarts = iBinBurstIn10sStarts;
parsedData.iBinBurstIn10sEnds = iBinBurstIn10sEnds;
parsedData.iBinBurstInOscStarts = iBinBurstInOscStarts;
parsedData.iBinBurstInOscEnds = iBinBurstInOscEnds;
parsedData.timeBurstStartsMs = timeBurstStartsMs;
parsedData.timeBurstEndsMs = timeBurstEndsMs;
parsedData.timeBurstIn10sStartsMs = timeBurstIn10sStartsMs;
parsedData.timeBurstIn10sEndsMs = timeBurstIn10sEndsMs;
parsedData.timeBurstInOscStartsMs = timeBurstInOscStartsMs;
parsedData.timeBurstInOscEndsMs = timeBurstInOscEndsMs;
parsedData.timeBurstStartsSec = timeBurstStartsSec;
parsedData.timeBurstEndsSec = timeBurstEndsSec;
parsedData.autoCorr = autoCorr;
parsedData.acf = acf;
parsedData.acfFiltered = acfFiltered;
parsedData.acfFilteredOfInterest = acfFilteredOfInterest;
parsedData.indPeaks = indPeaks;
parsedData.indTroughs = indTroughs;
parsedData.ampPeaks = ampPeaks;
parsedData.ampTroughs = ampTroughs;
parsedData.halfPeriodsToMultiple = halfPeriodsToMultiple;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [averageDistance, halfPeriodsToMultiple] = ...
                compute_average_distance_to_a_multiple(values, period)

[~, halfPeriodsToMultiple] = ...
    arrayfun(@(x) find_nearest_multiple(period, x, ...
                                    'RelativeToHalfBase', true), values);

averageDistance = mean(halfPeriodsToMultiple);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ampTroughs, indTroughs] = find_troughs_from_peaks(vec, indPeaks)
%% Finds troughs in between given peak indices

nPeaks = numel(indPeaks);

if nPeaks < 2
    % No troughs
    ampTroughs = [];
    indTroughs = [];
else
    % Left peak indices
    indLeftPeak = indPeaks(1:(end-1));

    % Right peak indices
    indRightPeak = indPeaks(2:end);

    % Use the minimums in each interval
    [ampTroughs, indTroughsRel] = ...
        arrayfun(@(x, y) min(vec(x:y)), indLeftPeak, indRightPeak);

    % Compute the original indices
    indTroughs = arrayfun(@(x, y) x + y - 1, indTroughsRel, indLeftPeak);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iBinLastOfLastBurst, spikeCountsEachBurst, ...
                nSpikesPerBurst] = ...
                compute_burst_statistics(spikeCounts, ...
                                        iBinBurstStarts, iBinBurstEnds)

% Determine the last bin of the last burst
if ~isempty(iBinBurstEnds)
    iBinLastOfLastBurst = iBinBurstEnds(end);
else
    iBinLastOfLastBurst = NaN;
end

% Compute the number of spikes in each burst window
if ~isempty(iBinBurstStarts) && ~isempty(iBinBurstEnds)
    spikeCountsEachBurst = arrayfun(@(x, y) sum(spikeCounts(x:y)), ...
                            iBinBurstStarts, iBinBurstEnds);
else
    spikeCountsEachBurst = [];
end

% Compute average number of spikes per burst
if ~isempty(spikeCountsEachBurst)
    nSpikesPerBurst = mean(spikeCountsEachBurst);
else
    nSpikesPerBurst = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fig, ax, lines, markers, raster] = ...
                plot_spike_detection(tVec, vVec, slopes, idxSpikes, ...
                                    baseSlopeNoise, slopeThreshold, ...
                                    vMin, vMax, vRange, slopeMin, slopeMax, ...
                                    figHandle, figTitle)
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
if ~isempty(figHandle)
    fig = figure(figHandle);
else
    fig = figure('Visible', 'off');
end
clf; 

% Plot the slope trace
ax(1) = subplot(3, 1, 1);
cla; hold on
lines(1) = plot(tVec(1:(end-1)), slopes, 'k');
lines(4) = plot_horizontal_line(baseSlopeNoise, 'Color', 'b', 'LineStyle', '--');
lines(5) = plot_horizontal_line(slopeThreshold, 'Color', 'g', 'LineStyle', '--');
if ~isempty(idxSpikes)
    markers(1) = plot(tVec(idxSpikes - 1), slopes(idxSpikes - 1), 'rx', 'LineWidth', 2);
else
    markers(1) = gobjects(1);
end
ylim(yLimits1);
ylabel('Slope (V/s)');
title('Detection of peaks in the slope vector');

% Plot the original trace
ax(2) = subplot(3, 1, 2);
cla; hold on
lines(2) = plot(tVec, vVec, 'k');
if ~isempty(idxSpikes)
    markers(2) = plot(tVec(idxSpikes), vVec(idxSpikes), 'rx', 'LineWidth', 2);
else
    markers(2) = gobjects(1);
end
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

function [histBars, histFig] = ...
                plot_spike_histogram(spikeCounts, edgesSec, durationWindows, ...
                                oscDurationSec, nSpikesInOsc, ...
                                xLimitsHist, yLimitsHist, figTitleBase)


% Plot figure
histFig = figure('Visible', 'off');
hold on;
[histBars, histFig] = ...
    plot_histogram([], 'Counts', spikeCounts, 'Edges', edgesSec, ...
                    'XLimits', xLimitsHist, 'YLimits', yLimitsHist, ...
                    'XLabel', 'Time (seconds)', ...
                    'YLabel', 'Spike Count per 10 ms', ...
                    'FigTitle', ['Spike histogram for ', figTitleBase], ...
                    'FigHandle', histFig);
text(0.5, 0.95, sprintf('Oscillation Duration = %.2g seconds', ...
    oscDurationSec), 'Units', 'normalized');
text(0.5, 0.9, sprintf('Number of spikes in oscillation = %d', ...
    nSpikesInOsc), 'Units', 'normalized');
plot_horizontal_line(0, 'XLimits', durationWindows, ...
                    'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [autoCorrFig, acfFig, acfLine1, acfLine2, acfFilteredLine] = ...
                plot_autocorrelogram(autoCorr, acf, acfFiltered, indPeaks, ...
                                    indTroughs, ampPeaks, ampTroughs, ...
                                    binWidthSec, nBins, halfNBins, ...
                                    oscIndex1, oscIndex2, oscIndex3, oscIndex4, ...
                                    oscPeriod1Ms, oscPeriod2Ms, ...
                                    oscDurationSec, nSpikesInOsc, ...
                                    xLimitsAutoCorr, yLimitsAutoCorr, ...
                                    xLimitsAcfFiltered, yLimitsAcfFiltered, ...
                                    yOscDur, figTitleBase)
                                
% Create time values 
if nBins > 1
    tAcfTemp = create_time_vectors(nBins - 1, 'SamplingIntervalSec', binWidthSec, ...
                                'TimeUnits', 's');
    tAcf = [0; tAcfTemp(1:halfNBins)];
    tAutoCorr = [-flipud(tAcfTemp); 0; tAcfTemp];
    timePeaksSec = (indPeaks - indPeaks(1)) * binWidthSec;
    timeTroughsSec = (indTroughs - indPeaks(1)) * binWidthSec;
else
    tAcf = NaN(size(acf));
    tAutoCorr = NaN(size(autoCorr));
    timePeaksSec = NaN(size(ampPeaks));
    timeTroughsSec = NaN(size(ampTroughs));
end

xLimitsOscDur = [0, oscDurationSec];

% Plot the autocorrelogram
autoCorrFig = figure('Visible', 'off');
acfLine1 = plot(tAutoCorr, autoCorr);
xlim(xLimitsAutoCorr);
ylim(yLimitsAutoCorr);
xlabel('Lag (s)');
ylabel('Spike rate squared (Hz^2)');
title(['Autocorrelation for ', figTitleBase]);

% Plot the autocorrelation function
acfFig = figure('Visible', 'off');
hold on;
acfLine2 = plot(tAcf, acf, 'k');
acfFilteredLine = plot(tAcf, acfFiltered, 'g', 'LineWidth', 1);
plot(timePeaksSec, ampPeaks, 'ro', 'LineWidth', 2);
plot(timeTroughsSec, ampTroughs, 'bx', 'LineWidth', 2);
plot_horizontal_line(yOscDur, 'XLimits', xLimitsOscDur, ...
                    'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
text(0.5, 0.98, sprintf('Oscillatory Index 4 = %.2g', oscIndex4), ...
    'Units', 'normalized');
text(0.5, 0.94, sprintf('Oscillatory Index 3 = %.2g', oscIndex3), ...
    'Units', 'normalized');
text(0.5, 0.90, sprintf('Oscillatory Index 2 = %.2g', oscIndex2), ...
    'Units', 'normalized');
text(0.5, 0.86, sprintf('Oscillatory Index 1 = %.2g', oscIndex1), ...
    'Units', 'normalized');
text(0.5, 0.82, sprintf('Oscillation Period 2 = %.3g ms', oscPeriod2Ms), ...
    'Units', 'normalized');
text(0.5, 0.78, sprintf('Oscillation Period 1 = %.3g ms', oscPeriod1Ms), ...
    'Units', 'normalized');
text(0.5, 0.74, sprintf('Total spike count = %g', nSpikesInOsc), ...
    'Units', 'normalized');
text(0.5, 0.70, sprintf('Oscillation Duration = %.2g seconds', ...
    oscDurationSec), 'Units', 'normalized');
xlim(xLimitsAcfFiltered);
ylim(yLimitsAcfFiltered);
xlabel('Lag (s)');
ylabel('Spike rate squared (Hz^2)');
title(['Autocorrelation function for ', figTitleBase]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_all_zooms(fig, outFolder, figPathBase, zoomWin1, zoomWin2, zoomWin3)
%% Save the figure as .fig and 4 zooms as .png
% TODO: Make this more general

% Get the figure
figure(fig);

% Create subdirectories for different zooms
check_subdir(outFolder, {'full', 'zoom1', 'zoom2', 'zoom3'});

% Save the full figure
saveas(fig, fullfile(outFolder, 'full', [figPathBase, '_full']), 'png');

% Zoom #1
xlim(zoomWin1);
saveas(fig, fullfile(outFolder, 'zoom1', [figPathBase, '_zoom1']), 'png');

% Zoom #2
xlim(zoomWin2);
saveas(fig, fullfile(outFolder, 'zoom2', [figPathBase, '_zoom2']), 'png');

% Zoom #3
xlim(zoomWin3);
saveas(fig, fullfile(outFolder, 'zoom3', [figPathBase, '_zoom3']), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vector = prepare_for_plot_horizontal_line(starts, ends)
%% Put in the form [start1, end1, start2, end2, ..., startn, endn]

if isempty(starts) || isempty(ends)
    form1 = [0; 0];
else
    % Put in the form [start1, start2, ..., startn;
    %                   end1,   end2,  ...,  endn]
    form1 = transpose([starts, ends]);
end

% Reshape as a column vector
vector = reshape(form1, [], 1);

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

oscWindow = transpose([stimStartMs, timeOscEndMs]);
[hLines, eventTimes, yEnds, yTicksTable] = ...
    plot_raster(spikeTimesMs, 'DurationWindow', oscWindow, ...
                'LineWidth', 0.5);
vertLine = plot_vertical_line(mean(stimStartMs), 'Color', 'g', ...
                                'LineStyle', '--');
save_all_zooms(figs(nVectors + 1), figPathBaseThis, ...
                mean(stimStartMs), mean(detectStartMs), mean(firstSpikeMs));

function save_all_zooms(fig, figPathBase, stimStartMs, detectStartMs, firstSpikeMs)
%% Save the figure as .fig and 4 zooms as .png
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

% Compute the minimum spikes per bin in the last burst
minSpikesPerBinLastBurst = ceil(minSpikeRateLastBurstHz * binWidthSec);
% Find the bins with number of spikes greater than minSpikesPerBinLastBurst
binsManyManySpikes = find(spikeCounts > minSpikesPerBinLastBurst);
% Find the time of oscillation end in ms
if isempty(binsManyManySpikes)
    timeOscEndMs = stimStartMs;
else
    iBin = numel(binsManyManySpikes) + 1;
    lastBurstBin = [];
    while isempty(lastBurstBin) && iBin > 1
        % Decrement the bin number
        iBin = iBin - 1;
        % Find the last bin left with 
        %   number of spikes greater than minSpikesPerBinLastBurst
        idxBinLast = binsManyManySpikes(iBin);
        % Compute the maximum number of bins between last two bursts
        maxIbiBins = floor(maxInterBurstIntervalMs / binWidthMs);
        % Compute the last bin index that is within maxIbiBins of 
        %   the last bin with many many spikes
        idxBinMax = min(idxBinLast + maxIbiBins, nBins);
        % Determine whether each bin is within maxIbiBins of 
        %   the last bin with many many spikes
        withinIBI = false(nBins, 1);
        withinIBI(idxBinLast:idxBinMax) = true;
        % Determine whether each bin has number of spikes at least a threshold
        isManySpikes = (spikeCounts >= minSpikesPerBinLastBurst);
        % Find the last consecutive bin with number of spikes greater than threshold
        %   within maxIbiBins of the last bin with many many spikes
        %   Note: First bin must be true
        lastBurstBin = find([false; isManySpikes(1:end-1)] & isManySpikes & ...
                            withinIBI, 1, 'last');
    end
    % If still not found, the last burst bin is the one with many many spikes
    if isempty(lastBurstBin)
        lastBurstBin = idxBinLast;
    end
    % Compute the time of oscillation end in ms
    timeOscEndMs = histLeftMs + lastBurstBin * binWidthMs;
end

% Compute the minimum spikes per bin in the last burst
minSpikesPerBinInBurst = ceil(minSpikeRateInBurstHz * binWidthSec);
% Determine whether each bin passes the number of spikes criterion
isInBurst = spikeCounts >= minSpikesPerBinInBurst;
% Determine whether each bin and its previous minBinsInBurst
%   consecutive bins all pass the number of spikes criterion
isLastBinInBurst = isInBurst;
previousInBurst = isInBurst;
for i = 1:minBinsInBurst
    % Whether the previous ith bin passes the number of spikes criterion
    previousInBurst = [false; previousInBurst(1:(end-1))];

    % Whether the previous i bins all pass the number of spikes criterion
    isLastBinInBurst = isLastBinInBurst & previousInBurst;
end
% Find the last bins of each burst
iBinBurstEnds = find(isLastBinInBurst);

% Record the amplitude of the primary peak
ampPeak1 = acfFiltered(1);
% Find the index and amplitude of the secondary peak
% TODO: Use all peaks
[peakAmp, peakInd] = ...
    findpeaks(acfFiltered, 'MinPeakProminence', minRelProm * ampPeak1);
idxPeak2 = peakInd(1);
ampPeak2 = peakAmp(1);
% Find the amplitude of the first trough
[troughNegAmp, troughInd] = findpeaks(-acfFiltered);
idxTrough1 = troughInd(1);
ampTrough1 = -troughNegAmp(1);
% Compute the average amplitude of first two peaks
ampPeak12 = mean([ampPeak1, ampPeak2]);
% Compute the oscillatory index
oscIndex3 = (ampPeak12 - ampTrough1) / ampPeak12;
% Compute the oscillation period
oscPeriod2Ms = idxPeak2 * binWidthMs;
parsedParams.ampPeak1 = ampPeak1;
parsedParams.idxPeak2 = idxPeak2;
parsedParams.ampPeak2 = ampPeak2;
parsedParams.idxTrough1 = idxTrough1;
parsedParams.ampTrough1 = ampTrough1;
parsedParams.ampPeak12 = ampPeak12;
timePeak1Sec = 0;
timeTrough1Sec = idxTrough1 * binWidthSec;
timePeak2Sec = idxPeak2 * binWidthSec;
plot(timePeak1Sec, ampPeak1, 'ro', 'LineWidth', 2);
plot(timeTrough1Sec, ampTrough1, 'bx', 'LineWidth', 2);
plot(timePeak2Sec, ampPeak2, 'ro', 'LineWidth', 2);
xlim([0, 7]);

% Take just the positive side
acf = autoCorr(nBins:end);

% Create figure path base
figPathBase = fullfile(outFolder, [fileBase, '_spike_detection']);
figPathBase = [fileBase, '_spike_detection'];

minDelaySamples = 2000;
minDelayMs = 200;

check_subdir(outFolder, {rasterDir, autoCorrDir, acfDir, ...
                spikeHistDir, spikeDetectionDir});

% Compute x and y limits
acfOfInterest = acf(1:floor(7/binWidthSec));
maxAcf = max(acfOfInterest);
yLimits = compute_axis_limits({acfOfInterest, 0}, 'y', 'Coverage', 90);
yOscDur = -(maxAcf * 0.025);
xLimitsOscDur = [0, oscDurationMs - autoCorrDelayMs] / MS_PER_S;
xLimitsAcfFiltered = [0, max(timePeaksSec(end), xLimitsOscDur(2)) + 1];
% xLimitsAcfFiltered = [0, 7];

% Convert to seconds
spikeTimesSec = cellfun(@(x) x/MS_PER_S, spikeTimesMs, ...
                        'UniformOutput', false);

% Record the delay for the autocorrelogram
autoCorrDelayMs = histLeftMs - stimStartMs;
autoCorrDelaySec = autoCorrDelayMs / MS_PER_S;
parsedParams.autoCorrDelayMs = autoCorrDelayMs;
parsedParams.autoCorrDelaySec = autoCorrDelaySec;
autoCorrDelayMs = NaN;
xLimitsOscDur = [0, oscDurationSec - autoCorrDelaySec];

xLimitsAcfFiltered = [0, max(timePeaksSec(end), xLimitsOscDur(2)) + 1];

acfOfInterest = acf(1:floor(7/binWidthSec));
maxAcf = max(acfOfInterest);
yLimits = compute_axis_limits({acfOfInterest, 0}, 'y', 'Coverage', 90);
yOscDur = -(maxAcf * 0.025);

acfFilteredRight = nanmean(bestRightForAll) + 1.96 * nanstderr(bestRightForAll);
bestUpperLimit = nanmean(largestAcfValues) + 1.96 * nanstderr(largestAcfValues);

xLimitsOscDur = [histLeftSec, timeOscEndSec];

% Find the index and amplitude of the peaks
if numel(acfFiltered) > 3
    [peakAmp, peakInd] = ...
        findpeaks(acfFiltered, 'MinPeakProminence', minRelProm * ampPeak1);

    % Record all peak indices and amplitudes
    indPeaks = [1; peakInd];
    ampPeaks = [ampPeak1; peakAmp];
else
    indPeaks = 1;
    ampPeaks = ampPeak1;
end

% Compute the number of peaks
nPeaks = numel(indPeaks);

% Find the indices and amplitudes of the troughs in between each pair of peak
[ampTroughs, indTroughs] = find_troughs_from_peaks(acfFiltered, indPeaks);

% Compute the average amplitudes between adjacent peaks
ampAdjPeaks = mean([ampPeaks(1:(end-1)), ampPeaks(2:end)], 2);

% Compute the oscillatory index
oscIndex3 = mean((ampAdjPeaks - ampTroughs) ./ ampAdjPeaks);

if nPeaks > 1
    oscPeriod2Ms = (indPeaks(2) - indPeaks(1)) * binWidthMs;
else
    oscPeriod2Ms = 0;
end

% Compute the lags between adjacent peaks in ms
lagsBetweenPeaksMs = diff(indPeaks) * binWidthMs;

% Compute the oscillatory index 
%   Note: This is one over the coefficient of variation 
%           of the lag differences between adjacent peaks
if numel(lagsBetweenPeaksMs) < 2
    oscIndex3 = NaN;
else
    oscIndex3 = 1 ./ compute_stats(lagsBetweenPeaksMs, 'cov');
end

parsedData.lagsBetweenPeaksMs = lagsBetweenPeaksMs;
    lagsBetweenPeaksMs = [];

if nPeaks > 1
    [~, iPeak] = max(ampPeaks(2:end));
    oscPeriod2Bins = indPeaks(iPeak + 1) - indPeaks(1);
else
    oscPeriod2Bins = 0;
end

[~, iPeak] = max(ampPeaks(2:end));
oscPeriod1Bins = indPeaks(iPeak + 1) - indPeaks(1);

% Check if total spike count is correct
if nSpikesTotal2 ~= nSpikesTotal
    error('Code logic error!');
end

plot_traces(tVecs, vVecs, 'Verbose', false, ...
            'PlotMode', 'parallel', 'SubplotOrder', 'list', ...
            'YLimit', [-4, 4], ...
            'XLabel', xLabel, 'LinkAxesOption', 'y', ...
            'TraceLabels', 'suppress', ...
            'FigTitle', figTitle, 'FigHandle', figs(1), ...
            'Color', 'k');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%