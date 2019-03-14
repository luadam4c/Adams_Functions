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
%                   - 'SetBoundaries': vector of set boundaries
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
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_time_info.m
%       cd/movingaveragefilter.m
%       cd/plot_raster.m
%       cd/plot_horizontal_line.m
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
% 

%% Hard-coded parameters
rasterDir = 'rasters';
autoCorrDir = 'autocorrelograms';
acfDir = 'smoothed_autocorrelograms';
spikeHistDir = 'spike_histograms';
spikeDetectionDir = 'spike_detection';
measuresDir = 'measures';
measuresToPlot = {'oscIndex', 'oscDurationSec', ...
                    'oscPeriodMs', 'spikeCountTotal'};

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
addParameter(iP, 'SetBoundaries', setBoundariesDefault, ...
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
setBoundaries = iP.Results.SetBoundaries;
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
figs = gobjects(1 + nMeasures, 1);

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
                                fileBase, titleBase, setBoundaries);
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

    binWidthSec = parsedParams.binWidthSec;
    histLeftSec = parsedParams.histLeftSec;
    timeOscEndSec = parsedParams.timeOscEndSec;
    maxInterBurstIntervalSec = parsedParams.maxInterBurstIntervalSec;
    oscDurationSec = parsedParams.oscDurationSec;
    spikeCountTotal = parsedParams.spikeCountTotal;
    figTitleBase = parsedParams.figTitleBase;
    figPathBase = parsedParams.figPathBase;

    % Find appropriate x limits
    histLeft = min(histLeftSec);
    % histRight = nanmean(timeOscEndSec) + 1.96 * stderr(timeOscEndSec) + ...
    %                 1.5 * max(maxInterBurstIntervalSec);
    histRight = 10;
    xLimitsHist = [histLeft, histRight];

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

    parfor iVec = 1:nVectors
        [histBars, histFig] = ...
            plot_spike_histogram(spikeCounts{iVec}, edgesSec{iVec}, ...
                                histLeftSec(iVec), timeOscEndSec(iVec), ...
                                oscDurationSec(iVec), spikeCountTotal(iVec), ...
                                xLimitsHist, yLimitsHist, figTitleBase{iVec});

        saveas(histFig, fullfile(outFolderHist, ...
                        [figPathBase{iVec}, '_spike_histogram']), 'png');
        close all force hidden
    end
end

%% Plot autocorrelograms
%if plotFlag
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
    oscIndexOld = parsedParams.oscIndexOld;
    oscIndex = parsedParams.oscIndex;
    oscPeriodMs = parsedParams.oscPeriodMs;
    oscDurationSec = parsedParams.oscDurationSec;
    spikeCountTotal = parsedParams.spikeCountTotal;
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
    outFolderAutoCorrFunc = fullfile(outFolder, acfDir);
    check_dir(outFolderAutoCorrFunc);

    parfor iVec = 1:nVectors
        [autoCorrFig, acfFig] = ...
            plot_autocorrelogram(autoCorr{iVec}, acf{iVec}, acfFiltered{iVec}, ...
                indPeaks{iVec}, indTroughs{iVec}, ...
                ampPeaks{iVec}, ampTroughs{iVec}, ...
                binWidthSec(iVec), nBins(iVec), halfNBins(iVec), ...
                oscIndexOld(iVec), oscIndex(iVec), ...
                oscPeriodMs(iVec), oscDurationSec(iVec), ...
                spikeCountTotal(iVec), ...
                xLimitsAutoCorr, yLimitsAutoCorr, ...
                xLimitsAcfFiltered, yLimitsAcfFiltered, ...
                yOscDur, figTitleBase{iVec});

        saveas(autoCorrFig, fullfile(outFolderAutoCorr, ...
                [figPathBase{iVec}, '_autocorrelogram']), 'png');
        saveas(acfFig, fullfile(outFolderAutoCorrFunc, ...
                [figPathBase{iVec}, '_smoothed_autocorrelogram']), 'png');

        close all force hidden
    end
%end

%% Plot raster plot
if plotFlag
    fprintf('Plotting raster plot for %s ...\n', fileBase);

    % Modify the figure base
    figBaseRaster = [fileBase, '_raster'];

    % Extract the spike times
    spikeTimesSec = parsedData.spikeTimesSec;
    stimStartSec = parsedParams.stimStartSec;
    detectStartSec = parsedParams.detectStartSec;
    firstSpikeSec = parsedParams.firstSpikeSec;
    timeOscEndSec = parsedParams.timeOscEndSec;

    % Oscillation window
    oscWindow = transpose([stimStartSec, timeOscEndSec]);

    % Burst windows
    % TODO burstWindows = 

    % Create output directory
    outFolderRaster = fullfile(outFolder, rasterDir);
    check_dir(outFolderRaster);

    % Create figure and plot
    figs(1) = figure('Visible', 'off');
    clf
    [hLines, eventTimes, yEnds, yTicksTable] = ...
        plot_raster(spikeTimesSec, 'DurationWindow', oscWindow, ...
                    'LineWidth', 0.5);
    vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                    'LineStyle', '--');
    if ~isempty(setBoundaries)
        yBoundaries = nVectors - setBoundaries + 1;
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
    save_all_zooms(figs(1), outFolderRaster, ...
                    figBaseRaster, zoomWin1, zoomWin2, zoomWin3);
end

%% Plot time series of measures
%if plotFlag
    fprintf('Plotting time series of measures for %s ...\n', fileBase);    

    % Create output directory and subdirectories for each measure
    outFolderMeasures = fullfile(outFolder, measuresDir);
    check_dir(outFolderMeasures);
    check_subdir(outFolderMeasures, measuresToPlot);

    % Create full figure paths
    figPathsMeasures = fullfile(outFolderMeasures, measuresToPlot, ...
                                strcat(fileBase, '_', measuresToPlot));

    % Create custom figure titles
    figTitlesMeasures = strcat(measuresToPlot, ' for ', titleBase);

    % Plot table
    figs(2:nMeasures + 1) = ...
        plot_table(parsedParams, 'VariableNames', measuresToPlot, ...
                    'XLabel', 'Time (min)', 'FigNames', figPathsMeasures, ...
                    'FigTitles', figTitlesMeasures, ...
                    'XBoundaries', setBoundaries, ...
                    'RemoveOutliers', true);
%end

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;
varargout{3} = figs;

fprintf('%s analyzed! ...\n\n', fileBase);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
                parse_multiunit_helper(iVec, vVec, tVec, siMs, ...
                                idxStimStart, stimStartMs, baseWindow, ...
                                fileBase, figTitleBase, setBoundaries)

% Parse a single multiunit recording

% Hard-coded constants
MS_PER_S = 1000;

% Hard-coded parameters
signal2Noise = 4; %2.5
minDelayMs = 25;
binWidthMs = 10;
filterWidthMs = 100;
minRelProm = 0.02;
minSpikeRateInBurstHz = 100;
minBurstLengthMs = 20;
maxInterBurstIntervalMs = 2000;

%% Preparation
% Compute the minimum delay in samples
minDelaySamples = ceil(minDelayMs ./ siMs);

% Compute the bin width in seconds
binWidthSec = binWidthMs ./ MS_PER_S;

% Compute the number of set boundaries
nBoundaries = numel(setBoundaries);

% Determine which set number this sweep belongs to
if nBoundaries > 0
    % For the first n - 1 sets, use find
    setNumber = find(setBoundaries > iVec, 1, 'first');

    % For the last set, use numel(setBoundaries) + 1
    if isempty(setNumber)
        setNumber = numel(setBoundaries) + 1;
    end
else
    setNumber = NaN;
end

% Create set names
if nBoundaries == 2
    if setNumber == 1
        setName = 'baseline';
    elseif setNumber == 2
        setName = 'washon';
    elseif setNumber == 3
        setName = 'washoff';
    end
else
    setName = '';
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
spikeCountTotal = numel(idxSpikes);

% Index of first spike
if spikeCountTotal == 0
    idxFirstSpike = NaN;
else
    idxFirstSpike = idxSpikes(1);
end

% Store spike times
if spikeCountTotal == 0
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

%% Compute the spike histogram
if spikeCountTotal == 0
    spikeCounts = [];
    edgesMs = [];
    nBins = 0;
    halfNBins = 0;
    histLeftMs = NaN;
    iBinLastOfLastBurst = NaN;
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
    for i = 1:minBinsInBurst
        % Compute spike counts from the previous ith bin
        spikeCountsPrev = [false; spikeCountsPrev(1:(end-1))];

        % Add the spike counts in the previous i bins
        spikeCountsWin = spikeCountsWin + spikeCountsPrev;
    end

    % Determine whether each sliding window passes the number of spikes criterion
    isInBurst = spikeCountsWin >= minSpikesPerWindowInBurst;

    % Find the last bins of each burst
    %   Note: this sliding window in burst but the next sliding window not in burst
    iBinLastInBurst = find(isInBurst & [~isInBurst(2:end); true]);

    % Find the last bin of the last burst, using maxIbiBins
    if isempty(iBinLastInBurst)
        iBinLastOfLastBurst = NaN;
    else
        % Compute the inter-burst intervals in bins
        ibiBins = diff(iBinLastInBurst);

        % Find the first inter-burst interval greater than maxIbiBins
        iBurstLast = find(ibiBins > maxIbiBins, 1, 'first');

        % Determine the last bin of the last burst
        if isempty(iBurstLast)
            % All bursts are close enough together
            iBinLastOfLastBurst = iBinLastInBurst(end);
        else
            % Actually (iBurstLast - 1) + 1
            iBinLastOfLastBurst = iBinLastInBurst(iBurstLast);
        end
    end

    % Find the time of oscillation end in ms
    if isnan(iBinLastOfLastBurst)
        timeOscEndMs = stimStartMs;
    else
        % Compute the time of oscillation end in ms
        timeOscEndMs = histLeftMs + iBinLastOfLastBurst * binWidthMs;
    end

    % Compute the oscillation duration in ms
    oscDurationMs = timeOscEndMs - stimStartMs;
end

%% Compute the average spikes per burst
% Compute the burst windows in bins
% TODO

% Compute the burst windows in ms
% TODO

% Compute the number of spikes in each burst window
% TODO

% Compute average number of spikes per burst
% TODO

%% Compute the autocorrelogram
if spikeCountTotal == 0
    oscIndexOld = NaN;
    oscIndex = NaN;
    oscPeriodMs = NaN;
    autoCorr = [];
    acf = [];
    acfFiltered = [];
    indPeaks = [];
    indTroughs = [];
    ampPeaks = [];
    ampTroughs = [];
    lagsBetweenPeaksMs = [];
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

    % Restrict the autocorrelation function to oscillation duration
    acfFilteredOfInterest = acfFiltered(1:(1+oscDurationBins));

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

    % Find the indices and amplitudes of the troughs in between each pair of peak
    [ampTroughs, indTroughs] = ...
        find_troughs_from_peaks(acfFilteredOfInterest, indPeaks);

    % Compute the average amplitudes between adjacent peaks
    ampAdjPeaks = mean([ampPeaks(1:(end-1)), ampPeaks(2:end)], 2);

    % Compute the old oscillatory index
    %   Note: This is the average of all oscillatory indices as defined
    %           by Sohal's paper between adjacent peaks
    oscIndexOld = mean((ampAdjPeaks - ampTroughs) ./ ampAdjPeaks);

    % Compute the lags between adjacent peaks in ms
    lagsBetweenPeaksMs = diff(indPeaks) * binWidthMs;

    % Compute the oscillatory index 
    %   Note: This is one over the coefficient of variation 
    %           of the lag differences between adjacent peaks
    if numel(lagsBetweenPeaksMs) < 2
        oscIndex = NaN;
    else
        oscIndex = 1 ./ compute_stats(lagsBetweenPeaksMs, 'cov');
    end

    % Compute the oscillation period
    %   Note: This is the period between the first peak 
    %           and the next largest peak
    if nPeaks > 1
        [~, iPeak] = max(ampPeaks(2:end));
        oscPeriodMs = (indPeaks(iPeak + 1) - indPeaks(1)) * binWidthMs;
    else
        oscPeriodMs = 0;
    end
end

%% For plotting later
% Modify the figure base
figPathBase = [fileBase, '_trace', num2str(iVec)];
figTitleBase = [figTitleBase, '\_trace', num2str(iVec)];

% Convert to seconds
[stimStartSec, detectStartSec, firstSpikeSec, ...
    histLeftSec, timeOscEndSec, oscDurationSec, ...
    maxInterBurstIntervalSec, spikeTimesSec, edgesSec] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            stimStartMs, detectStartMs, firstSpikeMs, ...
            histLeftMs, timeOscEndMs, oscDurationMs, ...
            maxInterBurstIntervalMs, spikeTimesMs, edgesMs);

%% Store in outputs
parsedParams.setNumber = setNumber;
parsedParams.setName = setName;
parsedParams.signal2Noise = signal2Noise;
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
parsedParams.spikeCountTotal = spikeCountTotal;
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
parsedParams.iBinLastOfLastBurst = iBinLastOfLastBurst;
parsedParams.timeOscEndMs = timeOscEndMs;
parsedParams.timeOscEndSec = timeOscEndSec;
parsedParams.oscDurationMs = oscDurationMs;
parsedParams.oscDurationSec = oscDurationSec;
parsedParams.oscIndexOld = oscIndexOld;
parsedParams.oscIndex = oscIndex;
parsedParams.oscPeriodMs = oscPeriodMs;
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
parsedData.autoCorr = autoCorr;
parsedData.acf = acf;
parsedData.acfFiltered = acfFiltered;
parsedData.acfFilteredOfInterest = acfFilteredOfInterest;
parsedData.indPeaks = indPeaks;
parsedData.indTroughs = indTroughs;
parsedData.ampPeaks = ampPeaks;
parsedData.ampTroughs = ampTroughs;
parsedData.lagsBetweenPeaksMs = lagsBetweenPeaksMs;

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
                plot_spike_histogram(spikeCounts, edgesSec, ...
                                histLeftSec, timeOscEndSec, ...
                                oscDurationSec, spikeCountTotal, ...
                                xLimitsHist, yLimitsHist, figTitleBase)

% Compute things
xLimitsOscDur = [xLimitsHist(1), timeOscEndSec];

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
text(0.5, 0.9, sprintf('Total number of spikes = %d', ...
    spikeCountTotal), 'Units', 'normalized');
plot_horizontal_line(0, 'XLimits', xLimitsOscDur, ...
                    'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [autoCorrFig, acfFig, acfLine1, acfLine2, acfFilteredLine] = ...
                plot_autocorrelogram(autoCorr, acf, acfFiltered, indPeaks, ...
                                    indTroughs, ampPeaks, ampTroughs, ...
                                    binWidthSec, nBins, halfNBins, ...
                                    oscIndexOld, oscIndex, oscPeriodMs, ...
                                    oscDurationSec, spikeCountTotal, ...
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
text(0.5, 0.95, sprintf('Oscillatory Index = %g', oscIndex), ...
    'Units', 'normalized');
text(0.5, 0.91, sprintf('Old oscillatory Index = %g', oscIndexOld), ...
    'Units', 'normalized');
text(0.5, 0.87, sprintf('Total spike count = %g', spikeCountTotal), ...
    'Units', 'normalized');
text(0.5, 0.83, sprintf('Oscillation Period = %g ms', oscPeriodMs), ...
    'Units', 'normalized');
text(0.5, 0.79, sprintf('Oscillation Duration = %.2g seconds', ...
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
iBinLastInBurst = find(isLastBinInBurst);

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
oscIndex = (ampPeak12 - ampTrough1) / ampPeak12;
% Compute the oscillation period
oscPeriodMs = idxPeak2 * binWidthMs;
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
oscIndex = mean((ampAdjPeaks - ampTroughs) ./ ampAdjPeaks);

if nPeaks > 1
    oscPeriodMs = (indPeaks(2) - indPeaks(1)) * binWidthMs;
else
    oscPeriodMs = 0;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%