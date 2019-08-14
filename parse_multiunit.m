function varargout = parse_multiunit (vVecsOrSlice, varargin)
%% Parses multiunit recordings: detect spikes, computes spike histograms and autocorrelograms
% Usage: [parsedParams, parsedData, phaseBoundaries, fileBase, figs] = parse_multiunit (vVecsOrSlice, siMs (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       parse_multiunit(vVecs, siMs, 'PulseVectors', iVecs);
%       parse_multiunit('20190217_slice3');
%       parse_multiunit('20190217_slice3', 'SaveResults', true);
%       parse_multiunit('20190217_slice3', 'PlotRaw', true);
%       [parsedParams, parsedData, phaseBoundaries, fileBase, figs] = parse_multiunit('20190217_slice3');
% Outputs:
%       parsedParams - parsed parameters, a table with columns:
%                       phaseNumber
%                       phaseName
%                       signal2Noise
%                       minDelayMs
%                       minDelaySamples
%                       binWidthMs
%                       binWidthSec
%                       filterWidthMs
%                       minRelProm
%                       minSpikeRateInBurstHz
%                       minBurstLengthMs
%                       maxInterBurstIntervalMs
%                       maxInterBurstIntervalSec
%                       siMs
%                       idxStimStart
%                       stimStartMs
%                       stimStartSec
%                       baseWindow
%                       baseSlopeNoise
%                       slopeThreshold
%                       idxDetectStart
%                       detectStartMs
%                       detectStartSec
%                       nSpikesTotal
%                       idxFirstSpike
%                       firstSpikeMs
%                       firstSpikeSec
%                       vMin
%                       vMax
%                       vRange
%                       slopeMin
%                       slopeMax
%                       slopeRange
%                       nBins
%                       halfNBins
%                       histLeftMs
%                       histLeftSec
%                       nSpikesPerBurst
%                       nSpikesPerBurstIn10s
%                       nSpikesPerBurstInOsc
%                       nSpikesIn10s
%                       nSpikesInOsc
%                       nBurstsTotal
%                       nBurstsIn10s
%                       nBurstsInOsc
%                       iBinLastOfLastBurst
%                       iBinLastOfLastBurstIn10s
%                       iBinLastOfLastBurstInOsc
%                       timeOscEndMs
%                       timeOscEndSec
%                       oscDurationMs
%                       oscDurationSec
%                       oscIndex1
%                       oscIndex2
%                       oscIndex3
%                       oscIndex4
%                       oscPeriod1Ms
%                       oscPeriod2Ms
%                       minOscPeriod2Bins
%                       maxOscPeriod2Bins
%                       figPathBase
%                       figTitleBase
%                   specified as a table
%       parsedData  - parsed parameters, a table with columns:
%                       tVec
%                       vVec
%                       vVecFilt
%                       slopes
%                       idxSpikes
%                       spikeTimesMs
%                       spikeTimesSec
%                       spikeCounts
%                       edgesMs
%                       edgesSec
%                       spikeCountsEachBurst
%                       spikeCountsEachBurstIn10s
%                       spikeCountsEachBurstInOsc
%                       iBinBurstStarts
%                       iBinBurstEnds
%                       iBinBurstIn10sStarts
%                       iBinBurstIn10sEnds
%                       iBinBurstInOscStarts
%                       iBinBurstInOscEnds
%                       timeBurstStartsMs
%                       timeBurstEndsMs
%                       timeBurstIn10sStartsMs
%                       timeBurstIn10sEndsMs
%                       timeBurstInOscStartsMs
%                       timeBurstInOscEndsMs
%                       timeBurstStartsSec
%                       timeBurstEndsSec
%                       autoCorr
%                       acf
%                       acfFiltered
%                       acfFilteredOfInterest
%                       indPeaks
%                       indTroughs
%                       ampPeaks
%                       ampTroughs
%                       halfPeriodsToMultiple
%                   specified as a table
%       phaseBoundaries     
%                   - phase boundaries
%                   specified as a numeric vector
%       fileBase    - file base
%                   specified as a character vector
%       figs        - figure handles
%                   specified as a Figure object handle array
% Arguments:
%       vVecsOrSlice - original voltage vector(s) in mV
%                       or the slice name
%                   must be a numeric array or a cell array of numeric arrays
%                       or a string scalar or a character vector
%       siMs        - (opt) sampling interval in ms
%                   must be a positive vector
%                   default == 0.1 ms
%       varargin    - 'PlotAllFlag': whether to plot everything
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotCombinedFlag': whether to plot raw data, 
%                           spike density and oscillation duration together
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeDetectionFlag': whether to plot spike detection
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeDensityFlag': whether to plot spike density
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotContourFlag': whether to plot a contour plot
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeHistogramFlag': whether to plot spike histograms
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotAutoCorrFlag': whether to plot autocorrelegrams
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotRawFlag': whether to plot raw traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotRasterFlag': whether to plot raster plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotMeasuresFlag': whether to plot time series 
%                                           of measures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveMatFlag': whether to save combined data
%                                           as matfiles
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveResultsFlag': whether to save parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'OutFolder': directory to place output files
%                   must be a string scalar or a character vector
%                   default == same as inFolder
%                   - 'FileBase': base of filename (without extension)
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%                   - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == detect from pulse vector
%                   - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == parsed from .abf files
%                   - 'PhaseBoundaries': vector of phase boundaries
%                   must be a numeric vector
%                   default == detected from .abf file names
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (not used)
%                   - 'BaseWindows': baseline window(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == start to stimStart
%                   - 'RelSnrThres2Max': relative signal to noise threshold
%                                           as a proportion of maximum
%                   must be empty or a numeric vector
%                   default == 0.1
%                   - 'Signal2Noise': signal-to-noise ratio
%                   must be a empty or a positive scalar
%                   default == use compute_default_signal2noise.m
%                   - 'ResolutionMs': time difference between each data point
%                                       for the spike density plot
%                   must be a positive scalar
%                   default == 5 ms
%                   - 'FiltFreq': cutoff frequency(ies) (Hz or normalized) 
%                                   for a bandpass filter
%                   must be a numeric a two-element vector
%                   default == [100, 1000]
%                   - 'MinDelayMs': minimum delay after stim start (ms)
%                   must be a positive scalar
%                   default == 25 ms
%                   - 'BinWidthMs': bin width (ms)
%                   must be a positive scalar
%                   default == 10 ms
%                   - 'MinBurstLengthMs': minimum burst length (ms)
%                   must be a positive scalar
%                   default == 20 ms
%                   - 'MaxFirstInterBurstIntervalMs': maximum inter-burst interval (ms)
%                                                   between the first two bursts
%                   must be a positive scalar
%                   default == 2000 ms
%                   - 'MaxInterBurstIntervalMs': maximum inter-burst interval (ms)
%                   must be a positive scalar
%                   default == 1000 ms
%                   - 'MinSpikeRateInBurstHz': minimum spike rate in a burst (Hz)
%                   must be a positive scalar
%                   default == 100 Hz
%                   - 'FilterWidthMs': filter width (ms) for 
%                                       moving average filter when computing
%                                       smoothed autocorrelogram
%                   must be a positive scalar
%                   default == 100 ms
%                   - 'MinRelProm': minimum relative prominence
%                   must be a positive scalar
%                   default == 0.02
%                   - 'NSweepsLastOfPhase': number of sweeps at 
%                                           the last of a phase
%                   must be a positive integer scalar
%                   default == 10
%                   - 'NSweepsToAverage': number of sweeps to average
%                   must be a positive integer scalar
%                   default == 5
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'maxRange2Mean'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 40%
%
% Requires:
%       cd/alternate_elements.m
%       cd/argfun.m
%       cd/check_dir.m
%       cd/check_subdir.m
%       cd/combine_data_from_same_slice.m
%       cd/compute_autocorrelogram.m
%       cd/compute_axis_limits.m
%       cd/compute_default_signal2noise.m
%       cd/compute_spike_density.m
%       cd/compute_spike_histogram.m
%       cd/compute_time_window.m
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_subplots.m
%       cd/create_time_vectors.m
%       cd/detect_spikes_multiunit.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_nearest_multiple.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%       cd/ispositivescalar.m
%       cd/match_time_info.m
%       cd/parse_stim.m
%       cd/plot_horizontal_line.m
%       cd/plot_raster.m
%       cd/plot_table.m
%       cd/save_all_zooms.m
%       cd/set_default_flag.m
%       cd/transform_vectors.m
%       ~/Downloaded_Functions/rgb.m
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
% 2019-05-03 Moved code to detect_spikes_multiunit.m
% 2019-05-06 Expanded plot flags
% 2019-05-16 Added spike density computation and plot
% 2019-05-16 Changed maxInterBurstIntervalMs to 1500
% 2019-05-16 Changed signal2Noise to 2.5 
% 2019-06-02 Added compute_default_signal2noise.m
% 2019-06-03 Moved code to save_all_zooms.m
% 2019-06-10 Compartmentalized plot code
% 2019-06-10 Added plotCombinedFlag
% 2019-07-24 Added saveResultsFlag
% 2019-07-25 Now returns phaseBoundaries and fileBase as 3rd and 4th arguments
% 2019-08-04 Now makes subplots maximally fit the figure
% 2019-08-04 Now makes 'Trace #' the y label for raw plots
% 2019-08-06 Made parameters optional arguments
% 2019-08-09 Now saves contour plots as epsc2 instead of eps
% 2019-08-13 Expanded combined plot
% 2019-08-13 Now uses alternate_elements.m

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};
plotTypeMeasures = 'bar'; %'tuning';
yAmountToStagger = [];                  % y amount to stagger for the raw plots
% yAmountToStagger = 10;
zoomWinRelStimStartSec = [-1; 20];      % for zoom window 1
% zoomWinRelStimStartSec = [-1; 10];
zoomWinRelDetectStartSec = [-0.2; 2];   % for zoom window 2
zoomWinRelFirstSpikeSec = [0; 0.1];     % for zoom window 3
sweepDurationSec = 60;                  % sweep duration in seconds
rawDir = 'raw';
rasterDir = 'rasters';
autoCorrDir = 'autocorrelograms';
acfDir = 'autocorrelation_functions';
spikeHistDir = 'spike_histograms';
spikeDensityDir = 'spike_density';
spikeDetectionDir = 'spike_detections';
measuresDir = 'measures';
combinedDir = 'combined';
contourDir = 'contour';
measuresToPlot = {'oscIndex1', 'oscIndex2', 'oscIndex3', 'oscIndex4', ...
                    'oscPeriod1Ms', 'oscPeriod2Ms', ...
                    'oscDurationSec', ...
                    'nSpikesTotal', 'nSpikesIn10s', 'nSpikesInOsc', ...
                    'nBurstsTotal', 'nBurstsIn10s', 'nBurstsInOsc', ...
                    'nSpikesPerBurst', 'nSpikesPerBurstIn10s', ...
                    'nSpikesPerBurstInOsc'};
paramsSuffix = '_params.csv';
resultsSuffix = '_parsed.mat';
varsNeeded = {'sliceBase', 'vVecsSl', 'siMsSl', 'iVecsSl', ...
                'phaseBoundaries', 'phaseStrs'};

%% Default values for optional arguments
siMsDefault = 0.1;                      % 0.1 ms by default
plotAllFlagDefault = false;
plotCombinedFlagDefault = [];
plotSpikeDetectionFlagDefault = [];
plotSpikeDensityFlagDefault = [];
plotContourFlagDefault = [];
plotSpikeHistogramFlagDefault = [];
plotAutoCorrFlagDefault = [];
plotRawFlagDefault = [];
plotRasterFlagDefault = [];
plotMeasuresFlagDefault = [];
saveMatFlagDefault = true;
saveResultsFlagDefault = false;
directoryDefault = pwd;
inFolderDefault = '';                   % set later
outFolderDefault = '';                  % set later
fileBaseDefault = 'unnamed';            % set later
stimStartMsDefault = [];                % set later
pulseVectorsDefault = [];               % don't use pulse vectors by default
phaseBoundariesDefault = [];   	        % no phase boundaries by default
tVecsDefault = [];                      % set later
baseWindowsDefault = [];                % set later
relSnrThres2MaxDefault = [];            % set in compute_default_signal2noise.m
signal2NoiseDefault = [];               % set later
resolutionMsDefault = 5;                % 5 ms resolution by default

% Note: Must be consistent with compute_oscillation_duration.m
filtFreqDefault = [100, 1000];
minDelayMsDefault = 25;
binWidthMsDefault = 10;                 % use a bin width of 10 ms by default
minBurstLengthMsDefault = 20;           % bursts must be at least 20 ms by default
maxFirstInterBurstIntervalMsDefault = 2000;
                                    % first two bursts are no more than 
                                    %   2 seconds apart by default
maxInterBurstIntervalMsDefault = 1000;  % bursts are no more than 
                                        %   1 second apart by default
minSpikeRateInBurstHzDefault = 100;     % bursts must have a spike rate of 
                                        %   at least 100 Hz by default

filterWidthMsDefault = 100;
minRelPromDefault = 0.02;

% Note: must be consistent with plot_measures.m
nSweepsLastOfPhaseDefault = 10;         % select from last 10 values by default
nSweepsToAverageDefault = 5;            % select 5 values by default
selectionMethodDefault = 'maxRange2Mean';   
                                        % select using maxRange2Mean by default
maxRange2MeanDefault = 40;              % range is not more than 40% of mean 
                                        %   by default

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
addRequired(iP, 'vVecsOrSlice', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x) || ...
                isstring(x) || ischar(x), ...
                ['vVecsOrSlice must be either a numeric array, ', ...
                    'a cell array of numeric arrays, ', ...
                    'a string scalar or a character vector!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'siMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotCombinedFlag', plotCombinedFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeDetectionFlag', plotSpikeDetectionFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeDensityFlag', plotSpikeDensityFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotContourFlag', plotContourFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeHistogramFlag', plotSpikeHistogramFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAutoCorrFlag', plotAutoCorrFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotRawFlag', plotRawFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotRasterFlag', plotRasterFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMeasuresFlag', plotMeasuresFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveResultsFlag', saveResultsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
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
addParameter(iP, 'PhaseBoundaries', phaseBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'BaseWindows', baseWindowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['baseWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'RelSnrThres2Max', relSnrThres2MaxDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'ResolutionMs', resolutionMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FiltFreq', filtFreqDefault, ...
    @(x) assert(any(isnan(x)) || isnumeric(x), ...
                ['FiltFreq must be either NaN ', ...
                    'or a numeric array of 2 elements!']));
addParameter(iP, 'MinDelayMs', minDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
addParameter(iP, 'BinWidthMs', binWidthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MinBurstLengthMs', minBurstLengthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MaxFirstInterBurstIntervalMs', ...
                maxFirstInterBurstIntervalMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MaxInterBurstIntervalMs', maxInterBurstIntervalMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MinSpikeRateInBurstHz', minSpikeRateInBurstHzDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'FilterWidthMs', filterWidthMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'MinRelProm', minRelPromDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'NSweepsLastOfPhase', nSweepsLastOfPhaseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'NSweepsToAverage', nSweepsToAverageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, vVecsOrSlice, varargin{:});
siMs = iP.Results.siMs;
plotAllFlag = iP.Results.PlotAllFlag;
plotCombinedFlag = iP.Results.PlotCombinedFlag;
plotSpikeDetectionFlag = iP.Results.PlotSpikeDetectionFlag;
plotSpikeDensityFlag = iP.Results.PlotSpikeDensityFlag;
plotContourFlag = iP.Results.PlotContourFlag;
plotSpikeHistogramFlag = iP.Results.PlotSpikeHistogramFlag;
plotAutoCorrFlag = iP.Results.PlotAutoCorrFlag;
plotRawFlag = iP.Results.PlotRawFlag;
plotRasterFlag = iP.Results.PlotRasterFlag;
plotMeasuresFlag = iP.Results.PlotMeasuresFlag;
saveMatFlag = iP.Results.SaveMatFlag;
saveResultsFlag = iP.Results.SaveResultsFlag;
directory = iP.Results.Directory;
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
fileBase = iP.Results.FileBase;
stimStartMs = iP.Results.StimStartMs;
pulseVectors = iP.Results.PulseVectors;
phaseBoundaries = iP.Results.PhaseBoundaries;
tVecs = iP.Results.tVecs;
baseWindows = iP.Results.BaseWindows;
relSnrThres2Max = iP.Results.RelSnrThres2Max;
signal2Noise = iP.Results.Signal2Noise;
resolutionMs = iP.Results.ResolutionMs;
filtFreq = iP.Results.FiltFreq;
minDelayMs = iP.Results.MinDelayMs;
binWidthMs = iP.Results.BinWidthMs;
minBurstLengthMs = iP.Results.MinBurstLengthMs;
maxFirstInterBurstIntervalMs = iP.Results.MaxFirstInterBurstIntervalMs;
maxInterBurstIntervalMs = iP.Results.MaxInterBurstIntervalMs;
minSpikeRateInBurstHz = iP.Results.MinSpikeRateInBurstHz;
filterWidthMs = iP.Results.FilterWidthMs;
minRelProm = iP.Results.MinRelProm;
nSweepsLastOfPhase = iP.Results.NSweepsLastOfPhase;
nSweepsToAverage = iP.Results.NSweepsToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;

%% Preparation
fprintf('Decide on the file base ...\n');

% Decide on the input directory
if isempty(inFolder)
    inFolder = directory;
end

% Decide on the output directory
if isempty(outFolder)
    outFolder = inFolder;
end

% Deal with the first argument
if ischar(vVecsOrSlice) || isstring(vVecsOrSlice)
    % The first argument is the slice name
    fileBase = vVecsOrSlice;

    % Data to be extracted
    toExtractData = true;
else
    % The first argument is the voltage vectors
    vVecs = vVecsOrSlice;

    % Data already provided
    toExtractData = false;
end

% Create the full path to the parameters file
paramsPath = fullfile(outFolder, [fileBase, paramsSuffix]);

% Create the full path to the results file
resultsPath = fullfile(outFolder, [fileBase, resultsSuffix]);

% Extract data only if results not provided
if toExtractData && ~isfile(resultsPath)
    % Construct the .mat file expected to exist
    regexpSliceMatFile = ['.*', fileBase, '.mat'];

    % Construct the .mat file expected
    [~, allMatPaths] = ...
        all_files('Directory', inFolder, 'RegExp', regexpSliceMatFile, ...
                    'SortBy', 'date', 'ForceCellOutput', true);

    % Load or combine data
    if numel(allMatPaths) > 1
        fprintf('Too many .mat files matching %s!\n', regexpSliceMatFile);
        varargout{1} = [];
        varargout{2} = [];
        varargout{3} = [];
        varargout{4} = fileBase;
        varargout{5} = gobjects(1);
        return;
    elseif numel(allMatPaths) == 1
        fprintf('Loading data for %s ...\n', fileBase);
        % Load data for each slice as a structure array
        allData = load(allMatPaths{1}, varsNeeded{:});
    else
        % Combine data from .abf files
        fprintf('Combining data for %s ...\n', fileBase);
        allData = ...
            combine_data_from_same_slice('SliceBase', fileBase, ...
                                        'SaveMatFlag', saveMatFlag, ...
                                        'VarsToSave', varsNeeded);
    end

    % Extract from the table
    %   Note: this will override whatever is provided in optional arguments
    vVecs = allData.vVecsSl;
    siMs = allData.siMsSl;
    pulseVectors = allData.iVecsSl;
    phaseBoundaries = allData.phaseBoundaries;
end

% Return if no data loaded
if isempty(vVecs)
    fprintf('No data for %s found!\n', fileBase);
    return
end

% Set default flags
fprintf('Setting default flags for %s ...\n', fileBase);
[plotSpikeDetectionFlag, plotSpikeDensityFlag, ...
plotSpikeHistogramFlag, plotAutoCorrFlag, ...
plotRawFlag, plotRasterFlag, plotMeasuresFlag, ...
plotCombinedFlag, plotContourFlag] = ...
    argfun(@(x) set_default_flag(x, plotAllFlag), ...
                plotSpikeDetectionFlag, plotSpikeDensityFlag, ...
                plotSpikeHistogramFlag, plotAutoCorrFlag, ...
                plotRawFlag, plotRasterFlag, plotMeasuresFlag, ...
                plotCombinedFlag, plotContourFlag);


fprintf('Initialize plotting parameters for %s ...\n', fileBase);
% Count the number of measures to plot
nMeasures = numel(measuresToPlot);

% Initialize figures array
figs = gobjects(nMeasures + 4, 1);

if ~isfile(resultsPath)
    fprintf('Match time information for %s ...\n', fileBase);
    % Count the number of vectors
    nVectors = count_vectors(vVecs);

    % Count the number of samples for each vector
    nSamples = count_samples(vVecs);

    % Match time vector(s) with sampling interval(s) and number(s) of samples
    [tVecs, siMs, ~] = match_time_info(tVecs, siMs, nSamples);
end

%% Do the job
if isfile(resultsPath)
    fprintf('Loading already parsed results from %s ...\n', resultsPath);
    load(resultsPath, 'parsedParams', 'parsedData', ...
                      'phaseBoundaries', 'fileBase');
else
    % Detect stimulation start time if not provided
    %   Otherwise find the corresponding index in the time vector
    fprintf('Detecting stimulation start for %s ...\n', fileBase);
    [stimParams, stimData] = ...
        parse_stim(pulseVectors, 'SamplingIntervalMs', siMs, ...
                        'StimStartMs', stimStartMs, 'tVecs', tVecs);

    idxStimStart = stimParams.idxStimStart;
    stimStartMs = stimParams.stimStartMs;

    % Compute the minimum delay in samples
    minDelaySamples = round(minDelayMs ./ siMs);

    % Find the starting index for detecting a spike
    idxDetectStart = idxStimStart + minDelaySamples;

    % Construct default baseline windows
    if isempty(baseWindows)
        fprintf('Constructing baseline window for %s ...\n', fileBase);
        baseWindows = compute_time_window(tVecs, 'TimeEnd', stimStartMs);
    end

    % Determine a slice-dependent signal-to-noise ratio if not provided
    %   Note: This assumes all sweeps have the same protocol
    if isempty(signal2Noise)
        fprintf('Determining signal-to-noise ratio for %s ...\n', fileBase);
        signal2Noise = compute_default_signal2noise(vVecs, ...
                            'siMs', siMs, 'tVecs', tVecs, ...
                            'IdxDetectStart', idxDetectStart, ...
                            'BaseWindows', baseWindows, ...
                            'FiltFreq', filtFreq, ...
                            'RelSnrThres2Max', relSnrThres2Max);
        fprintf('Signal-to-noise ratio to use: %g\n', signal2Noise);
    end

    % Force as a cell array of vectors
    [vVecs, tVecs, baseWindows] = ...
        argfun(@force_column_cell, vVecs, tVecs, baseWindows);

    % Parse all of them in a parfor loop
    fprintf('Parsing recording for %s ...\n', fileBase);
    parsedParamsCell = cell(nVectors, 1);
    parsedDataCell = cell(nVectors, 1);
    parfor iVec = 1:nVectors
%    for iVec = 1:nVectors
    %for iVec = 1:1
        [parsedParamsCell{iVec}, parsedDataCell{iVec}] = ...
            parse_multiunit_helper(iVec, vVecs{iVec}, tVecs{iVec}, siMs(iVec), ...
                                    idxStimStart(iVec), stimStartMs(iVec), ...
                                    baseWindows{iVec}, ...
                                    filtFreq, filterWidthMs, ...
                                    minDelayMs, binWidthMs, ...
                                    resolutionMs, signal2Noise, ...
                                    minBurstLengthMs, ...
                                    maxFirstInterBurstIntervalMs, ...
                                    maxInterBurstIntervalMs, ...
                                    minSpikeRateInBurstHz, minRelProm, ...
                                    fileBase, phaseBoundaries);
    end

    % Convert to a struct array
    %   Note: This removes all entries that are empty
    [parsedParamsStruct, parsedDataStruct] = ...
        argfun(@(x) [x{:}], parsedParamsCell, parsedDataCell);

    % Convert to a table
    [parsedParams, parsedData] = ...
        argfun(@(x) struct2table(x, 'AsArray', true), ...
                parsedParamsStruct, parsedDataStruct);

    % Save the parsed parameters table
    writetable(parsedParams, paramsPath);

    % Save the results
    if saveResultsFlag
        save(resultsPath, 'parsedParams', 'parsedData', ...
                          'phaseBoundaries', 'fileBase', '-v7.3');
    end
end

%% Prepare for plotting
% Determine zoom windows for multi-trace plots
if plotRawFlag || plotRasterFlag || plotSpikeDensityFlag || plotCombinedFlag
    % Retrieve params
    stimStartSec = parsedParams.stimStartSec;
    detectStartSec = parsedParams.detectStartSec;
    firstSpikeSec = parsedParams.firstSpikeSec;

    % Set zoom windows
    zoomWin1 = mean(stimStartSec) + zoomWinRelStimStartSec;
    zoomWin2 = mean(detectStartSec) + zoomWinRelDetectStartSec;
    meanFirstSpike = nanmean(firstSpikeSec);
    if ~isnan(meanFirstSpike)
        zoomWin3 = meanFirstSpike + zoomWinRelFirstSpikeSec;
    else
        zoomWin3 = zoomWinRelFirstSpikeSec;
    end

    % Combine zoom windows
    zoomWinsMulti = [zoomWin1, zoomWin2, zoomWin3];
end

%% Plot spike detection figures
if plotSpikeDetectionFlag
    fprintf('Plotting all spike detection plots for %s ...\n', fileBase);

    % Create output directory
    outFolderSpikeDetection = fullfile(outFolder, spikeDetectionDir);

    % Plot and save figures
    plot_all_spike_detections(parsedData, parsedParams, outFolderSpikeDetection);
end

%% Plot spike histograms
if plotSpikeHistogramFlag
    fprintf('Plotting all spike histograms for %s ...\n', fileBase);
    % Create output directory
    outFolderHist = fullfile(outFolder, spikeHistDir);

    % Plot and save figures
    plot_all_spike_histograms(parsedData, parsedParams, outFolderHist);
end

%% Plot autocorrelograms
if plotAutoCorrFlag
    fprintf('Plotting all autocorrelograms for %s ...\n', fileBase);
    % Create output directories
    outFolderAutoCorr = fullfile(outFolder, autoCorrDir);
    outFolderAcf = fullfile(outFolder, acfDir);

    % Plot and save figures
    plot_all_autocorrelograms(parsedData, parsedParams, ...
                                outFolderAutoCorr, outFolderAcf);
end

%% Plot raw traces
if plotRawFlag
    fprintf('Plotting raw traces for %s ...\n', fileBase);

    % Create a figure base
    figBaseRaw = fullfile(outFolder, rawDir, [fileBase, '_raw']);

    % Plot figure
    figs(1) = figure(1); clf
    plot_raw_multiunit(parsedData, parsedParams, ...
                        phaseBoundaries, fileBase, ...
                        'YAmountToStagger', yAmountToStagger);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(1), figBaseRaw, zoomWinsMulti);
end

%% Plot raster plot
if plotRasterFlag
    fprintf('Plotting raster plot for %s ...\n', fileBase);

    % Create a figure base
    figBaseRaster = fullfile(outFolder, rasterDir, [fileBase, '_raster']);

    % Plot figure
    figs(2) = figure(2); clf
    plot_raster_multiunit(parsedData, parsedParams, ...
                            phaseBoundaries, fileBase);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(2), figBaseRaster, zoomWinsMulti);
end

%% Plot spike density plot
if plotSpikeDensityFlag
    fprintf('Plotting spike density plot for %s ...\n', fileBase);

    % Create a figure base
    figBaseSpikeDensity = fullfile(outFolder, spikeDensityDir, ...
                                    [fileBase, '_spike_density']);

    % Plot figure
    figs(3) = figure(3); clf
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 phaseBoundaries, fileBase, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(3), figBaseSpikeDensity, zoomWinsMulti);
end

%% Plot combined plots
if plotCombinedFlag
    fprintf('Plotting a combined plot for %s ...\n', fileBase);    

    % Hard-Coded Parameters
    iTraceToSample = 21;
    MS_PER_S = 1000;
    barWidth2Range = 1/10;

    % Create output directory and subdirectories for each measure
    outFolderCombined = fullfile(outFolder, combinedDir);
    check_dir(outFolderCombined);

    % Create a figure base
    figBaseCombined = fullfile(outFolderCombined, [fileBase, '_combined']);

    % Extract the oscillation durations and periods
    oscDurationSec = parsedParams.oscDurationSec;
    oscPeriod2Ms = parsedParams.oscPeriod2Ms;

    % Extract sample data
    [parsedParamsStruct, parsedDataStruct] = ...
        argfun(@table2struct, parsedParams, parsedData);
    sampleDataStruct = parsedDataStruct(iTraceToSample);
    sampleParamsStruct = parsedParamsStruct(iTraceToSample);
    tVec = sampleDataStruct.tVec;
    vVec = sampleDataStruct.vVec;
    vVecFilt = sampleDataStruct.vVecFilt;
    idxSpikes = sampleDataStruct.idxSpikes;
    vRange = sampleParamsStruct.vRange;
    vMax = sampleParamsStruct.vMax;
    vMin = sampleParamsStruct.vMin;
    figTitleBase = sampleParamsStruct.figTitleBase;

    % Compute for spike detection plot
    barWidth = vRange * barWidth2Range;
    yMid = vMax + barWidth;
    yLimits = compute_axis_limits([vMin, yMid], 'y', 'Coverage', 100);

    % Compute the baseline average and indices selected for this field
    [baselineAverageOscDur, indSelectedOscDur] = ...
        compute_phase_average(oscDurationSec, ...
                    'PhaseBoundaries', phaseBoundaries, ...
                    'PhaseNumber', 1, ...
                    'NLastOfPhase', nSweepsLastOfPhase, ...
                    'NToAverage', nSweepsToAverage, ...
                    'SelectionMethod', selectionMethod, ...
                    'MaxRange2Mean', maxRange2Mean);
    [baselineAverageOscPeriod, indSelectedOscPeriod] = ...
        compute_phase_average(oscPeriod2Ms, ...
                    'PhaseBoundaries', phaseBoundaries, ...
                    'PhaseNumber', 1, ...
                    'NLastOfPhase', nSweepsLastOfPhase, ...
                    'NToAverage', nSweepsToAverage, ...
                    'SelectionMethod', selectionMethod, ...
                    'MaxRange2Mean', maxRange2Mean);

    % Create a new figure with 9 x 9 subplots
    close(figure(4));
    [figCombined, axCombined] = create_subplots(3, 3, 'FigNumber', 4);

%{
    % Plot raw data with zoomWin1
    axes(axCombined(1, 1));
    plot_raw_multiunit(parsedData, parsedParams, ...
                        phaseBoundaries, fileBase, ...
                        'YAmountToStagger', yAmountToStagger, ...
                        'XLimits', zoomWin1);

    % Plot spike density with zoomWin1
    axes(axCombined(1, 2));
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 phaseBoundaries, fileBase, ...
                                 'XLimits', zoomWin1, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Plot oscillation duration
    axes(axCombined(1, 3));
    plot_bar(oscDurationSec, ...
                'ForceVectorAsRow', false, ...
                'ReverseOrder', true, ...
                'BarDirection', 'horizontal', ...
                'PLabel', 'Time (min)', ...
                'ReadoutLabel', 'Oscillation Duration (s)', ...
                'ReadoutLimits', [0, range(zoomWin1)], ...
                'PBoundaries', phaseBoundaries, ...
                'RBoundaries', baselineAverageOscDur, ...
                'IndSelected', indSelectedOscDur);

    % Plot raw data with zoomWin2
    axes(axCombined(2, 1));
    plot_raw_multiunit(parsedData, parsedParams, ...
                        phaseBoundaries, fileBase, ...
                        'YAmountToStagger', yAmountToStagger, ...
                        'XLimits', zoomWin2);

    % Plot spike density with zoomWin2
    axes(axCombined(2, 2));
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 phaseBoundaries, fileBase, ...
                                 'XLimits', zoomWin2, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Plot oscillation period
    axes(axCombined(2, 3));
    plot_bar(oscPeriod2Ms, ...
                'ForceVectorAsRow', false, ...
                'ReverseOrder', true, ...
                'BarDirection', 'horizontal', ...
                'PLabel', 'Time (min)', ...
                'ReadoutLabel', 'Oscillation Period (ms)', ...
                'PBoundaries', phaseBoundaries, ...
                'RBoundaries', baselineAverageOscPeriod, ...
                'IndSelected', indSelectedOscPeriod);
%}

    % Plot spike detection
    axes(axCombined(3, 1)); hold on;
    plot(tVec, vVec, 'k');
    plot(tVec, vVecFilt, 'b');
    plot_raster(tVec(idxSpikes), 'YMid', yMid, 'BarWidth', barWidth, ...
                        'LineWidth', 0.5, 'Colors', {'Red'}, ...
                        'YLimits', 'suppress', 'YTickLocs', 'suppress', ...
                        'YTickLabels', 'suppress');
    xlim(zoomWin2 * MS_PER_S);
    ylim(yLimits);
    xlabel('Time (ms)');
    ylabel('Voltage (mV)');
    title(['Original voltage vector for ', figTitleBase]);

    % Plot spike histogram
    axes(axCombined(3, 2));
    plot_spike_histogram(sampleParamsStruct, sampleDataStruct, ...
                          'XLimits', zoomWin2);

    % Plot autocorrelogram
    axes(axCombined(3, 3));
    plot_autocorrelogram(sampleParamsStruct, sampleDataStruct);

    % Save the figure
    saveas(figCombined, figBaseCombined, 'png');

    % Output figure
    figs(4) = figCombined;
end

%% Plot contour plot
if plotContourFlag
    fprintf('Plotting contour plot for %s ...\n', fileBase);

    % Create a figure base
    check_subdir(outFolder, contourDir);
    figBaseContour = fullfile(outFolder, contourDir, ...
                                [fileBase, '_contour']);

    % Plot figure
    figs(5) = figure(5); clf
    figPosition = [300, 300, 1100, 300];
    xLimitsSeconds = [2.2, 10];
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                        phaseBoundaries, fileBase, ...
                        'XLimits', xLimitsSeconds, ...
                        'FigPosition', figPosition, ...
                        'PlotStimStart', false, ...
                        'BoundaryType', 'verticalBar', ...
                        'MaxNYTicks', 10);

    % Save the figure as an eps file
    save_all_figtypes(figs(5), figBaseContour, {'epsc2', 'png'});
end

%% Plot time series of measures
if plotMeasuresFlag
    fprintf('Plotting time series of measures for %s ...\n', fileBase);    

    % Create output directory and subdirectories for each measure
    outFolderMeasures = fullfile(outFolder, measuresDir);

    % Check if output directory exists
    check_dir(outFolderMeasures);
    check_subdir(outFolderMeasures, measuresToPlot);

    % Create full figure paths
    figPathsMeasures = fullfile(outFolderMeasures, measuresToPlot, ...
                                strcat(fileBase, '_', measuresToPlot));

    % Create custom figure titles
    titleBase = replace(fileBase, '_', '\_');
    figTitlesMeasures = strcat(measuresToPlot, [' for ', titleBase]);

    % Create new figure
    figure(6);

    % Plot table and save figures
    figs(5 + (1:nMeasures)) = ...
        plot_table(parsedParams, 'PlotType', plotTypeMeasures, ...
                    'VariableNames', measuresToPlot, ...
                    'PLabel', 'Time (min)', 'FigNames', figPathsMeasures, ...
                    'FigTitles', figTitlesMeasures, ...
                    'PBoundaries', phaseBoundaries, ...
                    'PlotSeparately', true, ...
                    'NLastOfPhase', nSweepsLastOfPhase, ...
                    'NToAverage', nSweepsToAverage, ...
                    'MaxRange2Mean', maxRange2Mean);
end

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;
varargout{3} = phaseBoundaries;
varargout{4} = fileBase;
varargout{5} = figs;

fprintf('%s analyzed! ...\n\n', fileBase);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
                parse_multiunit_helper(iVec, vVec, tVec, siMs, ...
                                idxStimStart, stimStartMs, baseWindow, ...
                                filtFreq, filterWidthMs, ...
                                minDelayMs, binWidthMs, ...
                                resolutionMs, signal2Noise, ...
                                minBurstLengthMs, ...
                                maxFirstInterBurstIntervalMs, ...
                                maxInterBurstIntervalMs, ...
                                minSpikeRateInBurstHz, minRelProm, ...
                                fileBase, phaseBoundaries)

% Parse a single multiunit recording

%% Hard-coded constants
MS_PER_S = 1000;

%% Preparation
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

% Compute the maximum and minimum times
maxTimeMs = tVec(end);
minTimeMs = tVec(1);

%% Detect spikes
% Detect spikes (bandpass filter before detection)
%   Note: This checks whether each inflection point 
%           crosses a sweep-dependent slope threshold, 
%           given a slice-dependent signal-to-noise ratio
[spikesParams, spikesData] = ...
    detect_spikes_multiunit(vVec, siMs, ...
                            'tVec', tVec, 'IdxStimStart', idxStimStart, ...
                            'FiltFreq', filtFreq, ...
                            'BaseWindow', baseWindow, ...
                            'MinDelayMs', minDelayMs, ...
                            'Signal2Noise', signal2Noise);

idxDetectStart = spikesParams.idxDetectStart;
detectStartMs = spikesParams.detectStartMs;
baseSlopeNoise = spikesParams.baseSlopeNoise;
slopeThreshold = spikesParams.slopeThreshold;
nSpikesTotal = spikesParams.nSpikesTotal;
idxFirstSpike = spikesParams.idxFirstSpike;
firstSpikeMs = spikesParams.firstSpikeMs;
vMin = spikesParams.vMin;
vMax = spikesParams.vMax;
vRange = spikesParams.vRange;
slopeMin = spikesParams.slopeMin;
slopeMax = spikesParams.slopeMax;
slopeRange = spikesParams.slopeRange;

vVecFilt = spikesData.vVecFilt;
slopes = spikesData.slopes;
idxSpikes = spikesData.idxSpikes;
spikeTimesMs = spikesData.spikeTimesMs;

%% Compute spike density
spikeDensityHz = ...
    compute_spike_density(spikeTimesMs, 'TimeWindow', [0, maxTimeMs], ...
                            'BinWidth', binWidthMs, ...
                            'Resolution', resolutionMs, ...
                            'TimeUnits', 'ms');

%% Compute the spike histogram, spikes per burst & oscillation duration
[spHistParams, spHistData] = ...
    compute_spike_histogram(spikeTimesMs, 'StimStartMs', stimStartMs, ...
            'BinWidthMs', binWidthMs, ...
            'MinBurstLengthMs', minBurstLengthMs, ...
            'MaxFirstInterBurstIntervalMs', maxFirstInterBurstIntervalMs, ...
            'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
            'MinSpikeRateInBurstHz', minSpikeRateInBurstHz);

nBins = spHistParams.nBins;
halfNBins = spHistParams.halfNBins;
histLeftMs = spHistParams.histLeftMs;
nBurstsTotal = spHistParams.nBurstsTotal;
nBurstsIn10s = spHistParams.nBurstsIn10s;
nBurstsInOsc = spHistParams.nBurstsInOsc;
iBinLastOfLastBurst = spHistParams.iBinLastOfLastBurst;
iBinLastOfLastBurstIn10s = spHistParams.iBinLastOfLastBurstIn10s;
iBinLastOfLastBurstInOsc = spHistParams.iBinLastOfLastBurstInOsc;
nSpikesPerBurst = spHistParams.nSpikesPerBurst;
nSpikesPerBurstIn10s = spHistParams.nSpikesPerBurstIn10s;
nSpikesPerBurstInOsc = spHistParams.nSpikesPerBurstInOsc;
nSpikesIn10s = spHistParams.nSpikesIn10s;
nSpikesInOsc = spHistParams.nSpikesInOsc;
timeOscEndMs = spHistParams.timeOscEndMs;
oscDurationMs = spHistParams.oscDurationMs;

spikeCounts = spHistData.spikeCounts;
edgesMs = spHistData.edgesMs;
iBinBurstStarts = spHistData.iBinBurstStarts;
iBinBurstEnds = spHistData.iBinBurstEnds;
iBinBurstIn10sStarts = spHistData.iBinBurstIn10sStarts;
iBinBurstIn10sEnds = spHistData.iBinBurstIn10sEnds;
iBinBurstInOscStarts = spHistData.iBinBurstInOscStarts;
iBinBurstInOscEnds = spHistData.iBinBurstInOscEnds;
spikeCountsEachBurst = spHistData.spikeCountsEachBurst;
spikeCountsEachBurstIn10s = spHistData.spikeCountsEachBurstIn10s;
spikeCountsEachBurstInOsc = spHistData.spikeCountsEachBurstInOsc;
timeBurstStartsMs = spHistData.timeBurstStartsMs;
timeBurstEndsMs = spHistData.timeBurstEndsMs;
timeBurstIn10sStartsMs = spHistData.timeBurstIn10sStartsMs;
timeBurstIn10sEndsMs = spHistData.timeBurstIn10sEndsMs;
timeBurstInOscStartsMs = spHistData.timeBurstInOscStartsMs;
timeBurstInOscEndsMs = spHistData.timeBurstInOscEndsMs;

%% Compute the autocorrelogram, oscillation period & oscillatory index
% TODO: compute_autocorrelogram.m
[autoCorrParams, autoCorrData] = ...
    compute_autocorrelogram(spikeTimesMs, ...
                            'StimStartMs', stimStartMs, ...
                            'BinWidthMs', binWidthMs, ...
                            'MinBurstLengthMs', minBurstLengthMs, ...
                            'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
                            'MinSpikeRateInBurstHz', minSpikeRateInBurstHz, ...
                            'FilterWidthMs', filterWidthMs, ...
                            'MinRelProm', minRelProm, ...
                            'SpikeHistParams', spHistParams, ...
                            'SpikeHistData', spHistData);
                            
oscIndex1 = autoCorrParams.oscIndex1;
oscIndex2 = autoCorrParams.oscIndex2;
oscIndex3 =  autoCorrParams.oscIndex3;
oscIndex4 = autoCorrParams.oscIndex4;
oscPeriod1Ms = autoCorrParams.oscPeriod1Ms;
oscPeriod2Ms = autoCorrParams.oscPeriod2Ms;
minOscPeriod2Bins = autoCorrParams.minOscPeriod2Bins;
maxOscPeriod2Bins = autoCorrParams.maxOscPeriod2Bins;

autoCorr = autoCorrData.autoCorr;
acf = autoCorrData.acf;
acfFiltered = autoCorrData.acfFiltered;
acfFilteredOfInterest = autoCorrData.acfFilteredOfInterest;
indPeaks = autoCorrData.indPeaks;
indTroughs = autoCorrData.indTroughs;
ampPeaks = autoCorrData.ampPeaks;
ampTroughs = autoCorrData.ampTroughs;
halfPeriodsToMultiple = autoCorrData.halfPeriodsToMultiple;

%% For plotting later
% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

% Modify the figure base and title base
figPathBase = [fileBase, '_trace', num2str(iVec)];
figTitleBase = [titleBase, '\_trace', num2str(iVec)];

% Convert to seconds
[siSeconds, maxTimeSec, minTimeSec, ...
    stimStartSec, detectStartSec, firstSpikeSec, ...
    histLeftSec, timeOscEndSec, oscDurationSec, ...
    maxInterBurstIntervalSec, spikeTimesSec, edgesSec, ...
    timeBurstStartsSec, timeBurstEndsSec] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            siMs, maxTimeMs, minTimeMs, ...
            stimStartMs, detectStartMs, firstSpikeMs, ...
            histLeftMs, timeOscEndMs, oscDurationMs, ...
            maxInterBurstIntervalMs, spikeTimesMs, edgesMs, ...
            timeBurstInOscStartsMs, timeBurstInOscEndsMs);

%% Store in outputs
parsedParams.fileBase = fileBase;
parsedParams.sweepNumber = iVec;
parsedParams.phaseNumber = phaseNumber;
parsedParams.phaseName = phaseName;
parsedParams.filtFreq = filtFreq;
parsedParams.signal2Noise = signal2Noise;
parsedParams.minDelayMs = minDelayMs;
parsedParams.binWidthMs = binWidthMs;
parsedParams.binWidthSec = binWidthSec;
parsedParams.filterWidthMs = filterWidthMs;
parsedParams.minRelProm = minRelProm;
parsedParams.minSpikeRateInBurstHz = minSpikeRateInBurstHz;
parsedParams.minBurstLengthMs = minBurstLengthMs;
parsedParams.maxInterBurstIntervalMs = maxInterBurstIntervalMs;
parsedParams.maxInterBurstIntervalSec = maxInterBurstIntervalSec;
parsedParams.resolutionMs = resolutionMs;
parsedParams.siMs = siMs;
parsedParams.minTimeMs = minTimeMs;
parsedParams.maxTimeMs = maxTimeMs;
parsedParams.siSeconds = siSeconds;
parsedParams.minTimeSec = minTimeSec;
parsedParams.maxTimeSec = maxTimeSec;
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
parsedData.vVecFilt = vVecFilt;
parsedData.slopes = slopes;
parsedData.idxSpikes = idxSpikes;
parsedData.spikeTimesMs = spikeTimesMs;
parsedData.spikeTimesSec = spikeTimesSec;
parsedData.spikeDensityHz = spikeDensityHz;
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

function [fig, ax, lines, markers, raster] = ...
                plot_spike_detection(tVec, vVec, vVecFilt, ...
                                    slopes, idxSpikes, ...
                                    baseSlopeNoise, slopeThreshold, ...
                                    vMin, vMax, vRange, slopeMin, slopeMax, ...
                                    figHandle, figTitle)
%% Plots the spike detection

% Hard-coded constants
barWidth2Range = 1/10;

% Compute the midpoint and bar width for the raster
% barWidth = vRange * barWidth2Range;
barWidth = 1;
% yMid = vMax + barWidth;
yMid = 9;

% Compute y axis limits
yLimits1 = compute_axis_limits([slopeMin, slopeMax], 'y', 'Coverage', 100);
% yLimits1 = [-15, 15];
yLimits2 = compute_axis_limits([vMin, vMax], 'y', 'Coverage', 100);
% yLimits2 = [-10, 10];
yLimits3 = compute_axis_limits([vMin, yMid], 'y', 'Coverage', 100);
% yLimits3 = [-10, 10];

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

% Plot the slope trace with peaks
ax(1) = subplot(3, 1, 1);
cla; hold on
lines(1) = plot(tVec(1:(end-1)), slopes, 'k');
lines(6) = plot_horizontal_line(baseSlopeNoise, 'Color', 'b', 'LineStyle', '--');
lines(7) = plot_horizontal_line(slopeThreshold, 'Color', 'g', 'LineStyle', '--');
if ~isempty(idxSpikes)
    markers(1) = plot(tVec(idxSpikes - 1), slopes(idxSpikes - 1), 'rx', 'LineWidth', 2);
else
    markers(1) = gobjects(1);
end
ylim(yLimits1);
ylabel('Slope (mV/s)');
title('Detection of peaks in the slope vector');

% Plot the original trace with maximum slope locations
ax(2) = subplot(3, 1, 2);
cla; hold on
lines(2) = plot(tVec, vVec, 'k');
lines(3) = plot(tVec, vVecFilt, 'b');
if ~isempty(idxSpikes)
    markers(2) = plot(tVec(idxSpikes), vVec(idxSpikes), 'rx', 'LineWidth', 2);
else
    markers(2) = gobjects(1);
end
ylim(yLimits2);
ylabel('Voltage (mV)');
title('Corresponding positions in the voltage vector');

% Plot the original trace with spike raster
ax(3) = subplot(3, 1, 3);
cla; hold on
lines(4) = plot(tVec, vVec, 'k');
lines(5) = plot(tVec, vVecFilt, 'b');
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
                plot_spike_histogram(spHistParams, spHistData, varargin)
%% Plots a spike histogram from the results of compute_spike_histogram

%                   - 'XLimits': x-axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == minimum and maximum edges of bins
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - Any other parameter-value pair for the plot_histogram() function

%% Hard-coded parameters
minNXTicks = 5;

%% Default values for optional arguments
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'spHistParams', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addRequired(iP, 'spHistData', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, spHistParams, spHistData, varargin{:});
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the plot_histogram() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Extract from spHistData
spikeCounts = spHistData.spikeCounts;
edgesSec = spHistData.edgesSec;
timeBurstStartsSec = spHistData.timeBurstStartsSec;
timeBurstEndsSec = spHistData.timeBurstEndsSec;

% Extract from spHistParams
oscDurationSec = spHistParams.oscDurationSec;
nSpikesInOsc = spHistParams.nSpikesInOsc;
figTitleBase = spHistParams.figTitleBase;

% Compute burst windows
burstWindows = alternate_elements(timeBurstStartsSec, timeBurstEndsSec, ...
                                    'ReturnNaNIfEmpty', true);

% Compute default x limits
if isempty(xLimits)
    % Extract parameters
    histLeftSec = spHistParams.histLeftSec;
    timeOscEndSec = spHistParams.timeOscEndSec;
    maxInterBurstIntervalSec = spHistParams.maxInterBurstIntervalSec;

    % Compute left and right ends of histogram to show
    histLeft = min(histLeftSec);
    histRight = timeOscEndSec + 1.5 * max(maxInterBurstIntervalSec);
    % histRight = 10;

    % Find appropriate x limits
    xLimits = [histLeft, histRight];
end

% Compute x tick locations
xTickLocs = linspace(xLimits(1), xLimits(2), minNXTicks);

% Compute default y limits
if isempty(yLimits)
    % Extract parameters
    binWidthSec = spHistParams.binWidthSec;

    % Find the last bin to show for all traces
    lastBinToShow = floor(range(xLimits) ./ binWidthSec) + 1;

    % Find appropriate y limits
    spikeCountsOfInterest = extract_subvectors(spikeCounts, ...
                            'IndexEnd', lastBinToShow);
    largestSpikeCount = apply_iteratively(@max, spikeCountsOfInterest);
    yLimits = [0, largestSpikeCount * 1.2];
end

%% Plot
% Plot figure
hold on;
[histBars, histFig] = ...
    plot_histogram([], 'Counts', spikeCounts, 'Edges', edgesSec, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'XTickLocs', xTickLocs, 'XLabel', 'Time (seconds)', ...
                    'YLabel', 'Spike Count per 10 ms', ...
                    'FigTitle', ['Spike histogram for ', figTitleBase], ...
                    otherArguments{:});
text(0.2, 0.95, sprintf('Oscillation Duration = %.2g seconds', ...
    oscDurationSec), 'Units', 'normalized');
text(0.2, 0.9, sprintf('Number of spikes in oscillation = %d', ...
    nSpikesInOsc), 'Units', 'normalized');
plot_horizontal_line(0, 'XLimits', burstWindows, ...
                    'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);
hold off;

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
%% Plots an autocorrelation function from the results of compute_autocorrelogram.m

%                   - 'XLimits': x-axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == minimum and maximum edges of bins
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - Any other parameter-value pair for the TODO() function

%% Hard-coded parameters

%% Default values for optional arguments
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'autoCorrParams', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addRequired(iP, 'autoCorrData', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, autoCorrParams, autoCorrData, varargin{:});
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

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

% Compute the x limits for the oscillation duration line
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

function plot_all_spike_detections(parsedData, parsedParams, outFolder)

% Retrieve data for plotting
tVec = parsedData.tVec;
vVec = parsedData.vVec;
vVecFilt = parsedData.vVecFilt;
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

% Check if output directory exists
check_dir(outFolder);

% Count the number of sweeps
nVectors = height(parsedParams);

parfor iVec = 1:nVectors
    % Plot spike detection
    [fig, ax, lines, markers, raster] = ...
        plot_spike_detection(tVec{iVec}, vVec{iVec}, vVecFilt{iVec}, ...
                            slopes{iVec}, idxSpikes{iVec}, ...
                            baseSlopeNoise(iVec), slopeThreshold(iVec), ...
                            vMin(iVec), vMax(iVec), vRange(iVec), ...
                            slopeMin(iVec), slopeMax(iVec), ...
                            [], figTitleBase{iVec});

    % Get the current figure path base
    figBaseThis = fullfile(outFolder, figPathBase{iVec});

    % Set zoom windows
    zoomWin1 = stimStartMs(iVec) + [0; 1e4];
%       zoomWin1 = stimStartMs(iVec) + [0; 2e4];
    zoomWin2 = detectStartMs(iVec) + [0; 2e3];
    if ~isnan(firstSpikeMs(iVec))
        zoomWin3 = firstSpikeMs(iVec) + [0; 60];
    else
        zoomWin3 = [0; 60];
    end            

    % Put zoom windows together
    zoomWins = [zoomWin1, zoomWin2, zoomWin3];

    % Save the figure zoomed to several x limits
    save_all_zooms(fig, figBaseThis, zoomWins);

    % Close all figures
    close all force hidden
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_spike_histograms(parsedData, parsedParams, outFolder)

%% Preparation
% Count the number of sweeps
nVectors = height(parsedParams);

% Retrieve data for plotting
spikeCounts = parsedData.spikeCounts;

binWidthSec = parsedParams.binWidthSec;
histLeftSec = parsedParams.histLeftSec;
timeOscEndSec = parsedParams.timeOscEndSec;
maxInterBurstIntervalSec = parsedParams.maxInterBurstIntervalSec;
figPathBase = parsedParams.figPathBase;

% Find appropriate x limits
histLeft = min(histLeftSec);
histRight = nanmean(timeOscEndSec) + 1.96 * stderr(timeOscEndSec) + ...
                1.5 * max(maxInterBurstIntervalSec);
% histRight = 10;
xLimitsHist = [histLeft, histRight];

% Find the last bin to show for all traces
lastBinToShow = floor((histRight - histLeft) ./ binWidthSec) + 1;

% Find appropriate y limits
spikeCountsOfInterest = extract_subvectors(spikeCounts, ...
                        'IndexEnd', lastBinToShow);
largestSpikeCount = apply_iteratively(@max, spikeCountsOfInterest);
yLimitsHist = [0, largestSpikeCount * 1.2];

% Check if output directory exists
check_dir(outFolder);

%% Do the job
% Convert to structure arrays
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@table2struct, parsedParams, parsedData);

% Plot histograms
parfor iVec = 1:nVectors
    thisParams = parsedParamsStruct(iVec);
    thisData = parsedDataStruct(iVec);

    histFig = figure('Visible', 'off');
    [histBars, histFig] = ...
        plot_spike_histogram(thisParams, thisData, ...
                                'XLimits', xLimitsHist, ...
                                'YLimits', yLimitsHist);

    saveas(histFig, fullfile(outFolder, ...
                    [figPathBase{iVec}, '_spike_histogram']), 'png');
    close all force hidden
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_autocorrelograms(parsedData, parsedParams, ...
                                    outFolderAutoCorr, outFolderAcf)

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

% Count the number of sweeps
nVectors = height(parsedParams);

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

% Check if output directory exists
check_dir(outFolderAutoCorr);
check_dir(outFolderAcf);

%% Do the job
% Convert to structure arrays
% TODO: USE THIS
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@table2struct, parsedParams, parsedData);

% Plot autocorrelograms
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_raw_multiunit (parsedData, parsedParams, ...
                                phaseBoundaries, fileBase, ...
                                varargin)

%       varargin    - 'YAmountToStagger': amount to stagger 
%                                           if 'plotmode' is 'stagger'
%                   must be a positive scalar
%                   default == uses the original y axis range
%                   - Any other parameter-value pair for plot_traces()

%% Hard-coded constants
MS_PER_S = 1000;

%% Default values for optional arguments
yAmountToStaggerDefault = [];   % set later  

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
% TODO

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YAmountToStagger', yAmountToStaggerDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                ['YAmountToStagger must be either a empty ', ...
                    'or a positive scalar!']));

% Read from the Input Parser
parse(iP, varargin{:});
yAmountToStagger = iP.Results.YAmountToStagger;

% Keep unmatched arguments for the plot_traces() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Extract parameters
stimStartSec = parsedParams.stimStartSec;
detectStartSec = parsedParams.detectStartSec;
firstSpikeSec = parsedParams.firstSpikeSec;
timeOscEndSec = parsedParams.timeOscEndSec;

% Extract parameters
tVecs = parsedData.tVec;
vVecs = parsedData.vVec;
vVecsFilt = parsedData.vVecFilt;

% Count the number of sweeps
nVectors = height(parsedParams);

% Convert time vector to seconds
tVecsSec = transform_vectors(tVecs, MS_PER_S, 'divide');

% Prepare for the plot
xLabel = 'Time (s)';
titleBase = replace(fileBase, '_', '\_');
figTitle = ['Raw and Filtered traces for ', titleBase];

% Compute the original y limits from data
bestYLimits = compute_axis_limits(vVecs, 'y', 'AutoZoom', true);

% Compute a default amount of y to stagger if not provided
if isempty(yAmountToStagger)
    yAmountToStagger = range(bestYLimits);
end

%% Plot
hold on
fig = gcf;
plot_traces(tVecsSec, vVecs, 'Verbose', false, ...
            'PlotMode', 'staggered', 'SubplotOrder', 'list', ...
            'YLimits', bestYLimits, 'YAmountToStagger', yAmountToStagger, ...
            'XLabel', xLabel, 'LinkAxesOption', 'y', ...
            'YLabel', 'Trace #', 'TraceLabels', 'suppress', ...
            'FigTitle', figTitle, 'FigHandle', fig, ...
            'Color', 'k', otherArguments{:});
plot_traces(tVecsSec, vVecsFilt, 'Verbose', false, ...
            'PlotMode', 'staggered', 'SubplotOrder', 'list', ...
            'YLimits', bestYLimits, 'YAmountToStagger', yAmountToStagger, ...
            'XLabel', xLabel, 'LinkAxesOption', 'y', ...
            'YLabel', 'Trace #', 'TraceLabels', 'suppress', ...
            'FigTitle', figTitle, 'FigHandle', fig, ...
            'Color', 'b', otherArguments{:});

% Plot stimulation start
vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                'LineStyle', '--', 'LineWidth', 0.5);

% Plot phase boundaries
if ~isempty(phaseBoundaries)
    yBoundaries = (nVectors - phaseBoundaries + 1) * yAmountToStagger;
    horzLine = plot_horizontal_line(yBoundaries, 'Color', 'g', ...
                                    'LineStyle', '--', 'LineWidth', 2);
end

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_raster_multiunit (parsedData, parsedParams, ...
                                phaseBoundaries, fileBase)
%% Plots a spike raster plot from parsed multiunit data
% TODO: Plot burst duration
% TODO: Plot oscillatory index

% Extract the spike times
spikeTimesSec = parsedData.spikeTimesSec;
timeBurstStartsSec = parsedData.timeBurstStartsSec;
timeBurstEndsSec = parsedData.timeBurstEndsSec;

stimStartSec = parsedParams.stimStartSec;
timeOscEndSec = parsedParams.timeOscEndSec;

% Count the number of sweeps
nVectors = height(parsedParams);

% Convert oscillatory index to a window
% TODO

% Oscillation window
oscWindow = transpose([stimStartSec, timeOscEndSec]);

% Burst windows
burstWindows = ...
    cellfun(@(x, y) alternate_elements(x, y, 'ReturnNaNIfEmpty', true), ...
            timeBurstStartsSec, timeBurstEndsSec, ...
            'UniformOutput', false);

% Create colors
nSweeps = numel(spikeTimesSec);
colorsRaster = repmat({'Black'}, nSweeps, 1);

% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

% Create figure and plot
hold on
[hLines, eventTimes, yEnds, yTicksTable] = ...
    plot_raster(spikeTimesSec, 'DurationWindow', burstWindows, ...
                'LineWidth', 0.5, 'Colors', colorsRaster);
% [hLines, eventTimes, yEnds, yTicksTable] = ...
%     plot_raster(spikeTimesSec, 'DurationWindow', oscWindow, ...
%                 'LineWidth', 0.5, 'Colors', colorsRaster);

% Plot stimulation start
vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                'LineStyle', '--', 'LineWidth', 0.5);
% Plot phase boundaries
if ~isempty(phaseBoundaries)
    yBoundaries = nVectors - phaseBoundaries + 1;
    horzLine = plot_horizontal_line(yBoundaries, 'Color', 'g', ...
                                    'LineStyle', '--', 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('Trace #');
title(['Spike times for ', titleBase]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_spike_density_multiunit (parsedData, parsedParams, ...
                                        phaseBoundaries, fileBase, ...
                                        varargin)
%% Plots a spike density plot from parsed multiunit data

%% Hard-coded parameters
validBoundaryTypes = {'horizontalLines', 'verticalBar'};

% Default values for optional arguments
figPositionDefault = [];
xLimitsDefault = [];
plotStimStartDefault = true;
boundaryTypeDefault = 'horizontalLines';
maxNYTicksDefault = 20;

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigPosition', figPositionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PlotStimStart', plotStimStartDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BoundaryType', boundaryTypeDefault, ...
    @(x) any(validatestring(x, validBoundaryTypes)));
addParameter(iP, 'MaxNYTicks', maxNYTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
figPosition = iP.Results.FigPosition;
xLimits = iP.Results.XLimits;
plotStimStart = iP.Results.PlotStimStart;
boundaryType = validatestring(iP.Results.BoundaryType, validBoundaryTypes);
maxNYTicks = iP.Results.MaxNYTicks;

% Retrieve data for plotting
spikeDensityHz = parsedData.spikeDensityHz;

siSeconds = parsedParams.siSeconds;
minTimeSec = parsedParams.minTimeSec;
maxTimeSec = parsedParams.maxTimeSec;
stimStartSec = parsedParams.stimStartSec;

% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

% Plot as a heatmap
hold on
% TODO: plot_heat_map(spikeDensityHz);

% Compute stimulation start
meanStimStartSec = mean(stimStartSec);

% Count traces
nSweeps = numel(spikeDensityHz);

% Get the average sampling interval in seconds
siSeconds = mean(siSeconds);

% Set x and y end points
xEnds = [min(minTimeSec); max(maxTimeSec)];
yEnds = [1; nSweeps];

% Set x and y limits
if isempty(xLimits)
    xLimits = [xEnds(1) - 0.5 * siSeconds; xEnds(2) + 0.5 * siSeconds];
end
yLimits = [yEnds(1) - 0.5; yEnds(2) + 0.5];

% Decide on y ticks and labels
yTicks = create_indices('IndexEnd', nSweeps, 'MaxNum', maxNYTicks, ...
                        'AlignMethod', 'left');
yTickLabels = create_labels_from_numbers(nSweeps - yTicks + 1);

% Force as a matrix and transpose it so that
%   each trace is a row
spikeDensityMatrix = transpose(force_matrix(spikeDensityHz));

% Set a gray-scale color map
% colormap(flipud(gray));
% colormap(jet);
cm = create_colormap('ColorMapFunc', @gray, 'ReverseOrder', true, ...
                     'HighContrast', true);
colormap(cm);

% Generate plot
imagesc(xEnds, flipud(yEnds), spikeDensityMatrix);

% Set the y ticks and labels
yticks(yTicks);
yticklabels(yTickLabels);

% Plot stimulation start
if plotStimStart
    stimStartLine = plot_vertical_line(meanStimStartSec, 'Color', rgb('Green'), ...
                                    'LineStyle', '--', 'LineWidth', 0.5, ...
                                    'YLimits', yLimits);
end

% Plot phase boundaries
if ~isempty(phaseBoundaries)
    % Compute y values for phase boundaries
    yBoundaries = nSweeps - phaseBoundaries + 1;

    % Plot phase boundaries
    switch boundaryType
        case 'horizontalLines'
            boundaryLine = plot_horizontal_line(yBoundaries, ...
                                    'Color', rgb('Green'), ...
                                    'LineStyle', '--', 'LineWidth', 2, ...
                                    'XLimits', xLimits);
        case 'verticalBar'
            boundaryLine = plot_vertical_line(meanStimStartSec - 0.1, ...
                                    'Color', rgb('Green'), ...
                                    'LineStyle', '-', 'LineWidth', 3, ...
                                    'YLimits', yBoundaries);
        otherwise
            error('boundaryType unrecognized!');
    end
end

xlim(xLimits);
ylim(yLimits);
xlabel('Time (s)');
ylabel('Trace #');
title(['Spike density (Hz) for ', titleBase]);

% Show a scale bar
colorbar;

% Set figure position if requested
if ~isempty(figPosition)
    set(gcf, 'Position', figPosition);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute x limits for durations
%    durationWindows = force_column_cell(transpose([histLeft, timeOscEndSec]));
burstWindows = ...
    cellfun(@(x, y) alternate_elements(x, y, 'ReturnNaNIfEmpty', true), ...
            timeBurstStartsSec, timeBurstEndsSec, ...
            'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%