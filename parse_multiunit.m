function varargout = parse_multiunit (vVecsOrSlice, varargin)
%% Parses multiunit recordings: detect spikes, computes spike histograms and autocorrelograms
% Usage: [parsedParams, parsedData, phaseBoundaries, fileBase, figs] = parse_multiunit (vVecsOrSlice, siMs (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       parse_multiunit(vVecs, siMs, 'PulseVectors', iVecs);
%       parse_multiunit('20190217_slice3');
%       parse_multiunit('20190217_slice3', 'SaveResults', true);
%       parse_multiunit('20190217_slice3', 'PlotRaw', true);
%       [parsedParams, parsedData, phaseBoundaries, fileBase, figs] = parse_multiunit('20190217_slice3');
%
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
%       vVecsOrSlice - original voltage vector(s) in uV
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
%                   - 'FigFolder': directory to place figures
%                   must be a string scalar or a character vector
%                   default == fullfile(outFolder, [create_time_stamp, '_', ...
%                                   outDirSuffix])
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
%                   default == 60 ms
%                   - 'MaxFirstInterBurstIntervalMs': maximum inter-burst interval (ms)
%                               between stimulation start and the first burst
%                   must be a positive scalar
%                   default == 2000 ms
%                   - 'MaxInterBurstIntervalMs': maximum inter-burst interval (ms)
%                   must be a positive scalar
%                   default == 2000 ms
%                   - 'MinSpikeRateInBurstHz': minimum spike rate in a burst (Hz)
%                   must be a positive scalar
%                   default == 100 ms
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
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in save_all_figtypes() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/alternate_elements.m
%       cd/argfun.m
%       cd/check_dir.m
%       cd/check_subdir.m
%       cd/combine_multiunit_data.m
%       cd/compute_autocorrelogram.m
%       cd/compute_axis_limits.m
%       cd/compute_default_signal2noise.m
%       cd/compute_spike_density.m
%       cd/compute_spike_histogram.m
%       cd/compute_time_window.m
%       cd/compute_stats.m
%       cd/copy_into.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/create_subplots.m
%       cd/create_time_stamp.m
%       cd/create_time_vectors.m
%       cd/detect_spikes_multiunit.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_nearest_multiple.m
%       cd/force_column_cell.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%       cd/isfigtype.m
%       cd/ispositivescalar.m
%       cd/match_time_info.m
%       cd/parse_stim.m
%       cd/plot_autocorrelogram.m
%       cd/plot_bar.m
%       cd/plot_horizontal_line.m
%       cd/plot_raster.m
%       cd/plot_raw_multiunit.m
%       cd/plot_spike_density_multiunit.m
%       cd/plot_spike_histogram.m
%       cd/plot_table.m
%       cd/save_all_figtypes.m
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
% 2019-08-29 Added 'FigFolder' as an optional argument
% 2019-11-30 Added 'FigTypes' as an optional argument
% 2019-12-05 Fixed units of Voltage label to uV

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};
plotTypeMeasures = 'bar'; %'tuning';
figDirSuffix = 'individual_measures';
yAmountToStagger = [];                  % y amount to stagger for the raw plots
% yAmountToStagger = 10;
zoomWinRelStimStartSec = [-1; 20];      % for zoom window 1
% zoomWinRelStimStartSec = [-1; 10];
zoomWinRelDetectStartSec = [-0.2; 2];   % for zoom window 2
zoomWinRelFirstSpikeSec = [0; 0.1];     % for zoom window 3
% sweepDurationSec = 60;                  % sweep duration in seconds
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
channelTypes = {'voltage', 'current'};
channelUnits = {'uV', 'arb'};

% TODO: Make optional argument
backupSheets = true;

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
figFolderDefault = '';                  % set later
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
minBurstLengthMsDefault = 60;       % bursts must be at least 60 ms by default
maxFirstInterBurstIntervalMsDefault = 2000;
                                    % first burst is not more than 2 seconds 
                                    %   after stimulation start by default
maxInterBurstIntervalMsDefault = 2000;  % subsequent bursts are no more than 
                                        %   2 seconds apart by default
minSpikeRateInBurstHzDefault = 100; % bursts must have a spike rate of 
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

% Figure stuff
figTypesDefault = 'png';                % plot just as .png files by default 

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
addParameter(iP, 'FigFolder', figFolderDefault, ...
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
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
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
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

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
figFolder = iP.Results.FigFolder;
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
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

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

% Decide on the figure directory
if isempty(figFolder)
    % Create time stamp
    dateStamp = create_time_stamp('FormatOut', 'yyyymmdd');

    % Create figure folder name
    figFolder = fullfile(outFolder, [dateStamp, '_', figDirSuffix]);
end

% Check if output directories exist
check_dir({outFolder, figFolder});

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
if ~isfile(resultsPath)
    toExtractData = false;
end

% Extract data if needed
if toExtractData
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
            combine_multiunit_data('SliceBase', fileBase, ...
                                    'ChannelTypes', channelTypes, ...
                                    'ChannelUnits', channelUnits, ...
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
if toExtractData && isempty(vVecs)
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

% Match time information
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
    stimParams = parse_stim(pulseVectors, 'SamplingIntervalMs', siMs, ...
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

    % Backup the spreadsheet in the figure folder
    if backupSheets
        copy_into(paramsPath, figFolder);
    end

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
    figFolderSpikeDetection = fullfile(figFolder, spikeDetectionDir);

    % Plot and save figures
    plot_all_spike_detections(parsedData, parsedParams, ...
                                figFolderSpikeDetection, figTypes);
end

%% Plot spike histograms
if plotSpikeHistogramFlag
    fprintf('Plotting all spike histograms for %s ...\n', fileBase);
    % Create output directory
    figFolderHist = fullfile(figFolder, spikeHistDir);

    % Plot and save figures
    plot_all_spike_histograms(parsedData, parsedParams, ...
                                figFolderHist, figTypes);
end

%% Plot autocorrelograms
if plotAutoCorrFlag
    fprintf('Plotting all autocorrelograms for %s ...\n', fileBase);
    % Create output directories
    figFolderAutoCorr = fullfile(figFolder, autoCorrDir);
    figFolderAcf = fullfile(figFolder, acfDir);

    % Plot and save figures
    plot_all_autocorrelograms(parsedData, parsedParams, ...
                                figFolderAutoCorr, figFolderAcf, figTypes);
end

%% Plot raw traces
if plotRawFlag
    fprintf('Plotting raw traces for %s ...\n', fileBase);
    % Hard-coded parameters
    figExpansion = 2;

    % Create a figure base
    figBaseRaw = fullfile(figFolder, rawDir, [fileBase, '_raw']);

    % Plot figure
    figs(1) = figure(1); clf
    plot_raw_multiunit(parsedData, parsedParams, ...
                        'PlotFiltered', true, 'PlotMode', 'staggered', ...
                        'YAmountToStagger', yAmountToStagger, ...
                        'FigExpansion', figExpansion);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(1), figBaseRaw, zoomWinsMulti, 'FigTypes', figTypes);
end

%% Plot raster plot
if plotRasterFlag
    fprintf('Plotting raster plot for %s ...\n', fileBase);

    % Create a figure base
    figBaseRaster = fullfile(figFolder, rasterDir, [fileBase, '_raster']);

    % Plot figure
    figs(2) = figure(2); clf
    plot_raster_multiunit(parsedData, parsedParams, ...
                            phaseBoundaries, fileBase);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(2), figBaseRaster, zoomWinsMulti, 'FigTypes', figTypes);
end

%% Plot spike density plot
if plotSpikeDensityFlag
    fprintf('Plotting spike density plot for %s ...\n', fileBase);

    % Create a figure base
    figBaseSpikeDensity = fullfile(figFolder, spikeDensityDir, ...
                                    [fileBase, '_spike_density']);

    % Plot figure
    figs(3) = figure(3); clf
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Save the figure zoomed to several x limits
    save_all_zooms(figs(3), figBaseSpikeDensity, zoomWinsMulti, ...
                    'FigTypes', figTypes);
end

%% Plot combined plots
if plotCombinedFlag
    fprintf('Plotting a combined plot for %s ...\n', fileBase);    

    % Hard-Coded Parameters
    iTraceToSample = 21;
    MS_PER_S = 1000;
    vertBarWidth2Range = 1/10;

    % Create output directory and subdirectories for each measure
    figFolderCombined = fullfile(figFolder, combinedDir);
    check_dir(figFolderCombined);

    % Create a figure base
    figBaseCombined = fullfile(figFolderCombined, [fileBase, '_combined']);

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
    vertBarWidth = vRange * vertBarWidth2Range;
    yMid = vMax + vertBarWidth;
    yLimits = compute_axis_limits([vMin, yMid], 'y', 'Coverage', 100);

    % Create a new figure with 9 x 9 subplots
    close(figure(4));
    [figCombined, axCombined] = create_subplots(3, 3, 'FigNumber', 4);

    % Plot raw data with zoomWin1
    % TODO: Plot with oscillation duration bars
    subplot(axCombined(1));
    plot_raw_multiunit(parsedData, parsedParams, ...
                        'PlotFiltered', true, 'PlotMode', 'staggered', ...
                        'YAmountToStagger', yAmountToStagger, ...
                        'XLimits', zoomWin1);

    % Plot spike density with zoomWin1
    subplot(axCombined(2));
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 'XLimits', zoomWin1, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Plot oscillation duration
    subplot(axCombined(3));
    plot_bar(oscDurationSec, ...
            'PhaseBoundaries', phaseBoundaries, ...
            'ForceVectorAsRow', false, ...
            'ReverseOrder', true, ...
            'BarDirection', 'horizontal', ...
            'PLabel', 'Time (min)', ...
            'ReadoutLabel', 'Oscillation Duration (s)', ...
            'ReadoutLimits', [0, range(zoomWin1)]);

    % Plot raw data with zoomWin2
    subplot(axCombined(4));
    plot_raw_multiunit(parsedData, parsedParams, ...
                        'PlotFiltered', true, 'PlotMode', 'staggered', ...
                        'YAmountToStagger', yAmountToStagger, ...
                        'XLimits', zoomWin2);

    % Plot spike density with zoomWin2
    subplot(axCombined(5));
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                                 'XLimits', zoomWin2, ...
                                 'PlotStimStart', true, ...
                                 'BoundaryType', 'horizontalLines', ...
                                 'MaxNYTicks', 20);

    % Plot oscillation period
    subplot(axCombined(6));
    plot_bar(oscPeriod2Ms, ...
            'PhaseBoundaries', phaseBoundaries, ...
            'ForceVectorAsRow', false, ...
            'ReverseOrder', true, ...
            'BarDirection', 'horizontal', ...
            'PLabel', 'Time (min)', ...
            'ReadoutLabel', 'Oscillation Period (ms)');

    % Plot spike detection
    subplot(axCombined(7)); hold on;
    plot(tVec, vVec, 'k');
    plot(tVec, vVecFilt, 'b');
    plot_raster(tVec(idxSpikes), 'YMid', yMid, 'VertBarWidth', vertBarWidth, ...
                'LineWidth', 0.5, 'ColorMap', 'Red', 'PlotOnly', true);
    xlim(zoomWin2 * MS_PER_S);
    ylim(yLimits);
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    title(['Original voltage vector for ', figTitleBase]);

    % Plot spike histogram
    subplot(axCombined(8));
    plot_spike_histogram(sampleDataStruct, sampleParamsStruct, ...
                            'XLimits', zoomWin2);

    % Plot autocorrelogram
    subplot(axCombined(9));
    plot_autocorrelogram(sampleDataStruct, sampleParamsStruct, ...
                        'PlotType', 'acfFiltered', 'PlotFiltered', true, ...
                        'PlotPeaks', true, 'PlotTroughs', true, ...
                        'PlotDuration', true, 'PlotText', true);

    % Save the figure
    save_all_figtypes(figCombined, figBaseCombined, figTypes);

    % Output figure
    figs(4) = figCombined;
end

%% Plot contour plot
if plotContourFlag
    fprintf('Plotting contour plot for %s ...\n', fileBase);

    % Create a figure base
    check_subdir(figFolder, contourDir);
    figBaseContour = fullfile(figFolder, contourDir, ...
                                [fileBase, '_contour']);

    % Plot figure
    figs(5) = set_figure_properties('ClearFigure', true, ...
                                    'Width', 1100, 'Height', 300);
    xLimitsSeconds = [2.2, 20];
    plot_spike_density_multiunit(parsedData, parsedParams, ...
                        'XLimits', xLimitsSeconds, ...
                        'PlotStimStart', false, ...
                        'BoundaryType', 'verticalBars', ...
                        'MaxNYTicks', 10);

    % Save the figure
    save_all_figtypes(figs(5), figBaseContour, figTypes);
end

%% Plot time series of measures
if plotMeasuresFlag
    fprintf('Plotting time series of measures for %s ...\n', fileBase);    

    % Create output directory and subdirectories for each measure
    figFolderMeasures = fullfile(figFolder, measuresDir);

    % Check if output directory exists
    check_dir(figFolderMeasures);
    check_subdir(figFolderMeasures, measuresToPlot);

    % Create full figure paths
    figPathsMeasures = fullfile(figFolderMeasures, measuresToPlot, ...
                                strcat(fileBase, '_', measuresToPlot));

    % Create custom figure titles
    titleBase = replace(fileBase, '_', '\_');
    figTitlesMeasures = strcat(measuresToPlot, [' for ', titleBase]);

    % Create new figure
    figure(6);

    % Plot table and save figures
    handles = plot_table(parsedParams, 'PlotMode', 'separate', ...
                        'PlotType', plotTypeMeasures, ...
                        'PhaseBoundaries', phaseBoundaries, ...
                        'VarsToPlot', measuresToPlot, ...
                        'PLabel', 'Time (min)', ...
                        'FigNames', figPathsMeasures, ...
                        'FigTitles', figTitlesMeasures, ...
                        'NLastOfPhase', nSweepsLastOfPhase, ...
                        'NToAverage', nSweepsToAverage, ...
                        'SelectionMethod', selectionMethod, ...
                        'MaxRange2Mean', maxRange2Mean);
    figs(5 + (1:nMeasures)) = handles.fig;
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

function handles = plot_spike_detection(tVec, vVec, vVecFilt, ...
                                    slopes, idxSpikes, ...
                                    baseSlopeNoise, slopeThreshold, ...
                                    vMin, vMax, vRange, slopeMin, slopeMax, ...
                                    figHandle, figTitle)
%% Plots the spike detection
% TODO: Organize arguments to use parsedParams and parsedData instead
% TODO: Pull out as a function

% Hard-coded constants
vertBarWidth2Range = 1/10;

% Compute the midpoint and bar width for the raster
% vertBarWidth = vRange * vertBarWidth2Range;
vertBarWidth = 1;
% yMid = vMax + vertBarWidth;
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
ylabel('Voltage (uV)');
title('Corresponding positions in the voltage vector');

% Plot the original trace with spike raster
ax(3) = subplot(3, 1, 3);
cla; hold on
lines(4) = plot(tVec, vVec, 'k');
lines(5) = plot(tVec, vVecFilt, 'b');
raster = plot_raster(tVec(idxSpikes), 'YMid', yMid, 'VertBarWidth', vertBarWidth, ...
                    'LineWidth', 0.5, 'ColorMap', 'Red', ...
                    'YLimits', 'suppress', 'XLabel', 'suppress', ...
                    'YTickLocs', 'suppress', 'YTickLabels', 'suppress', ...
                    'FigTitle', 'suppress');
ylim(yLimits3);
xlabel('Time (ms)');
ylabel('Voltage (uV)');
title('Original voltage vector with spikes');

% Create an overarching title
suptitle(figTitle);

% Link the x axes
linkaxes(ax, 'x');

%% Save handles in a structure
handles.fig = fig;
handles.ax = ax;
handles.lines = lines;
handles.markers = markers;
handles.raster = raster;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_spike_detections(parsedData, parsedParams, ...
                                        figFolder, figTypes)

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
check_dir(figFolder);

% Count the number of sweeps
nVectors = height(parsedParams);

parfor iVec = 1:nVectors
    % Create a figure
    figHandle = set_figure_properties;
    
    % Plot spike detection
    plot_spike_detection(tVec{iVec}, vVec{iVec}, vVecFilt{iVec}, ...
                            slopes{iVec}, idxSpikes{iVec}, ...
                            baseSlopeNoise(iVec), slopeThreshold(iVec), ...
                            vMin(iVec), vMax(iVec), vRange(iVec), ...
                            slopeMin(iVec), slopeMax(iVec), ...
                            figHandle, figTitleBase{iVec});

    % Get the current figure path base
    figBaseThis = fullfile(figFolder, figPathBase{iVec});

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
    save_all_zooms(figHandle, figBaseThis, zoomWins, 'FigTypes', figTypes);

    % Close all figures
    close all force hidden
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_spike_histograms(parsedData, parsedParams, ...
                                        figFolder, figTypes)

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
check_dir(figFolder);

%% Do the job
% Convert to structure arrays
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@table2struct, parsedParams, parsedData);

% Plot histograms
% for iVec = 1:nVectors
parfor iVec = 1:nVectors
    thisParams = parsedParamsStruct(iVec);
    thisData = parsedDataStruct(iVec);

    % Plot the histogram
    histFig = figure('Visible', 'off');
    plot_spike_histogram(thisData, thisParams, 'XLimits', xLimitsHist, ...
                        'YLimits', yLimitsHist);

    % Save the figure
    save_all_figtypes(histFig, fullfile(figFolder, ...
                    [figPathBase{iVec}, '_spike_histogram']), figTypes);
    close all force hidden
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_all_autocorrelograms(parsedData, parsedParams, ...
                                    figFolderAutoCorr, figFolderAcf, figTypes)

% Retrieve data
acf = parsedData.acf;
indPeaks = parsedData.indPeaks;

binWidthSec = parsedParams.binWidthSec;
oscDurationSec = parsedParams.oscDurationSec;
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
check_dir(figFolderAutoCorr);
check_dir(figFolderAcf);

%% Do the job
% Convert to structure arrays
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@table2struct, parsedParams, parsedData);

% Plot autocorrelograms
parfor iVec = 1:nVectors
    thisParams = parsedParamsStruct(iVec);
    thisData = parsedDataStruct(iVec);

    autoCorrFig = figure('Visible', 'off');
    plot_autocorrelogram(thisData, thisParams, ...
            'PlotType', 'autocorrelogram', ...
            'XLimits', xLimitsAutoCorr, ...
            'YLimits', yLimitsAutoCorr);
    save_all_figtypes(autoCorrFig, fullfile(figFolderAutoCorr, ...
            [figPathBase{iVec}, '_autocorrelogram']), figTypes);

    acfFig = figure('Visible', 'off');
    plot_autocorrelogram(thisData, thisParams, ...
            'PlotType', 'acfFiltered', 'PlotFiltered', true, ...
            'PlotPeaks', true, 'PlotTroughs', true, ...
            'PlotDuration', true, 'PlotText', true, ...
            'XLimits', xLimitsAcfFiltered, ...
            'YLimits', yLimitsAcfFiltered, ...
            'BarYValue', yOscDur);
    save_all_figtypes(acfFig, fullfile(figFolderAcf, ...
            [figPathBase{iVec}, '_autocorrelation_function']), figTypes);

    close all force hidden
end

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
% TODO
% nSweeps = numel(spikeTimesSec);

% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

% Create figure and plot
hold on
hLines = plot_raster(spikeTimesSec, 'HorzBarWindows', burstWindows, ...
                    'LineWidth', 0.5, 'ColorMap', 'Black', ...
                    'XLabel', 'Time (s)', 'YLabel', 'Trace #', ...
                    'FigTitle', ['Spike times for ', titleBase]);
% [hLines, eventTimes, yEnds, yTicksTable] = ...
%     plot_raster(spikeTimesSec, 'HorzBarWindows', oscWindow, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute x limits for durations
%    durationWindows = force_column_cell(transpose([histLeft, timeOscEndSec]));
burstWindows = ...
    cellfun(@(x, y) alternate_elements(x, y, 'ReturnNaNIfEmpty', true), ...
            timeBurstStartsSec, timeBurstEndsSec, ...
            'UniformOutput', false);

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

colorsRaster = repmat({'Black'}, nSweeps, 1);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
