function [autoCorrParams, autoCorrData] = compute_autocorrelogram (spikeTimesMs, varargin)
%% Computes an autocorrelogram and compute the oscillatory index and period from an array of event times
% Usage: [autoCorrParams, autoCorrData] = compute_autocorrelogram (spikeTimesMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [autoCorrParams, autoCorrData] = compute_autocorrelogram(1000*rand(100, 1));
%
% Outputs:
%       autoCorrParams  - autocorrelogram parameters, with fields:
%                           oscIndex1 (Sohal 2003)
%                           oscIndex2
%                           oscIndex3
%                           oscIndex4 (ClC2, m3ha)
%                           oscPeriod1Ms (Sohal 2003)
%                           oscPeriod2Ms (ClC2, m3ha)
%                           minOscPeriod2Bins
%                           maxOscPeriod2Bins
%                       specified as a scalar structure
%       autoCorrData    - autocorrelogram data, with fields:
%                           autoCorr
%                           acf
%                           acfFiltered
%                           acfFilteredOfInterest
%                           indPeaks
%                           indTroughs
%                           ampPeaks
%                           ampTroughs
%                           halfPeriodsToMultiple
%                       specified as a scalar structure
%
% Arguments:
%       spikeTimesMs    - spike times in milliseconds
%                       must be a numeric vector
%       varargin    - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == 0 ms
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
%                   - 'FilterWidthMs': filter width (ms) for moving average filter
%                   must be a positive scalar
%                   default == 100 ms
%                   - 'MinRelProm': minimum relative prominence
%                   must be a positive scalar
%                   default == 0.02
%                   - 'SpikeHistParams': spike histogram parameters
%                   must be empty or a numeric vector
%                   default == use compute_spike_histogram(spikeTimesMs)
%                   - 'SpikeHistData': spike histogram data
%                   must be empty or a numeric vector
%                   default == use compute_spike_histogram(spikeTimesMs)
%                   - Any other parameter-value pair for the xcorr() function
%
% Requires:
%       cd/compute_spike_histogram.m
%       cd/create_error_for_nargin.m
%       cd/movingaveragefilter.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_network_analyze_spikes.m
%       cd/parse_multiunit.m

% File History:
% 2019-07-23 Moved from parse_multiunit.m
% 2019-08-01 Fixed missing return
% 2020-01-08 Changed default minBurstLengthMs from 10 to 60
% 2020-01-08 Changed default maxInterBurstIntervalMs from 1000 to 2000
% TODO: allow input to be non-event data as well
% 

% Hard-coded constants
MS_PER_S = 1000;

%% Default values for optional arguments
stimStartMsDefault = 0;             % stimulation start is at 0 ms by default
binWidthMsDefault = 10;             % use a bin width of 10 ms by default
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
spikeHistParamsDefault = [];        % set later
spikeHistDataDefault = [];          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'spikeTimesMs', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));


% Add parameter-value pairs to the Input Parser
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
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
addParameter(iP, 'SpikeHistParams', spikeHistParamsDefault, ...
    @(x) validateattributes(x, {'numeric', 'struct'}, {'2d'}));
addParameter(iP, 'SpikeHistData', spikeHistDataDefault, ...
    @(x) validateattributes(x, {'numeric', 'struct'}, {'2d'}));

% Read from the Input Parser
parse(iP, spikeTimesMs, varargin{:});
stimStartMs = iP.Results.StimStartMs;
binWidthMs = iP.Results.BinWidthMs;
minBurstLengthMs = iP.Results.MinBurstLengthMs;
maxFirstInterBurstIntervalMs = iP.Results.MaxFirstInterBurstIntervalMs;
maxInterBurstIntervalMs = iP.Results.MaxInterBurstIntervalMs;
minSpikeRateInBurstHz = iP.Results.MinSpikeRateInBurstHz;
filterWidthMs = iP.Results.FilterWidthMs;
minRelProm = iP.Results.MinRelProm;
spikeHistParams = iP.Results.SpikeHistParams;
spikeHistData = iP.Results.SpikeHistData;

% Keep unmatched arguments for the xcorr() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Compute spike histogram if not provided
if isempty(spikeHistParams) || isempty(spikeHistData)
    [spikeHistParams, spikeHistData] = ...
        compute_spike_histogram(spikeTimesMs, 'StimStartMs', stimStartMs, ...
            'BinWidthMs', binWidthMs, ...
            'MinBurstLengthMs', minBurstLengthMs, ...
            'MaxFirstInterBurstIntervalMs', maxFirstInterBurstIntervalMs, ...
            'MaxInterBurstIntervalMs', maxInterBurstIntervalMs, ...
            'MinSpikeRateInBurstHz', minSpikeRateInBurstHz);
end

% Extract spike histogram parameters
nBins = spikeHistParams.nBins;
halfNBins = spikeHistParams.halfNBins;
oscDurationMs = spikeHistParams.oscDurationMs;

% Extract spike histogram data
spikeCounts = spikeHistData.spikeCounts;

% Compute the bin width in seconds
binWidthSec = binWidthMs ./ MS_PER_S;

% Save parameters
autoCorrParams = spikeHistParams;
autoCorrParams.filterWidthMs = filterWidthMs;
autoCorrParams.minRelProm = minRelProm;

%% Do the job
if isempty(spikeTimesMs)
    autoCorrParams.oscIndex1 = 0;
    autoCorrParams.oscIndex2 = 0;
    autoCorrParams.oscIndex3 = NaN;
    autoCorrParams.oscIndex4 = 0;
    autoCorrParams.oscPeriod1Ms = 0;
    autoCorrParams.oscPeriod2Ms = 0;
    autoCorrParams.minOscPeriod2Bins = 0;
    autoCorrParams.maxOscPeriod2Bins = 0;

    autoCorrData.autoCorr = [];
    autoCorrData.acf = [];
    autoCorrData.acfFiltered = [];
    autoCorrData.acfFilteredOfInterest = [];
    autoCorrData.indPeaks = [];
    autoCorrData.indTroughs = [];
    autoCorrData.ampPeaks = [];
    autoCorrData.ampTroughs = [];
    autoCorrData.halfPeriodsToMultiple = [];
    return
end

% Compute an unnormalized autocorrelogram in Hz^2
autoCorr = xcorr(spikeCounts, 'unbiased', otherArguments{:}) ./ binWidthSec ^ 2;

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

%% Output results
autoCorrParams.oscIndex1 = oscIndex1;
autoCorrParams.oscIndex2 = oscIndex2;
autoCorrParams.oscIndex3 = oscIndex3;
autoCorrParams.oscIndex4 = oscIndex4;
autoCorrParams.oscPeriod1Ms = oscPeriod1Ms;
autoCorrParams.oscPeriod2Ms = oscPeriod2Ms;
autoCorrParams.minOscPeriod2Bins = minOscPeriod2Bins;
autoCorrParams.maxOscPeriod2Bins = maxOscPeriod2Bins;

autoCorrData.autoCorr = autoCorr;
autoCorrData.acf = acf;
autoCorrData.acfFiltered = acfFiltered;
autoCorrData.acfFilteredOfInterest = acfFilteredOfInterest;
autoCorrData.indPeaks = indPeaks;
autoCorrData.indTroughs = indTroughs;
autoCorrData.ampPeaks = ampPeaks;
autoCorrData.ampTroughs = ampTroughs;
autoCorrData.halfPeriodsToMultiple = halfPeriodsToMultiple;

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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%