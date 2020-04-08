function [histParams, histData] = compute_spike_histogram (spikeTimesMs, varargin)
%% Computes a spike histogram, detect bursts and compute the oscillation duration
% Usage: [histParams, histData] = compute_spike_histogram (spikeTimesMs, varargin)
% Explanation:
%       1. Detected spikes are grouped into 10-ms bins
%       2. 'Bursts' are detected by setting 
%               i. a minimum burst length (minBurstLengthMs) 
%           and ii. a minimum spike rate within a burst (minSpikeRateInBurstHz)
%       3. 'Oscillation duration' is defined by the maximum contiguous stretch 
%           of bursts starting from stimulation start, 
%           with inter-burst intervals no greater than maxInterBurstIntervalMs.
%
% Example(s):
%       [histParams, histData] = compute_spike_histogram(1000*rand(100, 1));
%
% Outputs:
%       histParams  - spike histogram parameters, with fields:
%                       stimStartMs
%                       binWidthMs
%                       binWidthSec
%                       minBurstLengthMs
%                       maxInterBurstIntervalMs
%                       minSpikeRateInBurstHz
%                       nBins
%                       halfNBins
%                       histLeftMs
%                       nBurstsTotal
%                       nBurstsIn10s
%                       nBurstsInOsc
%                       iBinLastOfLastBurst
%                       iBinLastOfLastBurstIn10s
%                       iBinLastOfLastBurstInOsc
%                       nSpikesPerBurst
%                       nSpikesPerBurstIn10s
%                       nSpikesPerBurstInOsc
%                       nSpikesIn10s
%                       nSpikesInOsc
%                       timeOscEndMs
%                       oscDurationMs
%                   specified as a scalar structure
%       histData    - spike histogram data, with fields:
%                       spikeCounts
%                       edgesMs
%                       iBinBurstStarts
%                       iBinBurstEnds
%                       iBinBurstIn10sStarts
%                       iBinBurstIn10sEnds
%                       iBinBurstInOscStarts
%                       iBinBurstInOscEnds
%                       spikeCountsEachBurst
%                       spikeCountsEachBurstIn10s
%                       spikeCountsEachBurstInOsc
%                       timeBurstStartsMs
%                       timeBurstEndsMs
%                       timeBurstIn10sStartsMs
%                       timeBurstIn10sEndsMs
%                       timeBurstInOscStartsMs
%                       timeBurstInOscEndsMs
%                   specified as a scalar structure
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
%                   - Any other parameter-value pair for 
%                       the histcounts() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/compute_bins.m
%
% Used by:
%       cd/compute_autocorrelogram.m
%       cd/compute_oscillation_duration.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/parse_multiunit.m

% File History:
% 2019-05-13 Pulled from parse_multiunit.m
% 2019-08-08 Added maxFirstInterBurstIntervalMs
% 2019-11-13 Added stimStartMs as a fixed edge for compute_bins.m
% 2019-11-13 Now does not consider an oscillation evoked 
%               if the starting bin is more than maxFirstIbiBins away 
%               from stimulation start
% 2020-01-08 Now saves all the fields required in plot_spike_histogram.m
% 2020-01-08 Changed default minBurstLengthMs from 10 to 60
% 2020-01-08 Changed default maxInterBurstIntervalMs from 1000 to 2000

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

% Read from the Input Parser
parse(iP, spikeTimesMs, varargin{:});
stimStartMs = iP.Results.StimStartMs;
binWidthMs = iP.Results.BinWidthMs;
minBurstLengthMs = iP.Results.MinBurstLengthMs;
maxFirstInterBurstIntervalMs = iP.Results.MaxFirstInterBurstIntervalMs;
maxInterBurstIntervalMs = iP.Results.MaxInterBurstIntervalMs;
minSpikeRateInBurstHz = iP.Results.MinSpikeRateInBurstHz;

% Keep unmatched arguments for the histcounts() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Compute the bin width in seconds
[stimStartSec, binWidthSec, maxInterBurstIntervalSec] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            stimStartMs, binWidthMs, maxInterBurstIntervalMs);

% Save parameters used
histParams.stimStartMs = stimStartMs;
histParams.stimStartSec = stimStartSec;
histParams.binWidthMs = binWidthMs;
histParams.binWidthSec = binWidthSec;
histParams.minBurstLengthMs = minBurstLengthMs;
histParams.maxFirstInterBurstIntervalMs = maxFirstInterBurstIntervalMs;
histParams.maxInterBurstIntervalMs = maxInterBurstIntervalMs;
histParams.maxInterBurstIntervalSec = maxInterBurstIntervalSec;
histParams.minSpikeRateInBurstHz = minSpikeRateInBurstHz;

%% Do the job
if isempty(spikeTimesMs)
    histParams.nBins = 0;
    histParams.halfNBins = 0;
    histParams.histLeftMs = NaN;
    histParams.histLeftSec = NaN;
    histParams.nBurstsTotal = 0;
    histParams.nBurstsIn10s = 0;
    histParams.nBurstsInOsc = 0;
    histParams.iBinLastOfLastBurst = NaN;
    histParams.iBinLastOfLastBurstIn10s = NaN;
    histParams.iBinLastOfLastBurstInOsc = NaN;
    histParams.nSpikesPerBurst = NaN;
    histParams.nSpikesPerBurstIn10s = NaN;
    histParams.nSpikesPerBurstInOsc = NaN;
    histParams.nSpikesIn10s = 0;
    histParams.nSpikesInOsc = 0;
    histParams.timeOscEndMs = stimStartMs;
    histParams.timeOscEndSec = stimStartSec;
    histParams.oscDurationMs = 0;    
    histParams.oscDurationSec = 0;

    histData.spikeCounts = [];
    histData.edgesMs = [];
    histData.edgesSec = [];
    histData.iBinBurstStarts = [];
    histData.iBinBurstEnds = [];
    histData.iBinBurstIn10sStarts = [];
    histData.iBinBurstIn10sEnds = [];
    histData.iBinBurstInOscStarts = [];
    histData.iBinBurstInOscEnds = [];
    histData.spikeCountsEachBurst = [];
    histData.spikeCountsEachBurstIn10s = [];
    histData.spikeCountsEachBurstInOsc = [];
    histData.timeBurstStartsMs = [];
    histData.timeBurstStartsSec = [];
    histData.timeBurstEndsMs = [];
    histData.timeBurstEndsSec = [];
    histData.timeBurstIn10sStartsMs = [];
    histData.timeBurstIn10sEndsMs = [];
    histData.timeBurstInOscStartsMs = [];
    histData.timeBurstInOscEndsMs = [];
    return
end

% Compute a spike histogram
[spikeCounts, edgesMs] = compute_bins(spikeTimesMs, 'BinWidth', binWidthMs, ...
                                'FixedEdges', stimStartMs, otherArguments{:});

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

% Compute the minimum number of spikes per sliding window if in a burst
minSpikesPerWindowInBurst = ceil(minSpikeRateInBurstHz * slidingWinSec);

% Compute the maximum number of bins between the stimulation start 
%   and the first burst
maxFirstIbiBins = floor(maxFirstInterBurstIntervalMs ./ binWidthMs);

% Compute the maximum number of bins between subsequent consecutive bursts
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

    % Count the number of inter-burst intervals
    nIBIs = numel(ibiBins);

    % Find the first inter-burst interval greater than maxIbiBins
    firstIbiTooLarge = find(ibiBins > maxIbiBins, 1, 'first');

    % If this happens to be the first inter-burst interval, 
    %   and that it is not greater than maxFirstIbiBins, search again
    if ~isempty(firstIbiTooLarge) && firstIbiTooLarge == 1 && ...
            ibiBins(1) <= maxFirstIbiBins
        if nIBIs < 2
            firstIbiTooLarge = [];
        else
            % Find the first inter-burst interval in the rest of the bins
            %   than is greater than maxIbiBins
            firstIbiTooLargeShifted = ...
                find(ibiBins(2:end) > maxIbiBins, 1, 'first');

            % Add back the first bin so that it is the correct bin count
            firstIbiTooLarge = firstIbiTooLargeShifted + 1;
        end
    end

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

    % If the starting bin is more than maxFirstIbiBins away from stimulation
    %   start, don't consider it an evoked oscillation
    if ~isempty(iBinBurstInOscStarts) && ...
            iBinBurstInOscStarts(1) > maxFirstIbiBins
        nBurstsInOsc = 0;
        iBinBurstInOscStarts = [];
        iBinBurstInOscEnds = [];
    end
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

% Convert to seconds
[histLeftSec, timeOscEndSec, oscDurationSec, ...
    edgesSec, timeBurstStartsSec, timeBurstEndsSec] = ...
    argfun(@(x) x ./ MS_PER_S, ...
            histLeftMs, timeOscEndMs, oscDurationMs, ...
            edgesMs, timeBurstStartsMs, timeBurstEndsMs);

%% Output results
histParams.nBins = nBins;
histParams.halfNBins = halfNBins;
histParams.histLeftMs = histLeftMs;
histParams.histLeftSec = histLeftSec;
histParams.nBurstsTotal = nBurstsTotal;
histParams.nBurstsIn10s = nBurstsIn10s;
histParams.nBurstsInOsc = nBurstsInOsc;
histParams.iBinLastOfLastBurst = iBinLastOfLastBurst;
histParams.iBinLastOfLastBurstIn10s = iBinLastOfLastBurstIn10s;
histParams.iBinLastOfLastBurstInOsc = iBinLastOfLastBurstInOsc;
histParams.nSpikesPerBurst = nSpikesPerBurst;
histParams.nSpikesPerBurstIn10s = nSpikesPerBurstIn10s;
histParams.nSpikesPerBurstInOsc = nSpikesPerBurstInOsc;
histParams.nSpikesIn10s = nSpikesIn10s;
histParams.nSpikesInOsc = nSpikesInOsc;
histParams.timeOscEndMs = timeOscEndMs;
histParams.timeOscEndSec = timeOscEndSec;
histParams.oscDurationMs = oscDurationMs;
histParams.oscDurationSec = oscDurationSec;

histData.spikeCounts = spikeCounts;
histData.edgesMs = edgesMs;
histData.edgesSec = edgesSec;
histData.iBinBurstStarts = iBinBurstStarts;
histData.iBinBurstEnds = iBinBurstEnds;
histData.iBinBurstIn10sStarts = iBinBurstIn10sStarts;
histData.iBinBurstIn10sEnds = iBinBurstIn10sEnds;
histData.iBinBurstInOscStarts = iBinBurstInOscStarts;
histData.iBinBurstInOscEnds = iBinBurstInOscEnds;
histData.spikeCountsEachBurst = spikeCountsEachBurst;
histData.spikeCountsEachBurstIn10s = spikeCountsEachBurstIn10s;
histData.spikeCountsEachBurstInOsc = spikeCountsEachBurstInOsc;
histData.timeBurstStartsMs = timeBurstStartsMs;
histData.timeBurstStartsSec = timeBurstStartsSec;
histData.timeBurstEndsMs = timeBurstEndsMs;
histData.timeBurstEndsSec = timeBurstEndsSec;
histData.timeBurstIn10sStartsMs = timeBurstIn10sStartsMs;
histData.timeBurstIn10sEndsMs = timeBurstIn10sEndsMs;
histData.timeBurstInOscStartsMs = timeBurstInOscStartsMs;
histData.timeBurstInOscEndsMs = timeBurstInOscEndsMs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iBinLastOfLastBurst, spikeCountsEachBurst, ...
                nSpikesPerBurst] = ...
                compute_burst_statistics (spikeCounts, ...
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

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%