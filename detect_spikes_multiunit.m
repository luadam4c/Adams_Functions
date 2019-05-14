function [spikesParams, spikesData] = detect_spikes_multiunit(vVec, siMs, varargin)
%% Detects spikes from a multiunit recording
% Usage: [spikesParams, spikesData] = detect_spikes_multiunit(vVec, siMs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       spikesParams- Used and detected parameters, with fields:
%                       siMs
%                       idxStimStart
%                       minDelayMs
%                       maxDelayMs
%                       signal2Noise
%                       baseWindow
%                       minDelaySamples
%                       idxDetectStart
%                       detectStartMs
%                       baseSlopeNoise
%                       slopeThreshold
%                       nSpikesTotal
%                       idxFirstSpike
%                       firstSpikeMs
%                       vMin
%                       vMax
%                       vRange
%                       slopeMin
%                       slopeMax
%                       slopeRange
%                   specified as a scalar structure
%       spikesData  - Detected spikes data, with fields:
%                       slopes
%                       isPeakSlope
%                       isSpike
%                       idxSpikes
%                       spikeTimesMs
%                   specified as a scalar structure
% Arguments:
%       vVec        - voltage vector
%                   must be a TODO
%       siMs        - sampling interval in ms
%                   must be a TODO
%       varargin    - 'FiltFreq': cutoff frequency(ies) (Hz or normalized) 
%                                   for a bandpass filter
%                   must be a numeric a two-element vector
%                   default == none
%                   - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == first half of the trace
%                   - 'IdxStimStart': index of stimulation start
%                   must be a positive scalar
%                   default == 1
%                   - 'MinDelayMs': minimum delay after stim start (ms)
%                   must be a positive scalar
%                   default == 25 ms
%                   - 'MaxDelayMs': maximum delay after stim start (ms)
%                   must be a positive scalar
%                   default == 10000 ms
%                   - 'Signal2Noise': signal-to-noise ratio
%                   must be a positive scalar
%                   default == 3
%                   - 'tVec': original time vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (not used)
%
% Requires:
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/create_logical_array.m
%
% Used by:
%       cd/compute_oscillation_duration.m
%       cd/parse_multiunit.m

% File History:
% 2019-05-03 Moved from parse_multiunit.m
% 2019-05-04 Added input parser
% 2019-05-14 Added 'FiltFreq' as an optional argument
% 2019-05-14 Added 'MaxDelayMs' as an optional argument
% TODO: Finish documentation
% 

%% Hard-coded parameters
MS_PER_S = 1000;

%% Default values for optional arguments
filtFreqDefault = NaN;          % set later
baseWindowDefault = [];         % set later
idxStimStartDefault = 1;    
minDelayMsDefault = 25;
maxDelayMsDefault =10000;       % 1000 ms or 10 seconds
signal2NoiseDefault = 3;        
tVecDefault = [];               % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
addRequired(iP, 'vVec', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FiltFreq', filtFreqDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['BaseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IdxStimStart', idxStimStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MinDelayMs', minDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
addParameter(iP, 'MaxDelayMs', maxDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'tVec', tVecDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vVec, siMs, varargin{:});
filtFreq = iP.Results.FiltFreq;
baseWindow = iP.Results.BaseWindow;
idxStimStart = iP.Results.IdxStimStart;
minDelayMs = iP.Results.MinDelayMs;
maxDelayMs = iP.Results.MaxDelayMs;
signal2Noise = iP.Results.Signal2Noise;
tVec = iP.Results.tVec;

%% Preparation
% Create time vectors
if isempty(tVec)
    nSamples = count_samples(vVec);
    
    tVec = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs);
end

% Compute the minimum delay in samples
minDelaySamples = round(minDelayMs ./ siMs);

% Compute the maximum delay in samples
maxDelaySamples = round(maxDelayMs ./ siMs);

% Find the starting index for detecting a spike
idxDetectStart = idxStimStart + minDelaySamples;

% Find the corresponding time
detectStartMs = tVec(idxDetectStart);

% Compute the number of samples
nSamples = numel(vVec);

% Bandpass filter if requested
if ~isnan(filtFreq)
    siSeconds = siMs / MS_PER_S;    
    vVecFilt = freqfilter(vVec, filtFreq, siSeconds, 'FilterType', 'band');
else
    vVecFilt = vVec;
end

%% Do the job
% Compute all instantaneous slopes in uV/ms == mV/s
slopes = diff(vVecFilt) ./ siMs;

% Compute a baseline slope noise in mV/s
baseSlopeNoise = compute_baseline_noise(slopes, tVec(1:(end-1)), baseWindow);

% Compute a slope threshold in mV/s
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

% Find the last index in the range of interest
idxEndOfInterest1 = min(numel(vVecFilt), idxDetectStart + maxDelaySamples);

% Query the maximum and range of vVec after detectStartMs
vVecTrunc = vVecFilt(idxDetectStart:idxEndOfInterest1);
vMin = min(vVecTrunc);
vMax = max(vVecTrunc);
vRange = vMax - vMin;

% Find the last index in the range of interest
idxEndOfInterest2 = min(numel(slopes), idxDetectStart + maxDelaySamples);

% Query the maximum and range of slope after detectStartMs
slopesTrunc = slopes(idxDetectStart:idxEndOfInterest2);
slopeMin = min(slopesTrunc);
slopeMax = max(slopesTrunc);
slopeRange = slopeMax - slopeMin;

%% Save in output
spikesParams.siMs = siMs;
spikesParams.idxStimStart = idxStimStart;
spikesParams.minDelayMs = minDelayMs;
spikesParams.maxDelayMs = maxDelayMs;
spikesParams.filtFreq = filtFreq;
spikesParams.signal2Noise = signal2Noise;
spikesParams.baseWindow = baseWindow;
spikesParams.minDelaySamples = minDelaySamples;
spikesParams.idxDetectStart = idxDetectStart;
spikesParams.detectStartMs = detectStartMs;
spikesParams.baseSlopeNoise = baseSlopeNoise;
spikesParams.slopeThreshold = slopeThreshold;
spikesParams.nSpikesTotal = nSpikesTotal;
spikesParams.idxFirstSpike = idxFirstSpike;
spikesParams.firstSpikeMs = firstSpikeMs;
spikesParams.vMin = vMin;
spikesParams.vMax = vMax;
spikesParams.vRange = vRange;
spikesParams.slopeMin = slopeMin;
spikesParams.slopeMax = slopeMax;
spikesParams.slopeRange = slopeRange;

spikesData.vVecFilt = vVecFilt;
spikesData.slopes = slopes;
spikesData.isPeakSlope = isPeakSlope;
spikesData.isSpike = isSpike;
spikesData.idxSpikes = idxSpikes;
spikesData.spikeTimesMs = spikeTimesMs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%