function [spikesParams, spikesData] = ...
                detect_spikes_multiunit(vVec, siMs, ...
                                    tVec, idxStimStart, minDelaySamples, ...
                                    signal2Noise, baseWindow);
%% Detects spikes from a multiunit recording
% Usage: [spikesParams, spikesData] = ...
%               detect_spikes_multiunit(vVec, siMs, ...
%                                   tVec, idxStimStart, minDelaySamples, ...
%                                   signal2Noise, baseWindow);
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       spikesParams- Used and detected parameters, with fields:
%                       siMs
%                       idxStimStart
%                       minDelaySamples
%                       signal2Noise
%                       baseWindow
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
%       tVec        - time vector
%                   must be a TODO
%       vVec        - voltage vector
%                   must be a TODO
%       varargin    - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == first half of the trace
%
% Requires:
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/create_logical_array.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-05-03 Moved from parse_multiunit.m
% 2019-05-04 Added input parser
% TODO: Documentation
% 

%% Hard-coded parameters

%% Default values for optional arguments
tVecDefault = [];               % set later
idxStimStartDefault = 1;    
minDelaySamplesDefault = 0;
signal2NoiseDefault = 3;        
baseWindowDefault = [];         % set later

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
addRequired(iP, 'reqarg1');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'tVec', tVecDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IdxStimStart', idxStimStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MinDelaySamples', minDelaySamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['BaseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vVec, siMs, varargin{:});
tVec = iP.Results.tVec;
idxStimStart = iP.Results.IdxStimStart;
minDelaySamples = iP.Results.MinDelaySamples;
signal2Noise = iP.Results.Signal2Noise;
baseWindow = iP.Results.BaseWindow;

%% Preparation
% Create time vectors
if isempty(tVec)
    tVec = create_time_vectors(nSamples);
end

% Find the starting index for detecting a spike
idxDetectStart = idxStimStart + minDelaySamples;

% Find the corresponding time
detectStartMs = tVec(idxDetectStart);

% Compute the number of samples
nSamples = numel(vVec);

%% Do the job
% Compute all instantaneous slopes in uV/ms == mV/s
slopes = diff(vVec) ./ siMs;

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

%% Save in output
spikesParams.siMs = siMs;
spikesParams.idxStimStart = idxStimStart;
spikesParams.minDelaySamples = minDelaySamples;
spikesParams.signal2Noise = signal2Noise;
spikesParams.baseWindow = baseWindow;
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