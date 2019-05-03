function [spikesParams, spikesData] = ...
                detect_spikes_multiunit(vVec, tVec, ...
                                    siMs, idxStimStart, minDelaySamples, ...
                                    signal2Noise, baseWindow);
%% Detects spikes from a multiunit recording
% Usage: [spikesParams, spikesData] = ...
%               detect_spikes_multiunit(vVec, tVec, ...
%                                   siMs, idxStimStart, minDelaySamples, ...
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
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/create_logical_array.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-05-03 Created by Adam Lu
% TODO: Input parser
% TODO: Documentation
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, reqarg1, varargin{:});
param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% TODO

%% Do the job
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