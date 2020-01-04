function signal2Noise = compute_default_signal2noise(data, varargin)
%% Computes a default signal-to-noise ratio
% Usage: signal2Noise = compute_default_signal2noise(data, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       signal2Noise    - computed signal to noise ratio
%                       specified as a numeric scalar
%
% Arguments:
%       data        - data where each column is a vector of samples
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'siMs': sampling interval(s) in milliseconds
%                   must be empty or a positive vector
%                   default == []
%                   - 'tVecs': time vector(s)
%                   must be empty or 
%                       a numeric array or a cell array of numeric arrays
%                   default == []
%                   - 'IdxStimStart': index of stimulation start
%                   must be empty or a numeric vector
%                   default == 1
%                   - 'IdxDetectStart': index of detection start
%                   must be empty or a numeric vector
%                   default == idxStimStart + minDelaySamples
%                   - 'BaseWindows': baseline window(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   default == start to stimStart
%                   - 'FiltFreq': the cutoff frequency(ies) (Hz or normalized) 
%                                   for the bandpass filter
%                   must be a two-element numeric vector
%                   default == [NaN, NaN]
%                   - 'MinDelayMs': minimum delay from stimulation start in ms
%                   must be empty or a numeric vector
%                   default == artifactLengthMs == 25 ms
%                   - 'RelSnrThres2Max': relative signal to noise threshold
%                                           as a proportion of maximum
%                   must be empty or a numeric vector
%                   default == 0.1
%
% Requires:
%       cd/compute_baseline_noise.m
%       cd/compute_time_window.m
%       cd/count_samples.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/force_matrix.m
%       cd/freqfilter.m
%       cd/match_time_info.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-06-02 Created by Adam Lu
% 2019-07-25 Pull out to its own function
% TODO: use in detect_spikes_multiunit.m

%% Hard-coded constants
MS_PER_S = 1000;

% Must be consistent with detect_spikes_multiunit.m   
artifactLengthMs = 25;
defaultRelSnrThres2Max = 0.1;

%% Default values for optional arguments
siMsDefault = [];              % set later
tVecsDefault = [];              % set later
idxStimStartDefault = 1;        % stim at first time point by default
idxDetectStartDefault = [];     % set later
baseWindowsDefault = [];        % set later
filtFreqDefault = [NaN, NaN];   % no bandpass filter by default 
minDelayMsDefault = [];         % set later
relSnrThres2MaxDefault = [];    % set later

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
addRequired(iP, 'data', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['data must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'siMs', siMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IdxStimStart', idxStimStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'IdxDetectStart', idxDetectStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'BaseWindows', baseWindowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['baseWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'FiltFreq', filtFreqDefault, ...
    @(x) isnumeric(x) && isvector(x) && numel(x) <= 2);
addParameter(iP, 'MinDelayMs', minDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RelSnrThres2Max', relSnrThres2MaxDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
siMs = iP.Results.siMs;
tVecs = iP.Results.tVecs;
idxStimStart = iP.Results.IdxStimStart;
idxDetectStart = iP.Results.IdxDetectStart;
baseWindows = iP.Results.BaseWindows;
filtFreq = iP.Results.FiltFreq;
minDelayMs = iP.Results.MinDelayMs;
relSnrThres2Max = iP.Results.RelSnrThres2Max;

%% Preparation
% Count the number of samples for each vector
nSamples = count_samples(data);

% Match time vector(s) with sampling interval(s) and number(s) of samples
[tVecs, siMs, ~] = match_time_info(tVecs, siMs, nSamples);

% Compute the average sampling interval in ms
siMsAvg = mean(siMs);

% Set default relative signal-2-noise ratio from threshold and maximum
if isempty(relSnrThres2Max)
    relSnrThres2Max = defaultRelSnrThres2Max;
end

% Set default minimum delay in ms
if isempty(minDelayMs)
    minDelayMs = artifactLengthMs;
end

% Find the starting index for detecting a spike
if isempty(idxDetectStart)
    % Compute the minimum delay in samples
    minDelaySamples = round(minDelayMs ./ siMs);

    % Find the starting index for detecting a spike
    idxDetectStart = idxStimStart + minDelaySamples;
end

% Construct default baseline windows
if isempty(baseWindows)
    % Convert to the time of stimulation start
    stimStartMs = extract_elements(tVecs, 'specific', 'Index', idxStimStart);

    % Compute baseline windows
    baseWindows = compute_time_window(tVecs, 'TimeEnd', stimStartMs);
end

%% Do the job
% Force as a matrix
tVecs = force_matrix(tVecs);
data = force_matrix(data);

% Bandpass filter if requested
if ~any(isnan(filtFreq))
    siSeconds = siMsAvg / MS_PER_S;    
    vVecsFilt = freqfilter(data, filtFreq, siSeconds, 'FilterType', 'band');
else
    vVecsFilt = data;
end

% Compute all instantaneous slopes in uV/ms == mV/s
slopes = diff(vVecsFilt) ./ siMsAvg;

% Compute baseline slope noise in mV/s
baseSlopeNoise = compute_baseline_noise(slopes, tVecs(1:(end-1), :), ...
                                        baseWindows);

% Compute the average baseline slope noise
avgBaseSlopeNoise = mean(baseSlopeNoise);

% Extract slopes after detection start
slopesAfterDetectStart = ...
    extract_subvectors(slopes, 'IndexStart', idxDetectStart);

% Compute the average maximum slope after detection start
avgSlopeAfterDetectStart = ...
    mean(extract_elements(slopesAfterDetectStart, 'max'));

% Compute a default signal-to-noise ratio
signal2Noise = 1 + relSnrThres2Max * ...
                (avgSlopeAfterDetectStart / avgBaseSlopeNoise - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%