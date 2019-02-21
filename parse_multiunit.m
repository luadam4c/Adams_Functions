function varargout = parse_multiunit (vVecs, siMs, varargin)
%% Parses multiunit recordings: detect spikes
% Usage: [parsedParams, parsedData] = parse_multiunit (vVecs, siMs, varargin)
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
%       varargin    - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == detect from pulse vector
%                   - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   
% Requires:
%       cd/count_samples.m TODO
%       cd/count_vectors.m TODO
%       cd/iscellnumeric.m TODO
%       cd/find_stim_start.m TODO
%       cd/create_logical_array.m
%       cd/extract_elements.m
%       cd/force_column_cell.m
%       cd/compute_baseline_noise.m
%       cd/create_error_for_nargin.m
%       cd/match_time_info.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-02-19 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
stimStartMsDefault = [];        % set later
pulseVectorsDefault = [];       % don't use pulse vectors by default
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
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PulseVectors', pulseVectorsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['PulseVectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vVecs, siMs, varargin{:});
stimStartMs = iP.Results.StimStartMs;
pulseVectors = iP.Results.PulseVectors;
tVecs = iP.Results.tVecs;

%% Preparation
% Count the number of vectors
nVectors = count_vectors(vVecs);

% Count the number of samples for each vector
nSamples = count_samples(vVecs);

% Match time vector(s) with sampling interval(s) and number(s) of samples
[tVecs, siMs, nSamples] = match_time_info(tVecs, siMs, nSamples);

%% Do the job
% Detect stimulation start time if not provided
%   Otherwise find the corresponding index in the time vector
if isempty(stimStartMs)
    % TODO: Make this a function find_stim_start.m
    if ~isempty(pulseVectors)
        % Parse the pulse vectors
        [pulseParams, pulseData] = ...
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
if isempty(baseWindows)
    % Get the starting time(s)
    timeStartMs = extract_elements(tVecs, 'first');

    % Use timeStartMs to stimStartMs by default
    baseWindows = transpose([timeStartMs, stimStartMs]);
end

% Compute baseline rms noise from window
baseNoises = compute_baseline_noise(vVecs, tVecs, baseWindows);

% Force as a cell array of vectors
vVecs = force_column_cell(vVecs);

% Parse all of them in a parfor loop
parsedParamsCell = cell(nVectors, 1);
parsedDataCell = cell(nVectors, 1);
%parfor iVec = 1:nVectors
for iVec = 1:1
    [parsedParamsCell{iVec}, parsedDataCell{iVec}] = ...
        parse_multiunit_helper(vVecs{iVec}, siMs(iVec), ...
                                idxStimStart(iVec), baseNoises(iVec));
end

% Convert to a struct array
%   Note: This removes all entries that are empty
[parsedParamsStruct, parsedDataStruct] = ...
    argfun(@(x) [x{:}], parsedParamsCell, parsedDataCell);

% Convert to a table
[parsedParams, parsedData] = ...
    argfun(@struct2table, parsedParamsStruct, parsedDataStruct);

%% Outputs
varargout{1} = parsedParams;
varargout{2} = parsedData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parsedParams, parsedData] = ...
                parse_multiunit_helper(vVec, siMs, idxStimStart, baseNoise)
% Parse a single multiunit recording

% Hard-coded parameters
signal2Noise = 3; %2
minDelaySamples = 2000;

% Compute a slope threshold
slopeThreshold = (baseNoise / siMs) * signal2Noise;

% Compute the number of samples
nSamples = numel(vVec);

% Compute all instantaneous slopes
slopes = diff(vVec) / siMs;

% Determine whether each slope is a local maximum
[~, indPeakSlopes] = findpeaks(slopes);
isPeakSlope = create_logical_array(indPeakSlopes, [nSamples, 1]);

% Create all indices minus 1
allIndices = transpose(1:nSamples);

% Detect spikes after idxStimStart + minDelaySamples
isSpike = [false; slopes > slopeThreshold] & ...
            [false; isPeakSlope] & ...
            allIndices > idxStimStart + minDelaySamples;
idxSpikes = find(isSpike);

figure(1); 
clf; hold on
plot(vVec, 'k');
plot(idxSpikes, vVec(idxSpikes), 'rx');
xlim([34100, 34700])
%xlim([idxStimStart + minDelaySamples, idxStimStart + 1e4])

figure(2); 
clf; hold on
plot(slopes, 'k');
plot_horizontal_line(slopeThreshold, 'Color', 'b', 'LineStyle', '--');
plot(idxSpikes, slopes(idxSpikes - 1), 'rx');
xlim([34100, 34700])
%xlim([idxStimStart + minDelaySamples, idxStimStart + 1e4])

% Sufficiently steep positive deflections


% Sufficiently steep negative deflections

% Store in outputs
parsedParams.signal2Noise = signal2Noise;
parsedParams.minDelaySamples = minDelaySamples;
parsedParams.idxStimStart = idxStimStart;
parsedParams.baseNoise = baseNoise;
parsedParams.slopeThreshold = slopeThreshold;
parsedData.vVec = vVec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%