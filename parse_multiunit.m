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
%                   default == find_first_jump(vVec0s)
%                   - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'tVecs': original time vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%                   
% Requires:
%       cd/count_vectors.m TODO
%       cd/iscellnumeric.m TODO
%       cd/find_stim_start.m TODO
%       cd/extract_elements.m TODO
%       cd/compute_baseline_noise.m TODO
%       cd/create_error_for_nargin.m
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

% TODO: Pull this out to create_default_time_info.m and use this in parse_lts.m too
%   [tVecs, siMs] = create_default_time_info(tVecs, siMs);
% Compute sampling interval(s) and create time vector(s)
if isempty(siMs) && ~isempty(tVecs)
    % Compute sampling interval(s) in ms
    siMs = compute_sampling_interval(tVecs);
elseif isempty(tVecs) && ~isempty(siMs)
    % Create time vector(s)
    tVecs = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                                    'TimeUnits', 'ms');
elseif isempty(tVecs) && isempty(siMs)
    error('One of siMs and tVec0s must be provided!');
end

% Make the number of vectors and 
[tVecs, siMs] = match_dimensions(tVecs, siMs);

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
        stimStartMs = extract_elements(tVecs, 'specific', 'Index', idxStimStart);
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
    baseWindows = [timeStartMs, stimStartMs];
end

% Compute baseline rms noise from window
baseNoises = compute_baseline_noise(vVec, tVecs, baseWindows);

% Parse all of them in a parfor loop
parsedParamsCell = cell(nVectors, 1);
parsedDataCell = cell(nVectors, 1);
%parfor iVec = 1:nVectors
for iVec = 1:nVectors
    [parsedParamsCell{iVec}, parsedDataCell{iVec}] = ...
        parse_multiunit_helper(vVecs{iVec}, idxStimStart(iVec), baseNoises(iVec));
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
                parse_multiunit_helper(vVec, idxStimStart, baseNoise)
% Parse a single multiunit recording

signal2Noise = 2; %3

slopeThreshold = baseNoise * signal2Noise;

% Sufficiently steep positive deflections

% Sufficiently steep negative deflections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%