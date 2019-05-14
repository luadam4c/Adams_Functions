function [stimParams, pulseParams, pulseData] = parse_stim (pulseVecs, varargin)
%% Detects the index and time of stimulation start from pulse vectors
% Usage: [stimParams, pulseParams, pulseData] = parse_stim (pulseVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       stimParams  - stimulation parameters, with fields:
%                       idxStimStart
%                       stimStartMs
%                       any other fields returned by parse_pulse.m
%                   specified as a scalar structure
%       pulseParams - see parse_pulse.m
%                   specified as a table
%       pulseData   - see parse_pulse.m
%                   specified as a table
% Arguments:
%       pulseVecs   - vectors containing a stimulation pulse
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'SamplingIntervalMs': sampling interval(s) in ms
%                   must be a positive vector
%                   default == []
%                   - 'StimStartMs': time of stimulation start (ms)
%                   must be a positive scalar
%                   default == detect from pulse vector
%                   - Any other parameter-value pair for 
%                       the parse_pulse() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/extract_elements.m
%       cd/find_closest.m
%       cd/merge_structs.m
%       cd/parse_pulse.m
%
% Used by:
%       cd/compute_oscillation_duration.m
%       cd/parse_multiunit.m

% File History:
% 2019-05-14 Pulled from parse_multiunit.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
samplingIntervalMsDefault = 0.1;% 0.1 ms by default
tVecsDefault = [];              % no time info by default
stimStartMsDefault = [];        % set later

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
addRequired(iP, 'pulseVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingIntervalMs', samplingIntervalMsDefault, ...
    @(x) isempty(x) || ispositivevector(x));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, pulseVecs, varargin{:});
siMs = iP.Results.SamplingIntervalMs;
tVecs = iP.Results.tVecs;
stimStartMs = iP.Results.StimStartMs;

% Keep unmatched arguments for the parse_pulse() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation

%% Do the job
% If stimStartMs already provided, use tVecs or siMs to find the index
if ~isempty(stimStartMs)
    % Find the indices of stimulation start
    if ~isempty(tVecs)
        % Use the indices of tVecs with values closest to stimStartMs
        % TODO: find_closest.m
        idxStimStart = find_closest(tVecs, stimStartMs);
    else
        % Assume tVecs start from 0 and use siMs
        idxStimStart = round(stimStartMs ./ siMs);
    end

    stimParams.idxStimStart = idxStimStart;
    stimParams.stimStartMs = stimStartMs;
    pulseParams = [];
    pulseData = [];
    return
end

% Detect stimulation start from the pulse vectors
if ~isempty(pulseVecs)
    % Parse the pulse vectors
    [pulseParams, pulseData] = ...
        parse_pulse(pulseVecs, 'SamplingIntervalMs', siMs, ...
                    otherArguments{:});

    % Use the indices after pulse starts for stimulation start
    idxStimStart = pulseParams{:, 'idxAfterStart'};

    % Use the time vectors
    if ~isempty(tVecs)
        stimStartMs = extract_elements(tVecs, 'specific', ...
                                        'Index', idxStimStart);
    else
        stimStartMs = idxStimStart * siMs;
    end
else
    error('One of stimStartMs and pulseVecs must be provided!');
end

%% Output results
stimParams.idxStimStart = idxStimStart;
stimParams.stimStartMs = stimStartMs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
