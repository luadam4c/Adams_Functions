function [parsedParams, parsedData] = parse_pulse (vectors, varargin)
%% Parses pulse widths, endpoints, amplitudes for vector(s) containing a single pulse
% Usage: [parsedParams, parsedData] = parse_pulse (vectors, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%       
% Outputs:
%       parsedParams    - a table containing the parsed parameters, including:
%                           nSamples
%                           pulseWidthSamples
%                           idxBeforeStart
%                           idxBeforeEnd
%                           idxAfterStart
%                           idxAfterEnd
%                           idxMidpoint
%                           baseValue
%                           pulseValue
%                           pulseAmplitude
%                       if siMs or timeVecs is provided, also:
%                           pulseWidthMs
%                           timeBeforeStartMs
%                           timeBeforeEndMs
%                           timeAfterStartMs
%                           timeAfterEndMs
%                           timeMidpointMs
%                       specified as a table
%       parsedData      - a table containing the parsed data, including:
%                           vectors
%                           indBase
%                           indPulse
%                       specified as a table
% Arguments:    
%       vectors     - vectors containing a pulse
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'SamplingIntervalMs': sampling interval(s) in ms
%                   Note: Only one of 'SamplingIntervalMs' or 'TimeVecs'
%                         can be provided.
%                   must be a positive vector
%                   default == []
%                   - 'TimeVecs': corresponding time vectors
%                   Note: Only one of 'SamplingIntervalMs' or 'TimeVecs'
%                         can be provided.
%                   must be a numeric array or a cell array of numeric vectors
%                   default == []
%                   - 'MinPulseAmplitude': minimum pulse amplitude
%                   must be a non-negative scalar
%                   default == 0
%
% Requires:
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/find_pulse_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/ispositivevector.m
%       cd/match_format_vector_sets.m
%
% Used by:    
%       cd/parse_stim.m
%       cd/find_passive_params.m
%       cd/identify_repetitive_pulses.m
%       cd/parse_pulse_response.m
%       cd/plot_pulse.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-12-15 Now allows vectors to have no pulse (will return NaNs)
% 2018-12-17 Now computes times if SamplingIntervalMs is provided
% 2025-09-02 Modified to accept 'TimeVecs' as an optional parameter
% 

%% Hard-coded parameters

%% Default values for optional arguments
samplingIntervalMsDefault = []; % no time information by default
timeVecsDefault = [];           % no time vectors by default
minPulseAmplitudeDefault = 0;

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SamplingIntervalMs', samplingIntervalMsDefault, ...
    @(x) isempty(x) || ispositivevector(x));
addParameter(iP, 'TimeVecs', timeVecsDefault, ...
    @(x) isempty(x) || isnumeric(x) || iscellnumeric(x));
addParameter(iP, 'MinPulseAmplitude', minPulseAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
siMs = iP.Results.SamplingIntervalMs;
timeVecs = iP.Results.TimeVecs;
minPulseAmplitude = iP.Results.MinPulseAmplitude;

% Check that 'SamplingIntervalMs' and 'TimeVecs' are not both provided
if ~isempty(siMs) && ~isempty(timeVecs)
    error('Usage error: ''SamplingIntervalMs'' and ''TimeVecs'' cannot both be provided.');
end

%% Preparation
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);

% If time vectors provided, match the vector sets
if ~isempty(timeVecs)
    % Force to a column cell array
    timeVecs = force_column_cell(timeVecs);

    % Match the number of vectors
    [timeVecs, vectors] = match_format_vector_sets(timeVecs, vectors);
end

% Count the number of samples for each vector
nSamples = count_samples(vectors);

%% Do the job
% Find indices for all the pulse endpoints
[idxBeforeStart, idxBeforeEnd, idxAfterStart, idxAfterEnd] = ...
    find_pulse_endpoints(vectors, 'MinPulseAmplitude', minPulseAmplitude);

% Find indices for all the pulse midpoints
idxMidpoint = round((idxAfterStart + idxBeforeEnd) ./ 2);

% Find all pulse widths in samples
pulseWidthSamples = idxBeforeEnd - idxBeforeStart;

% Find the baseline indices
indBase = arrayfun(@(x, y, z) find_baseline_indices(x, y, z), ...
                    idxBeforeStart, idxAfterEnd, nSamples, ...
                    'UniformOutput', false);

% Find the average baseline value (could be NaN)
baseValue = cellfun(@(x, y) mean(x(y)), vectors, indBase);

% Find the pulse indices
indPulse = arrayfun(@(x, y) find_pulse_indices(x, y), ...
                    idxAfterStart, idxBeforeEnd, ...
                    'UniformOutput', false);

% Find the average pulse value (could be NaN)
pulseValue = compute_stats(vectors, 'mean', 'Indices', indPulse);

% Find the pulse amplitudes
pulseAmplitude = pulseValue - baseValue;

%% Store results in output
parsedParams = table(nSamples, pulseWidthSamples, ...
                    idxBeforeStart, idxBeforeEnd, ...
                    idxAfterStart, idxAfterEnd, idxMidpoint, ...
                    baseValue, pulseValue, pulseAmplitude);
parsedData = table(vectors, indBase, indPulse);

%% Optional extra parameters
% If time information is provided, add the corresponding times
if ~isempty(timeVecs)
    % Compute the corresponding times from the time vectors
    timeBeforeStartMs = cellfun(@(tv, idx) get_time_at_idx(tv, idx), timeVecs, num2cell(idxBeforeStart));
    timeBeforeEndMs   = cellfun(@(tv, idx) get_time_at_idx(tv, idx), timeVecs, num2cell(idxBeforeEnd));
    timeAfterStartMs  = cellfun(@(tv, idx) get_time_at_idx(tv, idx), timeVecs, num2cell(idxAfterStart));
    timeAfterEndMs    = cellfun(@(tv, idx) get_time_at_idx(tv, idx), timeVecs, num2cell(idxAfterEnd));
    timeMidpointMs    = cellfun(@(tv, idx) get_time_at_idx(tv, idx), timeVecs, num2cell(idxMidpoint));
    pulseWidthMs      = timeBeforeEndMs - timeBeforeStartMs;
elseif ~isempty(siMs)
    % Compute the corresponding times from the sampling interval
    % TODO: Use convert_to_time.m
    % TODO: Does the time make sense if doesn't start from zero?
    [pulseWidthMs, timeBeforeStartMs, timeBeforeEndMs, ...
        timeAfterStartMs, timeAfterEndMs, timeMidpointMs] = ...
        argfun(@(x) x .* siMs, ...
            pulseWidthSamples, idxBeforeStart, idxBeforeEnd, ...
            idxAfterStart, idxAfterEnd, idxMidpoint);
end

if ~isempty(timeVecs) || ~isempty(siMs)
    % Place in table
    timeParams = table(pulseWidthMs, timeBeforeStartMs, timeBeforeEndMs, ...
                        timeAfterStartMs, timeAfterEndMs, timeMidpointMs);

    % Join to parsedParams
    parsedParams = horzcat(parsedParams, timeParams);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indBase = find_baseline_indices(idxBeforeStart, idxAfterEnd, nSamples)
%% Returns the baseline indices

% Decide on the baseline indices based on whether a pulse was found
if isnan(idxBeforeStart) || isnan(idxAfterEnd) || isnan(nSamples)
    % No pulse is found, the entire trace is baseline
    indBase = 1:nSamples;
else
    % Everything outside the pulse is baseline
    indBase = [1:idxBeforeStart, idxAfterEnd:nSamples];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indPulse = find_pulse_indices(idxAfterStart, idxBeforeEnd)
%% Returns the pulse indices

% Decide on the pulse indices based on whether a pulse was found
if isnan(idxAfterStart) || isnan(idxBeforeEnd)
    % No pulse is found
    indPulse = [];
else
    % Use the restrictive pulse endpoints
    indPulse = idxAfterStart:idxBeforeEnd;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeValue = get_time_at_idx(timeVector, index)
%% Safely gets a time value from a time vector at a specific index
if isnan(index) || isempty(index) || index < 1 || index > numel(timeVector)
    timeValue = NaN;
else
    timeValue = timeVector(index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

pulseValue = cellfun(@(x, y) mean(x(y)), vectors, indPulse);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
