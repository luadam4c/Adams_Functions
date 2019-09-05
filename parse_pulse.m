function [parsedParams, parsedData] = parse_pulse (vectors, varargin)
%% Parses pulse widths, endpoints, amplitudes for vector(s) containing a pulse
% Usage: [parsedParams, parsedData] = parse_pulse (vectors, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
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
%                       if siMs is provided, also:
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
%                   must be a positive vector
%                   default == []
%
% Requires:
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/find_pulse_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/ispositivevector.m
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
% 

%% Hard-coded parameters

%% Default values for optional arguments
samplingIntervalMsDefault = []; % no time information by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
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

% Read from the Input Parser
parse(iP, vectors, varargin{:});
siMs = iP.Results.SamplingIntervalMs;

%% Preparation
% Force vectors to be a column cell array
vectors = force_column_cell(vectors);

% Count the number of samples for each vector
nSamples = count_samples(vectors);

%% Do the job
% Find indices for all the pulse endpoints
[idxBeforeStart, idxBeforeEnd, idxAfterStart, idxAfterEnd] = ...
    find_pulse_endpoints(vectors);

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
% If siMs provided, add the corresponding times
if ~isempty(siMs)
    % Compute the corresponding times
    % TODO: Use convert_to_time.m
    % TODO: Does the time make sense if doesn't start from zero?
    [pulseWidthMs, timeBeforeStartMs, timeBeforeEndMs, ...
        timeAfterStartMs, timeAfterEndMs, timeMidpointMs] = ...
        argfun(@(x) x .* siMs, ...
            pulseWidthSamples, idxBeforeStart, idxBeforeEnd, ...
            idxAfterStart, idxAfterEnd, idxMidpoint);

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

%{
OLD CODE:

pulseValue = cellfun(@(x, y) mean(x(y)), vectors, indPulse);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
