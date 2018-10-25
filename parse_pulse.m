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
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/count_samples.m
%       cd/find_pulse_endpoints.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%
% Used by:    
%       cd/find_passive_params.m
%       cd/plot_pulse.m

% File History:
% 2018-10-10 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

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

% Read from the Input Parser
parse(iP, vectors, varargin{:});

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
indBase = arrayfun(@(x, y, z) [1:x, y:z], ...
                    idxBeforeStart, idxAfterEnd, nSamples, ...
                    'UniformOutput', false);

% Find the average baseline value
baseValue = cellfun(@(x, y) mean(x(y)), vectors, indBase);

% Find the pulse indices
indPulse = arrayfun(@(x, y) [x:y], idxAfterStart, idxBeforeEnd, ...
                    'UniformOutput', false);

% Find the average pulse value
pulseValue = cellfun(@(x, y) mean(x(y)), vectors, indPulse);

% Find the pulse amplitudes
pulseAmplitude = pulseValue - baseValue;

%% Store results in output
parsedParams = table(nSamples, pulseWidthSamples, ...
                    idxBeforeStart, idxBeforeEnd, ...
                    idxAfterStart, idxAfterEnd, idxMidpoint, ...
                    baseValue, pulseValue, pulseAmplitude);
parsedData = table(vectors, indBase, indPulse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%