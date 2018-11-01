function residuals = compute_residuals (simVectors, realVectors, varargin)
%% Computes residual vector(s) from simulated and recorded vectors
% Usage: residuals = compute_residuals (simVectors, realVectors, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       residuals   - residual vectors
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
% Arguments:    
%       simVectors  - simulated vectors
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       realVectors - recorded vectors
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vectors
%                   must be a numeric array or a cell array of numeric arrays
%
% Requires:
%       cd/iscellnumeric.m
%       cd/match_format_vectors.m
%
% Used by:    
%       cd/m3ha_run_neuron_once.m

% File History:
% 2018-10-23 Created by Adam Lu
% 2018-10-28 Now uses match_format_vectors.m and accepts numeric arrays
%               with multiple columns
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'simVectors', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['simVectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'realVectors', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['realVectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, simVectors, realVectors, varargin{:});

%% Preparation
% Make sure simVectors and realVectors are both cell arrays of the same length
%   if one of them has multiple vectors
[simVectors, realVectors] = ...
    match_format_vectors(simVectors, realVectors, 'ForceCellOutputs', false);

%% Compute residuals
if iscell(simVectors) && iscell(realVectors)
    % Return residual vectors
    residuals = cellfun(@(x, y) x - y, simVectors, realVectors, ...
                        'UniformOutput', false);
elseif isnumeric(simVectors) && isnumeric(realVectors)
    % Return the residual vector
    residuals = simVectors - realVectors;
else
    error('Error in code logic!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

residuals = cell(1, nSweeps);         % residuals (sim - real)
parfor iSwp = 1:nSweeps
    % Get the voltage vector from the simulated data
    simVectors = simVectors{iSwp}(:, 2);

    % Get the voltage vector from the real data
    realVectors = realData{iSwp}(:, 2);

    % Compute the residual
    residuals{iSwp} = simVectors - realVectors;
end

%       cd/match_array_counts.m
[simVectors, realVectors] = ...
    match_array_counts(simVectors, realVectors, 'ForceCellOutputs', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%