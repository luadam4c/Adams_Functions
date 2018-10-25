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
%                   must be a numeric vector or a cell array of numeric vectors
%       realVectors - recorded vectors
%                   must be a numeric vector or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumeric.m
%       cd/match_vector_counts.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-23 Created by Adam Lu
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
%   if one of them is a cell array
[simVectors, realVectors] = ...
    match_vector_counts(simVectors, realVectors, 'ForceCellOutputs', false);

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%