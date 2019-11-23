function newValues = choose_random_values(lowerBounds, upperBounds, varargin)
%% Chooses random values from bounds
% Usage: newValues = choose_random_values(lowerBounds, upperBounds, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       newValues   - new values
%                   specified as a numeric vector
% Arguments:
%       lowerBounds - lower bounds for new values
%                   must be a numeric vector
%       upperBounds - upper bounds for new values
%                   must be a numeric vector
%       varargin    - 'IsLog': whether the parameter should be log-scaled
%                   must be a binary vector
%                   default == false(nValues, 1)
%
% Used by:
%       cd/m3ha_neuron_create_new_initial_params

% File History:
% 2018-12-11 Created by Adam Lu
% TODO: Use force_column_vector.m ('IgnoreArrays', true) and match_row_count.m 
%           for the first two arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
isLogDefault = [];              % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'lowerBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'upperBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IsLog', isLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));

% Read from the Input Parser
parse(iP, lowerBounds, upperBounds, varargin{:});
isLogs = iP.Results.IsLog;

% Check lengths of lower and upper bounds vectors
if length(lowerBounds) ~= length(upperBounds)
    error("There must be an equal number of lower and upper bounds!\n");
end

%% Preparation
% Decide on default isLogs
if isempty(isLogs)
    % None of them log-scaled by default
    isLogs = false(size(lowerBounds));
end

%% Do the job
newValues = arrayfun(@(x, y, z) choose_random_values_helper(x, y, z), ...
                        lowerBounds, upperBounds, isLogs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newValue = choose_random_values_helper(LB, UB, isLog)
%% Choose a random value from lower and upper bounds and whether the value should be log-scaled

if isLog
    % Linear interpolation in log space
    newValue = exp(log(LB) + rand() * (log(UB) - log(LB)));
else
    % Linear interpolation
    newValue = LB + rand() * (UB - LB);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of values
nValues = length(lowerBounds);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
