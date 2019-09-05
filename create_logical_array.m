function logicalArray = create_logical_array (indices, varargin)
%% Creates a logical array from indices for true and dimensions
% Usage: logicalArray = create_logical_array (indices, dimensions, varargin)
% Explanation:
%       TODO
% Example(s):
%       create_logical_array([2, 4])
%       create_logical_array([2, 4], [5, 1])
% Outputs:
%       logicalArray- logical array created
%                   specified as a logical array
% Arguments:
%       indices     - indices for true
%                   must be a positive integer array
%       dimensions  - (opt) dimensions for the logicial array
%                   must be a positive integer vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/detect_spikes_multiunit.m

% File History:
% 2019-02-20 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
dimensionsDefault = [];

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
addRequired(iP, 'indices', ...
    @(x) validateattributes(x, {'numeric'}, {'3d', 'positive', 'integer'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'dimensions', dimensionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, indices, varargin{:});
dimensions = iP.Results.dimensions;

%% Preparation
% Set default dimensions if not provided
if isempty(dimensions)
    % Find the maximum index for true
    maxIndex = max(indices);

    % Create a column vector
    dimensions = [maxIndex, 1];
end

%% Do the job
% Initialize all elements as false
logicalArray = false(dimensions);

% Change the given indices to true
logicalArray(indices) = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%