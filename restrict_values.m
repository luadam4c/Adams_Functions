function values = restrict_values (values, varargin)
%% Restrict values
% Usage: values = restrict_values (values, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       restrict_values(magic(3), 'LowerBound', 5)
%       restrict_values({6:9, 1:3}, 'UpperBound', 4)
%
% Outputs:
%       values      - restricted values
%                   specified as an array
%
% Arguments:
%       values      - original values
%                   must be an array
%       varargin    - 'LowerBound': lower bound of values
%                   must be a numeric scalar
%                   default == []
%                   - 'UpperBound': upper bound of values
%                   must be a numeric scalar
%                   default == []
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_find_decision_point.m

% File History:
% 2020-05-13 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
lowerBoundDefault = [];
upperBoundDefault = [];

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
addRequired(iP, 'values');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LowerBound', lowerBoundDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'UpperBound', upperBoundDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, values, varargin{:});
lowerBound = iP.Results.LowerBound;
upperBound = iP.Results.UpperBound;

%% Do the job
if iscell(values)
    values = cellfun(@(x) restrict_values_helper(x, lowerBound, upperBound), ...
                    values, 'UniformOutput', false);
else
    values = restrict_values_helper(values, lowerBound, upperBound);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function values = restrict_values_helper (values, lowerBound, upperBound)
%% Restrict values in a non-cell array

% Replace values with lower bound
if ~isempty(lowerBound)
    values(values < lowerBound) = lowerBound;
end

% Replace values with upper bound
if ~isempty(upperBound)
    values(values > upperBound) = upperBound;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%