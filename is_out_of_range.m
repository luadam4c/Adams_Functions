function [isOutOfRangeAny, isOutOfRangeAll, isOutOfRangeEachElement] = ...
                is_out_of_range (values, rangeToCheck, varargin)
%% Check if any of the value(s) are out of range
% Usage: [isOutOfRangeAny, isOutOfRangeAll, isOutOfRangeEachElement] = ...
%               is_out_of_range (values, rangeToCheck, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isOutOfRangeAny - whether any value is out of range
%                       specified as a logical scalar
%       isOutOfRangeAll - whether all values are out of range
%                       specified as a logical scalar
%       isOutOfRangeEachElement - whether each value is out of range
%                               specified as a logical array
% Arguments:
%       values      - TODO: Description of values
%                   must be a TODO
%       rangeToCheck    - TODO: Description of values
%                       must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/annotation_in_plot.m

% File History:
% 2018-12-28 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'values'); %, ...
    % TODO: validation function %);
addRequired(iP, 'rangeToCheck'); %, ...
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, values, rangeToCheck, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Check if each value is out of range
isOutOfRangeEachElement = values < rangeToCheck(1) | values > rangeToCheck(2);

% Check if all of the values are out of range
isOutOfRangeAll = apply_iteratively(@all, isOutOfRangeEachElement);

% Check if any of the values are out of range
isOutOfRangeAny = apply_iteratively(@any, isOutOfRangeEachElement);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%