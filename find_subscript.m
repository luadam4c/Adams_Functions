function [row, col] = find_subscript (matrix, valueOrFunc, varargin)
%% Find the subscripts (row and column) of the element matching a specific value (or matching a scalar function of the matrix) in a matrix
% Usage: [row, col] = find_subscript (matrix, valueOrFunc, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [row, col] = find_subscript (magic(3), 3)
%       [row, col] = find_subscript (magic(3), @min)
%
% Outputs:
%       row         - row of element found
%                   specified as a positive integer scalar
%       col         - column of element found
%                   specified as a positive integer scalar
%
% Arguments:
%       matrix      - an matrix to apply the function iteratively
%                   must be a matrix
%       valueOrFunc - a value or a custom function
%                   must be a value or a function handle
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for apply_iteratively()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/apply_iteratively.m
%
% Used by:
%       cd/crosscorr_profile.m

% File History:
% 2021-05-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'matrix');
addRequired(iP, 'valueOrFunc');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, matrix, valueOrFunc, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the apply_iteratively() function
otherArguments = iP.Unmatched;

%% Do the job
% Locate the desired value or apply scalar function to matrix
if isa(valueOrFunc, 'function_handle')
    value = apply_iteratively(valueOrFunc, matrix, otherArguments);
else
    value = valueOrFunc;
end

% Find the first index with maximum correlation
idxLinear = find(matrix == value, 1, 'first');

% Find the index
[row, col] = ind2sub(size(matrix), idxLinear);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%