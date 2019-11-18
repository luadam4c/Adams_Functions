function vecOut = apply_to_nonnan_part (myFunction, vecIn, varargin)
%% Applies a function to just the non-NaN part of a vector
% Usage: vecOut = apply_to_nonnan_part (myFunction, vecIn, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       smooth([nan(10, 1); transpose(1:10)])
%       apply_to_nonnan_part(@smooth, [nan(10, 1); transpose(1:10)])
%
% Outputs:
%       vecOut      - output vector
%                   specified as a vector
%
% Arguments:
%       myFunction  - a custom function that takes a vector 
%                           as the first argument
%                       e.g., smooth(), medfilt1()
%                   must be a function handle
%       vecIn       - input vector
%                   must be a vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/movingaveragefilter.m

% File History:
% 2019-11-18 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'vecIn');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, myFunction, vecIn, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Initialize output with input
vecOut = vecIn;

% Determine whether each value is NaN
isNan = isnan(vecIn);

% If all values are NaNs, return NaNs
if all(isNan)
    return
end

% Get the starting index of non-NaN values
idxStart = find(~isNan, 1, 'first');

% Get the ending index of non-NaN values
idxEnd = find(~isNan, 1, 'last');

% Apply the function to just the indices of interest
vecOut(idxStart:idxEnd) = myFunction(vecIn(idxStart:idxEnd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%