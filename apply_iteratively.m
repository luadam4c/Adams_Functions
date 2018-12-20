function varargout = apply_iteratively (myFunction, array, varargin)
%% Applies a function iteratively to the first argument
% Usage: varargout = apply_iteratively (myFunction, array, varargin)
% Explanation:
%       Applies a function iteratively to an array.
%           The function must be able to take elements of the array as an argument
%               and return outputs that can be retaken as input
% Example(s):
%       a = apply_iteratively(@max, magic(3))
%       b = apply_iteratively(@min, {1:10, -10:5, 5:30})
%       c = apply_iteratively(@max, {1:10, -10:5, 5:30})
% Outputs:
%       varargout     - TODO: Description of varargout
%                   specified as a TODO
% Arguments:
%       myFunction  - a custom function
%                   must be a function handle
%       array       - an array to apply the function iteratively
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_traces.m

% File History:
% 2018-12-19 Created by Adam Lu
% TODO: Does this work for more than one output?
% 
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
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, myFunction, array, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Get the number of output arguments that will be returned by the function
nOutputsPossible = nargout(myFunction);

% Check the number of outputs requested
nargoutchk(0, nOutputsPossible);

% Get the number of outputs requested
nOutputsRequested = nargout;

% Initialize a cell array with the requested number of outputs
varargout = cell(nOutputsRequested, 1);

%% Do the job
if iscell(array)
    [varargout{:}] = myFunction(cellfun(myFunction, array));
else
    [varargout{:}] = myFunction(myFunction(array));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%