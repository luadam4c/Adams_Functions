function maxNumel = compute_maximum_numel (list, varargin)
%% Given a list of arrays, compute the maximum number of elements
% Usage: maxNumel = compute_maximum_numel (list, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       compute_maximum_numel(blab)
%       compute_maximum_numel(myStruct)
%       compute_maximum_numel(myCellNumeric3D)
%
% Outputs:
%       maxNumel    - maximum number of elements
%                   specified as a positive integer scalar
%
% Arguments:
%       list        - a list of arrays
%                   must be a cell array or a structure
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-01-05 Created by Adam Lu
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'list', ...
    @(x) validateattributes(x, {'cell', 'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, list, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Count the number of values in each array
if iscell(list)
    allNumels = cellfun(@numel, list);
elseif isstruct(list)
    allNumels = structfun(@numel, list);
end

% Compute the maximum number of values
maxNumel = max(allNumels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%