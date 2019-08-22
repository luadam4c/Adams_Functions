function numArray = cell2num (cellArray, varargin)
%% This is the reverse of num2cell, replacing empty entries with NaNs
% Usage: numArray = cell2num (cellArray, varargin)
% Explanation:
%       TODO
% Example(s):
%       cell2num({[], 2, [3, 4, 5], []})
%       cell2num({[], 2; [3, 4, 5], []})
% Outputs:
%       numArray    - numeric array
%                   specified as a numeric array
% Arguments:
%       cellArray   - cell array with at most one number per cell
%                   must be a cell array of numeric vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%
% Used by:
%       cd/compute_sampsizepwr.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-08-20 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'cellArray');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, cellArray, varargin{:});
param1 = iP.Results.param1;

%% Do the job
numArray = cellfun(@force_numeric_scalar, cellArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = force_numeric_scalar(value)

if isempty(value) || ~isnum(value)
    out = NaN;
elseif numel(value) > 1
    out = value(1);
else
    out = value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
