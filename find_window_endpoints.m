function [idxStart, idxEnd] = find_window_endpoints (timeWindow, timeVec, varargin)
%% Returns the start and end indices of a time window in a time vector
% Usage: [idxStart, idxEnd] = find_window_endpoints (timeWindow, timeVec, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       idxStart    - index of window start
%                   specified as a positive integer scalar
%       idxEnd      - index of window end
%                   specified as a positive integer scalar
% Arguments:
%       timeWindow  - time window
%                   must be a numeric vector or 2 elements
%       timeVec     - time vector
%                   must be a numeric vector
%
% Used by:    
%       cd/find_passive_params.m

% File History:
% 2018-10-09 Created by Adam Lu
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
addRequired(iP, 'timeWindow', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addRequired(iP, 'timeVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, timeWindow, timeVec, varargin{:});

%% Do the job
idxStart = find(timeVec >= timeWindow(1), 1);
idxEnd = find(timeVec <= timeWindow(2), 1, 'last');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%