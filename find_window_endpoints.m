function [idxStarts, idxEnds] = find_window_endpoints (timeWindow, timeVecs, varargin)
%% Returns the start and end indices of a time window in a time vector
% Usage: [idxStarts, idxEnds] = find_window_endpoints (timeWindow, timeVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       idxStarts   - index(ices) of window start
%                   specified as a positive integer vector
%       idxEnds     - index(ices) of window end
%                   specified as a positive integer vector
% Arguments:
%       timeWindow  - time window
%                   must be a numeric vector or 2 elements
%       timeVecs    - time vector(s)
%                   must be a numeric vector
%                       or a cell array of numeric vectors
%
% Requires:
%       cd/iscellnumericvector.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/find_passive_params.m

% File History:
% 2018-10-09 Created by Adam Lu
% 2018-10-25 Now accepts a cell array of time vectors
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
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) && isvector(x) || iscellnumericvector(x), ...
                ['timeVecs must be either a numeric vector', ...
                    'or a cell array of numeric vectors!']));

% Read from the Input Parser
parse(iP, timeWindow, timeVecs, varargin{:});

%% Do the job
if iscell(timeVecs)
    [idxStarts, idxEnds] = ...
        cellfun(@(x, y) find_window_endpoints_helper(x, y), ...
                timeWindow, timeVecs);
else
    [idxStarts, idxEnds] = find_window_endpoints_helper(timeWindow, timeVecs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxStart, idxEnd] = find_window_endpoints_helper (timeWindow, timeVec)
%% Find the endpoints for one time vector

idxStart = find(timeVec >= timeWindow(1), 1);
idxEnd = find(timeVec <= timeWindow(2), 1, 'last');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%