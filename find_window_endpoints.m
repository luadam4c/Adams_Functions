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
%                   must be a nondecreasing numeric vector with 2 elements
%       timeVecs    - time vector(s)
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'BoundaryMode': boundary mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'inclusive'   - the time endpoints approximates
%                                           must include the time window
%                       'leftadjust'  - the time endpoints 
%                                           the time window to the left
%                       'rightadjust' - the time endpoints approximates
%                                           the time window to the right
%                       'restrictive' - the time endpoints 
%                                           must be within the time window
%                   default == 'restrictive'
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
validBoundaryModes = {'inclusive', 'leftadjust', 'rightadjust', 'restrictive'};

%% Default values for optional arguments
boundaryModeDefault = 'restrictive';

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
    @(x) validateattributes(x, {'numeric'}, ...
            {'vector', 'nondecreasing', 'numel', 2}));
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) && isvector(x) || iscellnumericvector(x), ...
                ['timeVecs must be either a numeric vector', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BoundaryMode', boundaryModeDefault, ...
    @(x) any(validatestring(x, validBoundaryModes)));

% Read from the Input Parser
parse(iP, timeWindow, timeVecs, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);

%% Do the job
if iscell(timeVecs)
    [idxStarts, idxEnds] = ...
        cellfun(@(x, y) find_window_endpoints_helper(x, y, boundaryMode), ...
                timeWindow, timeVecs);
else
    [idxStarts, idxEnds] = ...
        find_window_endpoints_helper(timeWindow, timeVecs, boundaryMode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxStart, idxEnd] = ...
                find_window_endpoints_helper (timeWindow, timeVec, boundaryMode)
%% Find the endpoints for one time vector

% Get the time to start
timeStart = timeWindow(1);

% Get the time to end
timeEnd = timeWindow(2);

% Decide on endpoints based on boundary mode:
%   'inclusive'   - the time endpoints must include the time window
%   'leftadjust'  - the time endpoints approximates the time window to the left
%   'rightadjust' - the time endpoints approximates the time window to the right
%   'restrictive' - the time endpoints must be within the time window
switch boundaryMode
    case 'inclusive'
        % Find the last time pointbefore the start of the time window
        idxStart = find(timeVec <= timeStart, 1, 'last');

        % Find the first time point after the end of the time window
        idxEnd = find(timeVec >= timeEnd, 1, 'first');
    case 'leftadjust'
        % Find the last time point before the start of the time window
        idxStart = find(timeVec <= timeStart, 1, 'last');

        % Find the last time point before the end of the time window
        idxEnd = find(timeVec <= timeEnd, 1, 'last');
    case 'rightadjust'
        % Find the first time point after the start of the time window
        idxStart = find(timeVec >= timeStart, 1, 'first');

        % Find the first time point after the end of the time window
        idxEnd = find(timeVec >= timeEnd, 1, 'first');
    case 'restrictive'
        % Find the first time point after the start of the time window
        idxStart = find(timeVec >= timeStart, 1, 'first');

        % Find the last time point before the end of the time window
        idxEnd = find(timeVec <= timeEnd, 1, 'last');

        % If idxStart is greater than idxEnd, 
        %   the time window is between sample points, 
        %   so print message and return empty arrays
        if idxStart > idxEnd
            fprintf('Time window is in between sample points!\n\n');
            idxStart = [];
            idxEnd = [];
        end
    otherwise
        error('boundaryMode unrecognized!');
end

% If either is empty, the time window is out of range, 
%   so print message and return empty arrays
if isempty(idxStart) || isempty(idxEnd)
    fprintf('Time window is out of range!\n\n');
    idxStart = [];
    idxEnd = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%