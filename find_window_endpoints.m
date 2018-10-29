function endPoints = find_window_endpoints (timeWindows, timeVecs, varargin)
%% Returns the start and end indices of a time window in a time vector
% Usage: endPoints = find_window_endpoints (timeWindows, timeVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       endPoints1 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'inclusive')
%       endPoints2 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'leftadjust')
%       endPoints3 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'rightadjust')
%       endPoints4 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'restrictive')
%       endPoints5 = find_window_endpoints([1.5, 3.5], {1:5, 0:6})
%       endPoints6 = find_window_endpoints([0.5, 1.5; 2.5, 3.5], 0:6)
%       endPoints7 = find_window_endpoints({[0.5, 1.5], [2.5; 3.5]}, 0:6)
%       endPoints8 = find_window_endpoints({[], [2.5; 3.5]}, 0:6)
%
% Outputs:
%       endPoints   - index(ices) of window endpoints
%                   specified as a positive integer column vector
%                       or a cell array of positive integer column vectors
% Arguments:
%       timeWindows - time window(s); if empty, returns the first and last index
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%       timeVecs    - time vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'BoundaryMode': boundary mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'inclusive'   - the time endPoints approximates
%                                           must include the time window
%                       'leftadjust'  - the time endPoints 
%                                           the time window to the left
%                       'rightadjust' - the time endPoints approximates
%                                           the time window to the right
%                       'restrictive' - the time endPoints 
%                                           must be within the time window
%                   default == 'restrictive'
%
% Requires:
%       cd/iscellnumeric.m
%       cd/iscellnumericvector.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_subvectors.m
%       cd/find_passive_params.m

% File History:
% 2018-10-09 Created by Adam Lu
% 2018-10-25 Now accepts a cell array of time vectors
% 2018-10-27 Now returns first and last index if the window is empty
% 2018-10-27 Now returns endPoints as a single output
%               and is a cell array if multiple vectors are passed in
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
addRequired(iP, 'timeWindows', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['timeWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['timeVecs must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BoundaryMode', boundaryModeDefault, ...
    @(x) any(validatestring(x, validBoundaryModes)));

% Read from the Input Parser
parse(iP, timeWindows, timeVecs, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);

%% Preparation
% Match the formats of timeWindows and timeVecs so that cellfun can be used
[timeWindows, timeVecs] = ...
    match_format_vectors(timeWindows, timeVecs, 'ForceCellOutputs', false);

% If the time window is a nonempty numeric array, make sure it has two rows
if ~isempty(timeWindows) && isnumeric(timeWindows) && ...
    size(timeWindows, 1) ~= 2
    error('Time windows must have only two elements!');
end

%% Do the job
if iscell(timeVecs)
    endPoints = ...
        cellfun(@(x, y) find_window_endpoints_helper(x, y, boundaryMode), ...
                timeWindows, timeVecs, 'UniformOutput', false);
else
    endPoints = ...
        find_window_endpoints_helper(timeWindows, timeVecs, boundaryMode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_window_endpoints_helper (timeWindow, timeVec, ...
                                                    boundaryMode)
%% Find the endPoints for one time vector

% If the time window is empty, return the first and last indices
if isempty(timeWindow)
    endPoints = [1; length(timeVec)];
    return
end

% Get the time to start
timeStart = timeWindow(1);

% Get the time to end
timeEnd = timeWindow(2);

% Decide on endPoints based on boundary mode:
%   'inclusive'   - the time endPoints must include the time window
%   'leftadjust'  - the time endPoints approximates the time window to the left
%   'rightadjust' - the time endPoints approximates the time window to the right
%   'restrictive' - the time endPoints must be within the time window
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
    endPoints = [];
else
    endPoints = [idxStart; idxEnd];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%                   must be a nondecreasing numeric vector with 2 elements
@(x) validateattributes(x, {'numeric'}, ...
        {'vector', 'nondecreasing', 'numel', 2}));

function [idxStarts, idxEnds] = find_window_endpoints (timeWindows, timeVecs, varargin)

idxStart = [];
idxEnd = [];

%       idxStarts   - index(ices) of window start
%                   specified as a positive integer vector
%       idxEnds     - index(ices) of window end
%                   specified as a positive integer vector

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%