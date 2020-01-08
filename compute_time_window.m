function windows = compute_time_window (timeVecs, varargin)
%% Computes time windows from time vectors and given time end points
% Usage: windows = compute_time_window (timeVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_time_window({2:100, 50:200})
%       compute_time_window({2:100, 50:200}, 'TimeEnd', 150)
%       compute_time_window({2:100, 50:200}, 'TimeStart', 3)
%
% Outputs:
%       windows     - computed time windows
%                   specified as a numeric array with 2 rows
% Arguments:
%       timeVecs    - time vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'TimeStart': time of window start
%                   must be a numeric vector
%                   default == start times of time vectors
%                   - 'TimeEnd': time of window end
%                   must be a numeric vector
%                   default == end times of time vectors
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/compute_oscillation_duration.m
%       cd/detect_spikes_multiunit.m
%       cd/parse_multiunit.m

% File History:
% 2019-05-14 Created by Adam Lu
% TODO: Fix time window to start from tBase if timeStart is 1?

%% Hard-coded parameters

%% Default values for optional arguments
timeStartDefault = [];    % set later
timeEndDefault = [];      % set later

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
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['timeVecs must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, timeVecs, varargin{:});
timeStart = iP.Results.TimeStart;
timeEnd = iP.Results.TimeEnd;

%% Preparation
% Get the starting time(s)
if isempty(timeStart)
    timeStart = extract_elements(timeVecs, 'first');
end

if isempty(timeEnd)
    timeEnd = extract_elements(timeVecs, 'last');
end

% Match the vector formats as column vectors
[timeStart, timeEnd] = match_format_vectors(timeStart, timeEnd);

%% Do the job
% Use timeStart to stimStartMs by default
windows = transpose([timeStart, timeEnd]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%