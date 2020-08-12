function [endPoints, windows, timeInstants] = ...
                compute_running_windows (timeVec, windowSize, varargin)
%% Computes running windows based on time vectors
% Usage: [endPoints, windows, timeInstants] = ...
%               compute_running_windows (timeVec, windowSize, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [endPoints, windows, timeInstants] = compute_running_windows(0:1:10, 2)
%       [endPoints, windows, timeInstants] = compute_running_windows(0:1:10, [2; 4])
%
% Outputs:
%       endPoints   - time window endpoints in time vector
%                   specified as a TODO
%       windows     - time windows
%                   specified as a TODO
%       timeInstants    - time instants (right end points) for windows
%                       specified as a TODO 
%
% Arguments:
%       timeVec     - time vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       windowSize  - window size(s) for each time vector
%                   must be a positive numeric vectors
%       varargin    - 'Resolution': resolution for time instants
%                   must be a positive numeric vector
%                   default == windowSize / 2
%                   - Any other parameter-value pair for find_window_endpoints()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_row_vector.m
%       cd/iscellnumeric.m
%       cd/ispositivevector.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/parse_pleth_trace.m

% File History:
% 2020-08-12 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
resolutionDefault = [];         % set later

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
addRequired(iP, 'timeVec', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'windowSize', ...
    @(x) isempty(x) || ispositivevector(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Resolution', resolutionDefault, ...
    @(x) isempty(x) || ispositivevector(x));

% Read from the Input Parser
parse(iP, timeVec, windowSize, varargin{:});
resolution = iP.Results.Resolution;

% Keep unmatched arguments for the find_window_endpoints() function
otherArguments = iP.Unmatched;

%% Preparation
% Decide on resolution
if isempty(resolution)
    resolution = windowSize / 2;
end

% Match formats
if iscell(timeVec) || count_vectors(timeVec) > 1 || ...
        numel(windowSize) > 1 || numel(resolution) > 1
    % Force as a column cell array
    timeVec = force_column_cell(timeVec);

    % Match format of vectors
    [timeVec, windowSize, resolution] = ...
        match_format_vectors(timeVec, windowSize, resolution);
end
    
%% Do the job
if iscell(timeVec)
    [endPoints, windows, timeInstants] = ...
        cellfun(@(a, b, c) compute_running_windows_helper(a, b, c, ...
                                    otherArguments), ...
                timeVec, num2cell(windowSize), num2cell(resolution), ...
                'UniformOutput', false);
else
    [endPoints, windows, timeInstants] = ...
        compute_running_windows_helper(timeVec, windowSize, resolution, ...
                                        otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [endPoints, windows, timeInstants] = ...
                compute_running_windows_helper (timeVec, windowSize, ...
                                                resolution, otherArguments)

% Compute the number of time windows
nWindows = floor(timeVec(end) / resolution);

% Compute the time instants
timeInstants = resolution * (1:nWindows)';

% Force as a row vector
timeInstantsRow = force_row_vector(timeInstants);

% Compute running windows
windows = [timeInstantsRow - windowSize; timeInstantsRow];

% Compute endpoints for running windows
endPoints = find_window_endpoints(windows, timeVec, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%