function varargout = ...
                compute_running_windows (timeVec, windowSize, varargin)
%% Computes running windows based on time vectors
% Usage: [endPoints, timeInstants, windows] = ...
%               compute_running_windows (timeVec, windowSize, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [endPoints, timeInstants, windows] = compute_running_windows(1:10, 2)
%       [endPoints, timeInstants, windows] = compute_running_windows(1:10, 2, 'IsRegular', true)
%       [endPoints, timeInstants, windows] = compute_running_windows(1:10, [2; 4])
%
% Outputs:
%       endPoints   - time window endpoints in time vector
%                   specified as a TODO
%       timeInstants    - time instants (right end points) for windows
%                       specified as a TODO 
%       windows     - time windows
%                   specified as a TODO
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
%                   - 'IsRegular': whether time vector is regular
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for find_window_endpoints()
%
% Requires:
%       cd/compute_sampling_interval.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%       cd/ispositivevector.m
%       cd/match_format_vectors.m
%       cd/vecfun.m
%
% Used by:
%       cd/parse_pleth_trace.m

% File History:
% 2020-08-12 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
resolutionDefault = [];         % set later
isRegularDefault = false;       % don't assume regularity by default

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
addParameter(iP, 'IsRegular', isRegularDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, timeVec, windowSize, varargin{:});
resolution = iP.Results.Resolution;
isRegular = iP.Results.IsRegular;

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
    [varargout{1:nargout}] = ...
        cellfun(@(a, b, c) compute_running_windows_helper(a, b, c, ...
                                    isRegular, otherArguments), ...
                timeVec, num2cell(windowSize), num2cell(resolution), ...
                'UniformOutput', false);
else
    [varargout{1:nargout}] = ...
        compute_running_windows_helper(timeVec, windowSize, resolution, ...
                                        isRegular, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = compute_running_windows_helper (timeVec, windowSize, ...
                                        resolution, isRegular, otherArguments)

% Two methods
if isRegular
    % Compute the sampling interval
    samplingInterval = compute_sampling_interval(timeVec);

    % Compute the time window in samples
    timeWindowSamples = ceil(windowSize / samplingInterval);

    % Compute the resolution in samples
    resolutionSamples = ceil(resolution / samplingInterval);

    % Compute the number of time windows
    nWindows = floor(numel(timeVec) / resolutionSamples);

    % Compute the time instant indices
    indTimeInstants = resolutionSamples * (1:nWindows);

    % Compute running window endpoints
    endPoints = [indTimeInstants - timeWindowSamples; indTimeInstants];

    % Replace nonpositive values with 1
    endPoints(endPoints <= 0) = 1;

    % Extract time instants
    if nargout >= 2
        timeInstants = timeVec(indTimeInstants);
    end

    % Compute running windows
    % TODO: Make more efficient
    if nargout >= 3
        windows = vecfun(@(x) extract_subvectors(timeVec, 'Indices', x), ...
                        endPoints);
    end
else
    % Compute the number of time windows
    nWindows = floor(timeVec(end) / resolution);

    % Compute the time instants
    timeInstants = resolution * (1:nWindows);

    % Compute running windows
    windows = [timeInstants - windowSize; timeInstants];

    % Compute endpoints for running windows
    endPoints = find_window_endpoints(windows, timeVec, otherArguments);

    % Match vector formats
    if iscolumn(timeVec) && nargout >= 2
        timeInstants = force_column_vector(timeInstants);
    end
end

%% Outputs
varargout{1} = endPoints;
if nargout >= 2
    varargout{2} = timeInstants;
end
if nargout >= 3
    varargout{3} = windows;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%