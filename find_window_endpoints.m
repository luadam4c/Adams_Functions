function endPoints = find_window_endpoints (timeWindows, timeVecs, varargin)
%% Returns the start and end indices of a time window in a time vector
% Usage: endPoints = find_window_endpoints (timeWindows, timeVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       endPoints0 = find_window_endpoints([1.5, 3.5], 1:5)
%       endPoints0 = find_window_endpoints([2.5, 2.5], 1:5)
%       endPoints1 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'inclusive')
%       endPoints2 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'leftadjust')
%       endPoints3 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'rightadjust')
%       endPoints4 = find_window_endpoints([1.5, 3.5], 1:5, 'BoundaryMode', 'restrictive')
%       endPoints5 = find_window_endpoints([1.5, 3.5], {1:5, 0:6})
%       endPoints6 = find_window_endpoints([1.5, 3.5], magic(3))
%       endPoints7 = find_window_endpoints([0.5, 1.5; 2.5, 3.5], 0:6)
%       endPoints8 = find_window_endpoints({[0.5, 1.5], [2.5; 3.5]}, 0:6)
%       endPoints9 = find_window_endpoints({[0.5, 1.5], [2.5; 3.5]}, 0:2:4)
%       endPoints10 = find_window_endpoints({[], [2.5; 3.5]}, 0:6)
%       endPoints11 = find_window_endpoints([0.5, 1.5; 2.5, 3.5], [(0:6)', (1:7)'])
%       endPoints12 = find_window_endpoints([3.5, 1.5], 5:-1:1, 'BoundaryMode', 'inclusive')
%       endPoints8 = find_window_endpoints({[0.5, 1.5], [2.5; 3.5]})
%
% Outputs:
%       endPoints   - index(ices) of window endpoints
%                   specified as a 2-row numeric matrix
%                       or a cell array of 2-row numeric matrices
%
% Arguments:
%       timeWindows - time window(s); if empty, returns the first and last index
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%       timeVecs    - time vector(s), must be monotonic
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       varargin    - 'BoundaryMode': boundary mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'inclusive'   - the time endPoints 
%                                           must include the time window
%                       'leftadjust'  - the time endPoints approximates
%                                           the time window to the left
%                       'rightadjust' - the time endPoints approximates
%                                           the time window to the right
%                       'restrictive' - the time endPoints 
%                                           must be within the time window
%                   default == 'restrictive'
%                   - 'ForceMatrixOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'WarnFlag': whether to warn if no files found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/count_samples.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/iscellnumeric.m
%       cd/iscellnumericvector.m
%       cd/isnum.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_baseline_noise.m
%       cd/compute_lts_errors.m
%       cd/compute_running_windows.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/compute_time_average.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_passive_params.m
%       cd/read_neuron_outputs.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_ipsc.m
%       cd/plot_fitted_traces.m
%       cd/plot_spectrogram_multiunit.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m

% File History:
% 2018-10-09 Created by Adam Lu
% 2018-10-25 Now accepts a cell array of time vectors
% 2018-10-27 Now returns first and last index if the window is empty
% 2018-10-27 Now returns endPoints as a single output
%               and is a cell array if multiple vectors are passed in
% 2019-01-03 Now accepts a cell array of non-vector arrays and 
%               added 'ForceMatrixOutput' as an optional argument
% 2019-09-10 Added 'WarnFlag' as an optional flag
% 2019-11-14 Now allows time vectors to be decreasing
% 2020-01-02 Now returns matrices if inputs are all matrices
% 2020-04-21 Now returns empty if timeVec is empty
% 2020-08-12 Improved performance when there is one vectors and many windows

%% Hard-coded parameters
validBoundaryModes = {'inclusive', 'leftadjust', 'rightadjust', 'restrictive'};

%% Default values for optional arguments
boundaryModeDefault = 'restrictive';
forceMatrixOutputDefault = logical.empty; % set later
warnFlagDefault = true;         % warn if windows out of bounds by default

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
addRequired(iP, 'timeWindows', ...
    @(x) assert(isnum(x) || iscellnumeric(x), ...
                ['timeWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnum(x) || iscellnumeric(x), ...
                ['timeVecs must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BoundaryMode', boundaryModeDefault, ...
    @(x) any(validatestring(x, validBoundaryModes)));
addParameter(iP, 'ForceMatrixOutput', forceMatrixOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'WarnFlag', warnFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, timeWindows, timeVecs, varargin{:});
boundaryMode = validatestring(iP.Results.BoundaryMode, validBoundaryModes);
forceMatrixOutputUser = iP.Results.ForceMatrixOutput;
warnFlag = iP.Results.WarnFlag;

%% Preparation
% Decide whether to force output as a matrix if not provided
if isempty(forceMatrixOutputUser)
    % Force output as a matrix only if the original formats
    %   of both windows and vectors are matrices
    forceMatrixOutput = ~iscell(timeWindows) && ismatrix(timeWindows) && ...
                        ~iscell(timeVecs) && ismatrix(timeVecs);
else
    forceMatrixOutput = forceMatrixOutputUser;
end

% Match the formats of timeWindows and timeVecs so that cellfun can be used
if iscell(timeVecs)
    [timeWindows, timeVecs] = ...
        match_format_vector_sets(timeWindows, timeVecs, ...
                                'ForceCellOutputs', false);
end

% Force as column vectors
[timeWindows, timeVecs] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVector', true), ...
            timeWindows, timeVecs);

% If the time window is a nonempty numeric array, make sure it has two rows
if ~isempty(timeWindows) && isnumeric(timeWindows) && ...
    size(timeWindows, 1) ~= 2
    error('Time windows must have only two rows!');
end

%% Do the job
if iscellnumericvector(timeVecs)
    % Find end points for each vector
    endPoints = ...
        array_fun(@(x, y) find_window_endpoints_helper(x, y, ...
                                                boundaryMode, warnFlag), ...
                timeWindows, timeVecs, 'UniformOutput', false);

    % Force as a matrix if requested
    if forceMatrixOutput
        endPoints = force_matrix(endPoints);
    end
elseif iscell(timeVecs)
    % Find end points for each array
    endPoints = ...
        array_fun(@(x, y) find_window_endpoints(x, y, ...
                        'BoundaryMode', boundaryMode, ...
                        'ForceMatrixOutput', forceMatrixOutputUser), ...
                timeWindows, timeVecs, 'UniformOutput', false);
else
    endPoints = find_window_endpoints_helper(timeWindows, timeVecs, ...
                                                boundaryMode, warnFlag);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_window_endpoints_helper (timeWindow, timeVec, ...
                                                    boundaryMode, warnFlag)
%% Find the endPoints for time vector(s)

% Return empty if timeVec is empty
if isempty(timeVec)
    endPoints = [];
    return
end

% Compute the number of time points
nTimePoints = size(timeVec, 1);

% Compute the number of vectors
nVectors = size(timeVec, 2);

% Compute the number of windows
nWindows = size(timeWindow, 2);

% If the time window is empty, return the first and last indices
if isempty(timeWindow)
    endPoints = [1; nTimePoints];
    return
end

% Make sure time window is increasing
timeWindow = sort(timeWindow, 'ascend');

% Force time vector to be increasing
if timeVec(end, 1) < timeVec(1, 1)
    timeVec = flipud(timeVec);
    vectorFlipped = true;
else
    vectorFlipped = false;
end

% Get the time to start
timeStart = timeWindow(1, :);

% Get the time to end
timeEnd = timeWindow(2, :);

% Decide on first and last indices for each vector based on boundaryMode
%   Note: this should return row vectors
if nVectors > 1
    [idxStart, idxEnd] = ...
        array_fun(@(x) decide_on_endpoints(timeVec(:, x), ...
                        timeStart(x), timeEnd(x), boundaryMode, warnFlag), ...
                    1:nVectors, 'UniformOutput', true);
else
    [idxStart, idxEnd] = ...
        array_fun(@(a, b) decide_on_endpoints(timeVec, a, b, ...
                                            boundaryMode, warnFlag), ...
                    timeStart, timeEnd, 'UniformOutput', true);
end

% Create endpoints
if vectorFlipped
    endPoints = nTimePoints - [idxEnd; idxStart] + 1;
else
    endPoints = [idxStart; idxEnd];
end

% If all end points are NaN, return empty
if all(all(all(isnan(endPoints))))
    endPoints = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxStart, idxEnd] = ...
                decide_on_endpoints (timeVec, timeStart, timeEnd, ...
                                    boundaryMode, warnFlag)
%% Decide on endpoints for a single vector

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
            if warnFlag
                fprintf('Time window is in between sample points!\n\n');
            end
            idxStart = [];
            idxEnd = [];
        end
    otherwise
        error('boundaryMode unrecognized!');
end

% If either is empty, the time window is out of range, 
%   so print message and return NaNs
if isempty(idxStart) || isempty(idxEnd)
    if warnFlag
        fprintf('Time window is out of range!\n\n');
    end
    idxStart = NaN;
    idxEnd = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%                   must be a nondecreasing numeric vector with 2 elements
@(x) validateattributes(x, {'numeric'}, ...
        {'vector', 'nondecreasing', 'numel', 2}));

function [idxStarts, idxEnds] = ...
    find_window_endpoints (timeWindows, timeVecs, varargin)

idxStart = [];
idxEnd = [];

%       idxStarts   - index(ices) of window start
%                   specified as a positive integer vector
%       idxEnds     - index(ices) of window end
%                   specified as a positive integer vector

addRequired(iP, 'timeVecs', ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['timeVecs must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Collapse columns if they are the same
[timeWindow, timeVec] = argfun(@collapse_identical_vectors, timeWindow, timeVec);

% Count the number of samples
nSamples = count_samples(timeVec);
endPoints = transpose([ones(size(nSamples)), nSamples]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
