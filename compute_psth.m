function [counts, edges, relEventTimes] = ...
                compute_psth (eventTimes, stimTimes, varargin)
%% Computes a peri-stimulus time histogram
% Usage: [counts, edges, relEventTimes] = ...
%               compute_psth (eventTimes, stimTimes, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [counts, edges] = compute_psth(randi(100, 100, 1), 10:10:80)
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {10:10:80; 110:10:200};
%       [counts, edges] = compute_psth(eventTimes, stimTimes)
%
% Outputs:
%       counts      - counts for each bin
%                   specified as a numeric vector
%       edges       - edges for each bin
%                   specified as a numeric vector
%       relEventTimes   - all relative event times
%                       specified as a numeric vector
%
% Arguments:
%       eventTimes  - event times
%                   must be a numeric vector or a cell array of numeric vectors
%       stimTimes   - stimulus times
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'RelativeTimeWindow': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - 'Grouping': group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == pre- or post- stimulus
%                   - Any other parameter-value pair for compute_grouped_histcounts()
%
% Requires:
%       cd/argfun.m
%       cd/compute_grouped_histcounts.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_psth.m

% File History:
% 2019-09-07 Created by Adam Lu
% 2019-09-08 Added 'Grouping' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
relativeTimeWindowDefault = [];
groupingDefault = [];                   % set later

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
addRequired(iP, 'eventTimes', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['eventTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'stimTimes', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['stimTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RelativeTimeWindow', relativeTimeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Read from the Input Parser
parse(iP, eventTimes, stimTimes, varargin{:});
relativeTimeWindow = iP.Results.RelativeTimeWindow;
grouping = iP.Results.Grouping;

% Keep unmatched arguments for the compute_grouped_histcounts() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force as column vectors
[eventTimes, stimTimes] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            eventTimes, stimTimes);

% Make sure the vector numbers are identical and force as column cell arrays
[eventTimes, stimTimes] = ...
    match_format_vector_sets(eventTimes, stimTimes, 'ForceCellOutputs', true);

% Sort the times in ascending order
[eventTimes, stimTimes] = ...
    argfun(@(x) cellfun(@sort, x, 'UniformOutput', false), ...
            eventTimes, stimTimes);

% Compute the default relative time window for the PSTH
if isempty(relativeTimeWindow)
    % Compute the average inter-stimulus interval
    interStimInterval = compute_average_interval(stimTimes);

    % Use half of the average inter-stimulus interval on each side
    relativeTimeWindow = interStimInterval * 0.5 * [-1, 1];
end

%% Do the job
% Extract relative event times for each window
relEventTimesCellCell = ...
    cellfun(@(x, y) compute_relative_event_times(x, y, relativeTimeWindow), ...
            eventTimes, stimTimes, 'UniformOutput', false);

% Put all relative event times together
relEventTimesCell = vertcat(relEventTimesCellCell{:});
relEventTimes = vertcat(relEventTimesCell{:});

% Create a grouping vector with the pre-stimulus and post-stimulus times
%   as separate groups
%   Note: must be consistent with plot_psth.m
if isempty(grouping)
    grouping = ones(size(relEventTimes));
    grouping(relEventTimes < 0) = -1;
end

% Compute the peri-stimulus time histogram
%   Note: must be consistent with plot_psth.m
[counts, edges] = ...
    compute_grouped_histcounts(relEventTimes, 'Grouping', grouping, ...
                                'FixedEdges', 0, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanInterval = compute_average_interval(eventTimes)
%% Compute the average inter-event interval
% TODO: Pull out as its own function

% Construct a function that calculates the mean over a matrix
nanmeansq = @(x) nanmean(nanmean(x));

% Compute the average inter-event interval
if iscell(eventTimes)
    meanInterval = nanmeansq(cellfun(@(x) nanmeansq(diff(x)), eventTimes));
else
    meanInterval = nanmeansq(diff(eventTimes));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function relEventTimes = compute_relative_event_times(eventTimes, stimTimes, ...
                                                            relativeTimeWindow)
%% Computes the relative event times from event times and stimulus times
% TODO: Pull out as its own function

% Extract a time window for each stimulus time
windows = compute_time_windows(stimTimes, relativeTimeWindow);

% Extract the event times corresponding to each time window
eventTimesEachWindow = extract_subvectors(eventTimes, 'Windows', windows, ...
                                                    'ForceCellOutput', true);

% Compute the relative event times for each time window
relEventTimes = cellfun(@(x, y) x - y, eventTimesEachWindow, ...
                        num2cell(stimTimes), 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windows = compute_time_windows (centerTimes, relativeTimeWindow)
%% Computes time windows based on center times and a relative time window
% TODO: Either merge with compute_time_window.m or pull out as its own function

% Force the relative time window as a column vector
relativeTimeWindow = force_column_vector(relativeTimeWindow);

% Force the center times as a row vector
centerTimes = force_row_vector(centerTimes);

% Compute the time windows
%   Note: Each column corresponds to a time window
windows = repmat(centerTimes, size(relativeTimeWindow)) + ...
            repmat(relativeTimeWindow, size(centerTimes));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%