function [relEventTimes, relativeTimeWindow] = ...
                compute_relative_event_times (eventTimes, stimTimes, varargin)
%% Computes the relative event times from event times and stimulus times
% Usage: [relEventTimes, relativeTimeWindow] = ...
%               compute_relative_event_times (eventTimes, stimTimes, varargin)
% Explanation:
%       Computes the relative event times around stimulus times
%       For each stimulus time, a window around it is used
%       Note: Output is:
%               - a numeric vector if eventTimes is a vector 
%                       and stimTimes is a scalar
%               - a cell array of numeric vectors if eventTimes is a vector 
%                       and stimTimes is a vector
%               - a cell array of cell arrays of numeric vectors 
%                       if either eventTimes or stimTimes is a cell array
%
% Example(s):
%       eventTimes = randi(100, 100, 1);
%       stimTimes = 50;
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = randi(100, 100, 1);
%       stimTimes = 10:10:80;
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {50; 150};
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {10:10:80; 110:10:200};
%       relEventTimes = compute_relative_event_times(eventTimes, stimTimes)
%
% Outputs:
%       relEventTimes   - relative event times
%                       specified as a numeric vector,
%                           a cell array of numeric vectors or
%                           a cell array of cell arrays of numeric vectors
%
% Arguments:
%       eventTimes  - event times
%                       Note: If a cell array, each cell 
%                               will be matched with stimTimes
%                   must be a numeric vector or a cell array of numeric vectors
%       stimTimes   - stimulus times
%                       Note: If a cell array, each cell 
%                               will be matched with eventTimes
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'RelativeTimeWindow': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/force_row_vector.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%       cd/plot_relative_events.m
%
% Used by:
%       cd/compute_psth.m

% File History:
% 2019-09-11 Moved from compute_psth.m
% 2019-09-15 Now returns relative event time window used
% TODO: Add option to shift relative event times by stimDelay
% 

%% Hard-coded parameters

%% Default values for optional arguments
relativeTimeWindowDefault = [];

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

% Read from the Input Parser
parse(iP, eventTimes, stimTimes,  varargin{:});
relativeTimeWindow = iP.Results.RelativeTimeWindow;

%% Preparation
% Force as column vectors
[eventTimes, stimTimes] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            eventTimes, stimTimes);

% Make sure the vector numbers are identical and force as column cell arrays
[eventTimesCell, stimTimesCell] = ...
    match_format_vector_sets(eventTimes, stimTimes, 'ForceCellOutputs', true);

% Sort the times in ascending order
[eventTimesCell, stimTimesCell] = ...
    argfun(@(x) cellfun(@sort, x, 'UniformOutput', false), ...
            eventTimesCell, stimTimesCell);

% Compute the default relative time window
if isempty(relativeTimeWindow)
    % Compute the average inter-stimulus interval
    interStimInterval = compute_average_interval(stimTimes);

    if isnan(interStimInterval)
        % Compute the minimum interval that contains 
        %   all event times for each set
        minHalfWidthEachCell = ...
            cellfun(@(x, y) max(abs([max(x) - y, min(x) - y])), ...
                    eventTimesCell, stimTimesCell);

        % Use an interval that works for all sets
        windowHalfWidth = max(minHalfWidthEachCell);
    else
        % Use half of the average inter-stimulus interval on each side
        windowHalfWidth = interStimInterval * 0.5;
    end

    relativeTimeWindow = windowHalfWidth * [-1, 1];
end

%% Do the job
% Extract relative event times for each window
relEventTimesCellCell = ...
    cellfun(@(x, y) compute_relative_event_times_helper(x, y, ...
                                                    relativeTimeWindow), ...
            eventTimesCell, stimTimesCell, 'UniformOutput', false);

%% Output results
if ~iscell(eventTimes) && ~iscell(stimTimes)
    % Extract from the cell array if there is only one set of event times
    %   and one set of stim times
    if numel(relEventTimesCellCell) == 1
        relEventTimesCell = relEventTimesCellCell{1};

        if numel(relEventTimesCell) == 1
            relEventTimes = relEventTimesCell{1};
        else
            relEventTimes = relEventTimesCell;
        end
    else
        error('Not implemented yet!')
    end
else
    relEventTimes = relEventTimesCellCell;
end

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

function relEventTimes = ...
                compute_relative_event_times_helper(eventTimes, stimTimes, ...
                                                            relativeTimeWindow)
%% Computes the relative event times from event times and stimulus times

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
