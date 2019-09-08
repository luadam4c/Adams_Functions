function [counts, edges] = compute_psth (eventTimes, stimTimes, varargin)
%% Computes a peri-stimulus time histogram
% Usage: [counts, edges] = compute_psth (eventTimes, stimTimes, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       eventTimes = {1:2:99; 101:2:199};
%       stimTimes = {10:10:80; 110:10:200};
%       [counts, edges] = compute_psth(eventTimes, stimTimes)
%
% Outputs:
%       counts      - counts for each bin
%                   specified as a numeric vector
%       edges       - edges for each bin
%                   specified as a numeric vector
%
% Arguments:
%       eventTimes  - event times
%                   must be a numeric vector or a cell array of numeric vectors
%       stimTimes   - stimulus times
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'RelativeTimeWindow': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - Any other parameter-value pair for histcounts()
%
% Requires:
%       cd/create_error_for_nargin.m
%   argfun
% force_column_vector
% force_row_vector
% extract_subvectors
%   iscellnumeric
%       cd/match_format_vector_sets.m
%       cd/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-07 Created by Adam Lu
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
parse(iP, eventTimes, stimTimes, varargin{:});
relativeTimeWindow = iP.Results.RelativeTimeWindow;

% Keep unmatched arguments for the histcounts() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force as column vectors
[eventTimes, stimTimes] = argfun(@force_column_vector, eventTimes, stimTimes);

% Make sure the vector numbers are identical and force as column cell arrays
[eventTimes, stimTimes] = ...
    match_format_vector_sets(eventTimes, stimTimes, 'ForceCellOutputs', true);

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

% Compute the peri-stimulus time histogram
[counts, edges] = histcounts(relEventTimes, otherArguments{:});

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

function relEventTimes = compute_relative_event_times(eventTimes, stimTimes, ...
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

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%