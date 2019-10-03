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
%                   - 'StimDuration': stimulus duration for plotting
%                                       (stim always occur at 0)
%                   must be a positive scalar
%                   default == [] (not plotted)
%                   - 'Grouping': group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == pre- or post- stimulus
%                   - 'Edges': bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == automatic detection of 
%                   - Any other parameter-value pair for compute_grouped_histcounts()
%
% Requires:
%       cd/adjust_window_to_bounds.m
%       cd/compute_grouped_histcounts.m
%       cd/compute_relative_event_times.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_psth.m

% File History:
% 2019-09-07 Created by Adam Lu
% 2019-09-08 Added 'Grouping' as an optional argument
% 2019-09-15 Now trims the stimulus window so that 
%               it does not exceed relativeTimeWindow
% TODO: Add option to shift relative event times by stimDelay


%% Hard-coded parameters
% TODO: Make optional parameters
stimStart = 0;
stimWindow = [];

%% Default values for optional arguments
relativeTimeWindowDefault = [];
stimDurationDefault = 0;
groupingDefault = [];                   % set later
edgesDefault = [];              % set later

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
addParameter(iP, 'StimDuration', stimDurationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Read from the Input Parser
parse(iP, eventTimes, stimTimes, varargin{:});
relativeTimeWindow = iP.Results.RelativeTimeWindow;
stimDuration = iP.Results.StimDuration;
grouping = iP.Results.Grouping;
edges = iP.Results.Edges;

% Keep unmatched arguments for the compute_grouped_histcounts() function
otherArguments = iP.Unmatched;

%% Preparation
% Decide on stimulus window
if isempty(stimWindow)
    stimWindow = stimStart + [0, stimDuration];
end

% Force as column vectors
[eventTimes, stimTimes] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            eventTimes, stimTimes);

% Make sure the vector numbers are identical and force as column cell arrays
[eventTimes, stimTimes] = ...
    match_format_vector_sets(eventTimes, stimTimes, 'ForceCellOutputs', true);

%% Do the job
% Compute relative event times
[relEventTimesCellCell, relativeTimeWindow] = ...
    compute_relative_event_times(eventTimes, stimTimes, ...
                                    'RelativeTimeWindow', relativeTimeWindow);

% Trim the stimulus window so that it does not exceed relativeTimeWindow
%   bounds
stimWindow = adjust_window_to_bounds(stimWindow, relativeTimeWindow);

% Put all relative event times together
% TODO: Change this
relEventTimesCell = vertcat(relEventTimesCellCell{:});
relEventTimes = vertcat(relEventTimesCell{:});

% Create a grouping vector with the pre-stimulus and post-stimulus times
%   as separate groups
% TODO: How to assign default grouping
%   Note: must be consistent with plot_psth.m
if isempty(grouping)
    grouping = ones(size(relEventTimes));
    grouping(relEventTimes < 0) = -1;
end

% Compute the peri-stimulus time histogram
%   Note: must be consistent with plot_psth.m
[counts, edges] = ...
    compute_grouped_histcounts(relEventTimes, 'Grouping', grouping, ...
                                'Edges', edges, 'FixedEdges', stimWindow, ...
                                otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%