function [counts, edges, relativeEventTimes] = ...
                compute_psth (eventTimes, stimTimes, varargin)
%% Computes a peri-stimulus time histogram
% Usage: [counts, edges, relativeEventTimes] = ...
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
%       relativeEventTimes   - all relative event times
%                       specified as a numeric vector
%
% Arguments:
%       eventTimes  - event times
%                   must be a numeric vector or a cell array of numeric vectors
%       stimTimes   - stimulus times
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'RelativeEventTimes': relative event times
%                   must be a numeric vector or a cell array of numeric vectors
%                       or a cell array of cell arrays of numeric vectors
%                   default == not provided
%                   - 'RelativeTimeWindow': relative time window
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
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%       cd/iscellnumeric.m
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
relativeEventTimesDefault = [];
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
    @(x) assert(isempty(x) || isnum(x) || iscellnumeric(x), ...
                ['eventTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'stimTimes', ...
    @(x) assert(isempty(x) || isnum(x) || iscellnumeric(x), ...
                ['stimTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RelativeEventTimes', relativeEventTimesDefault, ...
    @(x) assert(isempty(x) || isnum(x) || iscellnumeric(x), ...
                ['RelativeEventTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'RelativeTimeWindow', relativeTimeWindowDefault, ...
    @(x) isempty(x) || isnum(x));
addParameter(iP, 'StimDuration', stimDurationDefault, ...
    @(x) isempty(x) || isnum(x));
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Read from the Input Parser
parse(iP, eventTimes, stimTimes, varargin{:});
relativeEventTimes = iP.Results.RelativeEventTimes;
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

% Decide on relative time window based on given edges of the histogram
if isempty(relativeTimeWindow) && ~isempty(edges)
    relativeTimeWindow = [edges(1), edges(end)];
end

%% Do the job
% Compute relative event times
%   Note: If eventTimes is a cell array, 
%           this will be a cell array of cell arrays;
%         Otherwise, this will be a cell array of numeric vectors (many stims)
%           or a numeric vector (one stim time)
if isempty(relativeEventTimes)
    [relativeEventTimes, relativeTimeWindow] = ...
        compute_relative_event_times(eventTimes, stimTimes, ...
                                    'RelativeTimeWindow', relativeTimeWindow);
end

% If relative time window is still empty, estimate it from data
if isempty(relativeTimeWindow)
    % Find the maximum and minimum times
    minTime = apply_iteratively(@min, relativeEventTimes);
    maxTime = apply_iteratively(@max, relativeEventTimes);

    % Construct the relative time window
    relativeTimeWindow = [minTime, maxTime];
end

% Trim the stimulus window so that it does not exceed relativeTimeWindow
%   bounds
stimWindow = adjust_window_to_bounds(stimWindow, relativeTimeWindow);

% Pool relative event times across stims for each file
if iscellnumeric(relativeEventTimes)
    % Many stims for one condition, just pool
    relativeEventTimes = vertcat(relativeEventTimes{:});
elseif iscell(relativeEventTimes)
    % Pool for each file
    relativeEventTimes = cellfun(@(x) vertcat(x{:}), ...
                            relativeEventTimes, 'UniformOutput', false);    
end

% Decide on the grouping vector
if isempty(grouping)
    if iscell(relativeEventTimes)
        % Create a grouping vector with each file as separate groups
        grouping = create_default_grouping('Stats', relativeEventTimes);
    else
        % Create a grouping vector with the pre-stimulus and post-stimulus times
        %   as separate groups
        grouping = ones(size(relativeEventTimes));
        grouping(relativeEventTimes < 0) = -1;
    end
end

% Compute the peri-stimulus time histogram
[counts, edges] = ...
    compute_grouped_histcounts(relativeEventTimes, 'Grouping', grouping, ...
                                'Edges', edges, 'FixedEdges', stimWindow, ...
                                otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%