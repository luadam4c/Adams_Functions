function handles = plot_psth (varargin)
%% Plots a peri-stimulus time histogram
% Usage: handles = plot_psth (relEventTimes (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [counts, edges] = compute_psth(randi(100, 100, 1), 10:10:80);
%       handles = plot_psth('Counts', counts, 'Edges', edges);
%       eventTimes = {randi(100, 100, 1); randi(100, 100, 1) + 100};
%       stimTimes = {10:10:80; 110:10:200};
%       handles = plot_psth('EventTimes', eventTimes, 'StimTimes', stimTimes);
%       handles = plot_psth('EventTimes', eventTimes, 'StimTimes', stimTimes, 'StimDuration', 3);
%
% Outputs:
%
% Arguments:
%       relEventTimes   - relative event times
%                       must be a vector of one the following types:
%                           'numeric', 'logical', 'datetime', 'duration'
%                       default == not provided
%       varargin    - 'Counts': bin counts, with each group 
%                                   being a different column
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_psth(eventTimes, stimTimes)
%                   - 'Edges': bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_psth(eventTimes, stimTimes)
%                   - 'EventTimes': event times
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == not provided
%                   - 'StimTimes': stimulus times
%                   must be a numeric vector or a cell array of numeric vectors
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
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'Relative Time From Stim'
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Event Count'
%                   - 'GroupingLabels': labels for the groupings, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Pre-Stim', 'Post-Stim'}
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Peri-Stimulus Time Histogram on ', ...
%                               create_time_stamp]
%                   - Any other parameter-value pair for plot_histogram()
%
% Requires:
%       cd/compute_grouped_histcounts.m
%       cd/compute_psth.m
%       cd/create_error_for_nargin.m
%       cd/create_time_stamp.m
%       cd/plot_histogram.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-07 Created by Adam Lu
% 2019-09-08 Added 'Grouping' as an optional argument
% 2019-09-08 Added 'StimDuration' as an optional argument
% 2019-09-08 Now plots vertical shade

%% Hard-coded parameters
% TODO: Make optional parameters
vertLineLineWidth = 2;
vertLineLineStyle = '-';
vertLineColor = 'k';
stimStart = 0;
stimWindow = [];

%% Default values for optional arguments
relEventTimesDefault = [];              % set later
countsDefault = [];                     % set later
edgesDefault = [];                      % set later
eventTimesDefault = [];
stimTimesDefault = [];
relativeTimeWindowDefault = [];
stimDurationDefault = 0;
groupingDefault = [];                   % set later
xLabelDefault = 'Relative Time From Stim';
yLabelDefault = 'Event Count';
groupingLabelsDefault = {};             % set later
figTitleDefault = '';                   % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'relEventTimes', relEventTimesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Counts', countsDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'EventTimes', eventTimesDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['eventTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'StimTimes', stimTimesDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['stimTimes must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'RelativeTimeWindow', relativeTimeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'StimDuration', stimDurationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
relEventTimes = iP.Results.relEventTimes;
counts = iP.Results.Counts;
edges = iP.Results.Edges;
eventTimes = iP.Results.EventTimes;
stimTimes = iP.Results.StimTimes;
relativeTimeWindow = iP.Results.RelativeTimeWindow;
stimDuration = iP.Results.StimDuration;
grouping = iP.Results.Grouping;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
groupingLabels = iP.Results.GroupingLabels;
figTitle = iP.Results.FigTitle;

% Keep unmatched arguments for the plot_histogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output structure
handles = struct;

% Make sure there is data
if isempty(relEventTimes) && (isempty(counts) || isempty(edges)) && ...
        (isempty(eventTimes) || isempty(stimTimes))
    disp('There is no data to plot!');
    return
end

% Decide on stimulus window
if isempty(stimWindow)
    stimWindow = stimStart + [0, stimDuration];
end

% Compute histogram counts if not already done
if isempty(counts)
    if isempty(relEventTimes)
        % Compute counts and edges from eventTimes and stimTimes
        [counts, edges] = ...
            compute_psth(eventTimes, stimTimes, ...
                        'StimDuration', stimDuration, ...
                        'RelativeTimeWindow', relativeTimeWindow);
    else
        % Compute a grouping vector from relEventTimes
        %   Note: must be consistent with compute_psth.m
        if isempty(grouping)
            grouping = ones(size(relEventTimes));
            grouping(relEventTimes < 0) = -1;
        end

        % Compute counts and edges from relEventTimes
        %   Note: must be consistent with compute_psth.m
        [counts, edges] = ...
            compute_grouped_histcounts(relEventTimes, 'Grouping', grouping, ...
                                'FixedEdges', stimWindow, otherArguments{:});
    end
else
    if isempty(edges)
        disp('Edges must be provided if counts are provided!');
        return
    end
end

% Create default grouping labels
if isempty(groupingLabels)
    if size(counts, 2) <= 2
        groupingLabels = {'Pre-Stim', 'Post-Stim'};
    end
end

% Create default title
if isempty(figTitle)
    figTitle = ['Peri-Stimulus Time Histogram on ', create_time_stamp];
end

%% Do the job
% Plot the histogram
[bars, fig] = plot_histogram('Counts', counts, 'Edges', edges, ...
                            'XLabel', xLabel, 'YLabel', yLabel, ...
                            'GroupingLabels', groupingLabels, ...
                            'FigTitle', figTitle, ...
                            otherArguments);

% Plot the stimulus start as a vertical line
vertLine = plot_vertical_line(stimStart, 'LineWidth', vertLineLineWidth, ...
                                'LineStyle', vertLineLineStyle, ...
                                'Color', vertLineColor);

% Plot stimulus window as a vertical shade
vertLine = plot_vertical_shade(stimWindow);

%% Output handles
handles.fig = fig;
handles.bars = bars;
handles.vertLine = vertLine;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%