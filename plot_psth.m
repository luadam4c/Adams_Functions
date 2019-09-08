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
%                   - Any other parameter-value pair for plot_histogram()
%
% Requires:
%       cd/compute_bins.m
%       cd/compute_psth.m
%       cd/create_error_for_nargin.m
%       cd/plot_histogram.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-07 Created by Adam Lu
% TODO: Add stimulus duration and vertical shade

%% Hard-coded parameters

%% Default values for optional arguments
relEventTimesDefault = [];              % set later
countsDefault = [];                     % set later
edgesDefault = [];                      % set later
eventTimesDefault = [];
stimTimesDefault = [];
relativeTimeWindowDefault = [];

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

% Read from the Input Parser
parse(iP, varargin{:});
relEventTimes = iP.Results.relEventTimes;
counts = iP.Results.Counts;
edges = iP.Results.Edges;
eventTimes = iP.Results.EventTimes;
stimTimes = iP.Results.StimTimes;
relativeTimeWindow = iP.Results.RelativeTimeWindow;

% Keep unmatched arguments for the plot_histogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Make sure there is data
if isempty(relEventTimes) && (isempty(counts) || isempty(edges)) && ...
        (isempty(eventTimes) || isempty(stimTimes))
    handles = struct;
    disp('There is no data to plot!');
    return
end

% Compute histogram if not already done
if isempty(counts) || isempty(edges)
    if isempty(relEventTimes)
        % Compute counts and edges from eventTimes and stimTimes
        [counts, edges] = compute_psth(eventTimes, stimTimes, ...
                                    'RelativeTimeWindow', relativeTimeWindow);
    else
        % Compute counts and edges from relEventTimes
        [counts, edges] = compute_bins(relEventTimes, 'FixedEdges', 0);
    end
end

%% Do the job
% Plot the histogram
[bars, fig] = plot_histogram('Counts', counts, 'Edges', edges, otherArguments);

% Plot a vertical line at zero
vertLine = plot_vertical_line(0, 'LineWidth', 2, 'LineStyle', '-', ...
                                'Color', 'r');

%% Output handles
handles.fig = fig;
handles.bars = bars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%