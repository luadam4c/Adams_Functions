function handles = compare_events_pre_post_stim (varargin)
%% Plots events (such as SWDs) relative to stim (such as gas pulses) (unfinished)
% Usage: handles = compare_events_pre_post_stim (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'raster'    - event raster
%                       'psth'      - peri-stimulus time histogram
%                       'chevron'   - chevron plot
%                   default == 'raster'
%                   - 'StimIndices': stimulation indices to restrict to
%                   must be a positive integer array or string recognized by
%                        the 'Pattern' option of extract_subvectors.m
%                           'odd'   - odd indices
%                           'even'  - even indices
%                   default == no restrictions
%                   - 'EventTableSuffix': Suffix for the event table
%                   must be a character vector or a string scalar 
%                   default == '_SWDs'
%                   - 'StimTableSuffix': Suffix for the stim table
%                   must be a character vector or a string scalar 
%                   default == '_pulses'
%                   - 'Directory': directory to look for event table
%                                   and stim table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'RelativeTimeWindowMin': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - 'StimDurationMin': stimulus duration for plotting
%                                       (stim always occur at 0)
%                   must be a positive scalar
%                   default == [] (not plotted)
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'YLimitsLog2Ratio': limits of y axis 
%                                           for the log2 ratio plot
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == TODO
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == TODO
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc2'}
%                   - Any other parameter-value pair for plot_violin()
%
% Requires:
%       cd/argfun.m
%       cd/apply_to_all_cells.m
%       cd/compute_relative_event_times.m
%       cd/extract_subvectors.m
%       cd/load_matching_sheets.m
%       cd/plot_violin.m
%   TODO
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/create_label_from_sequence.m
%       cd/set_figure_properties.m
%       cd/save_all_figtypes.m
%
% Used by:
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2020-06-26 Modified from plot_relative_events.m
% 

%% Hard-coded parameters
SEC_PER_MIN = 60;
validPlotTypes = {'raster', 'psth', 'chevron'};
plotLog2Ratio = true;

% TODO: Make optional arguments
stimStartLineColor = [0.5, 0.5, 0.5];
stimStartLineWidth = 1;
pathBase = '';
sheetType = 'csv';
figSuffix = '';
labels = {};

%% Default values for optional arguments
plotTypeDefault = 'raster';
stimIndicesDefault = [];        % take all stims by default
eventTableSuffixDefault = '_SWDs';
stimTableSuffixDefault = '_pulses';
directoryDefault = '';          % set later
relativeTimeWindowMinDefault = [];
stimDurationMinDefault = [];
yLimitsDefault = [];            % set later
yLimitsLog2RatioDefault = [];  % set later
figTitleDefault = '';           % set later
figNameDefault = '';            % set later
figTypesDefault = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'StimIndices', stimIndicesDefault, ...
    @(x) validateattributes(x, {'numeric', 'char'}, {'2d'}));
addParameter(iP, 'EventTableSuffix', eventTableSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'StimTableSuffix', stimTableSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RelativeTimeWindowMin', relativeTimeWindowMinDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'StimDurationMin', stimDurationMinDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimitsLog2Ratio', yLimitsLog2RatioDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
stimIndices = iP.Results.StimIndices;
eventTableSuffix = iP.Results.EventTableSuffix;
stimTableSuffix = iP.Results.StimTableSuffix;
directory = iP.Results.Directory;
relTimeWindowMin = iP.Results.RelativeTimeWindowMin;
avgStimDurationMin = iP.Results.StimDurationMin;
yLimits = iP.Results.YLimits;
yLimitsLog2Ratio = iP.Results.YLimitsLog2Ratio;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_violin() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output
handles = struct;

% Set default directory
if isempty(directory)
    directory = pwd;
end

% Set default figure title
if isempty(figTitle)
    if isempty(stimIndices)
        figTitle = 'Event property around all stims';
    elseif isnumeric(stimIndices)
        figTitle = ['Event property around stims ', ...
                    create_label_from_sequence(stimIndices)];
    elseif ischar(stimIndices)
        figTitle = ['Event property around ', stimIndices, ' stims'];
    end
end

%% Get relative event times
% Load matching stimulus and event tables
[stimTables, swdTables, distinctParts] = ...
    load_matching_sheets(stimTableSuffix, eventTableSuffix, ...
                        'Directory', directory);

% Set default labels for each raster
if isempty(labels)
    labels = distinctParts;
end

% Extract all start times in seconds
[stimStartTimesSec, swdStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Extract durations in seconds
[stimDurationsSec, swdDurationsSec] = ...
    argfun(@(x) cellfun(@(y) y.duration, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Restrict to certain stimulation windows if requested
if isempty(stimIndices)
    % Do nothing
elseif isnumeric(stimIndices)
    [stimStartTimesSec, stimDurationsSec] = ...
        argfun(@(x) extract_subvectors(x, 'Indices', stimIndices), ...
                stimStartTimesSec, stimDurationsSec);
elseif ischar(stimIndices)
    [stimStartTimesSec, stimDurationsSec] = ...
        argfun(@(x) extract_subvectors(x, 'Pattern', stimIndices), ...
                stimStartTimesSec, stimDurationsSec);
end

% Convert to minutes
[stimStartTimesMin, stimDurationsMin, ...
        swdStartTimesMin, swdDurationsMin] = ...
    argfun(@(x) cellfun(@(y) y / SEC_PER_MIN, x, 'UniformOutput', false), ...
            stimStartTimesSec, stimDurationsSec, ...
            swdStartTimesSec, swdDurationsSec);

%% Classify events
% Extract events of interest and compute the relative event times
%   Note: this should return a cell matrix:
%           Each column is a file
%           Each row is a stim
[relEventTimesMin, relTimeWindowMin, origInd] = ...
    compute_relative_event_times(swdStartTimesMin, stimStartTimesMin, ...
                            'RelativeTimeWindow', relTimeWindowMin, ...
                            'ForceMatrixOutput', true);

% Extract event information for the events of interest
eventDurationsSec = extract_subvectors(swdDurationsSec, 'Indices', origInd);

% Determine whether each event is after the stimulus
isAfterStim = apply_to_all_cells(@(x) x > 0, relEventTimesMin);

% Pool all events together
[isAfterStimPooled, relEventTimesMinPooled, eventDurationsSecPooled] = ...
    argfun(@(x) apply_over_cells(@vertcat, x), ...
            isAfterStim, relEventTimesMin, eventDurationsSec);

% Extract events before stim
[relEventTimesMinBefore, eventDurationsSecBefore] = ...
    argfun(@(x) x(~isAfterStimPooled), ...
            relEventTimesMinPooled, eventDurationsSecPooled);

% Extract events after stim
[relEventTimesMinAfter, eventDurationsSecAfter] = ...
    argfun(@(x) x(isAfterStimPooled), ...
            relEventTimesMinPooled, eventDurationsSecPooled);

%% Compare event properties
% Create figure
if ~isempty(figName)
    fig = set_figure_properties('AlwaysNew', true);
end

% Hold on
wasHold = hold_on;

% Place groups together
twoGroups = {eventDurationsSecBefore, eventDurationsSecAfter};
xTickLabels = {'Before', 'After'};

% Plot violin plot
violins = plot_violin({eventDurationsSecBefore, eventDurationsSecAfter}, ...
                'XTickLabels', xTickLabels, 'YLabel', 'Event Duration (sec)');

% Update y axis limits to include 0
yLimitsOrig = get(gca, 'YLim');
yLimitsNew = yLimitsOrig;
yLimitsNew(1) = 0;
set(gca, 'YLim', yLimitsNew);

% TODO: plot_significance.m
% Test difference
statsStruct = test_difference(twoGroups, 'IsPaired', false);
symbol = statsStruct.symbol;
pValue = statsStruct.pValue;
testFunction = statsStruct.testFunction;

% Plot text for difference
xLimitsOrig = get(gca, 'XLim');
yLimitsOrig = get(gca, 'YLim');
text(0.5, 0.9, symbol, 'Units', 'normalized', 'HorizontalAlignment', 'center');
text(0.5, 0.85, sprintf('p_{%s} = %g', testFunction, pValue), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center');

% Create title
title(figTitle);

% Hold off
hold_off(wasHold);

% Save figure
if ~isempty(figName)
    save_all_figtypes(fig, figName, figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
