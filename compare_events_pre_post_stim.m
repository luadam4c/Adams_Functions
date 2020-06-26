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
%       varargin    - 'PropertyName': property name to compare
%                   must be a string scalar or a character vector
%                   default == 'duration'
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
%       cd/create_label_from_sequence.m
%       cd/extract_subvectors.m
%       cd/hold_on.m
%       cd/hold_off.m
%       cd/read_matching_sheets.m
%       cd/plot_violin.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_difference.m
%
% Used by:
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2020-06-26 Modified from plot_relative_events.m
% 

%% Hard-coded parameters
SEC_PER_MIN = 60;

%% Default values for optional arguments
propertyNameDefault = 'duration';
stimIndicesDefault = [];        % take all stims by default
eventTableSuffixDefault = '_SWDs';
stimTableSuffixDefault = '_pulses';
directoryDefault = '';          % set later
relativeTimeWindowMinDefault = [];
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
addParameter(iP, 'PropertyName', propertyNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
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
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
propertyName = iP.Results.PropertyName;
stimIndices = iP.Results.StimIndices;
eventTableSuffix = iP.Results.EventTableSuffix;
stimTableSuffix = iP.Results.StimTableSuffix;
directory = iP.Results.Directory;
relTimeWindowMin = iP.Results.RelativeTimeWindowMin;
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

% Decide on property label
switch propertyName
    case 'duration'
        propertyLabel = 'Event Duration (sec)';
    case 'peakFrequency'
        propertyLabel = 'Peak Frequency (Hz)';
    otherwise
        propertyLabel = propertyName;
end

% Set default figure title
if isempty(figTitle)
    if isempty(stimIndices)
        figTitle = [propertyLabel, ' around all stims'];
    elseif isnumeric(stimIndices)
        figTitle = [propertyLabel, ' around stims #', ...
                    create_label_from_sequence(stimIndices)];
    elseif ischar(stimIndices)
        figTitle = [propertyLabel ' around ', stimIndices, ' stims'];
    end
end

%% Get relative event times
% Load matching stimulus and event tables
[stimTables, eventTables, distinctParts] = ...
    read_matching_sheets(stimTableSuffix, eventTableSuffix, ...
                        'Directory', directory);

% Extract all start times in seconds
[stimStartTimesSec, eventStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, eventTables);

% Extract event properties
if ~is_field(eventTables{1}, propertyName)
    fprintf('There is no column named ''%s''!!\n', propertyName);
    handles = struct;
    return
else
    eventProperties = ...
        cellfun(@(y) y.(propertyName), eventTables, 'UniformOutput', false);
end

% Restrict to certain stimulation windows if requested
if isempty(stimIndices)
    % Do nothing
elseif isnumeric(stimIndices)
    stimStartTimesSec = ...
        extract_subvectors(stimStartTimesSec, 'Indices', stimIndices);
elseif ischar(stimIndices)
    stimStartTimesSec = ...
        extract_subvectors(stimStartTimesSec, 'Pattern', stimIndices);
end

% Convert to minutes
[stimStartTimesMin, swdStartTimesMin] = ...
    argfun(@(x) cellfun(@(y) y / SEC_PER_MIN, x, 'UniformOutput', false), ...
            stimStartTimesSec, eventStartTimesSec);

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
eventProperties = extract_subvectors(eventProperties, 'Indices', origInd);

% Determine whether each event is after the stimulus
isAfterStim = apply_to_all_cells(@(x) x > 0, relEventTimesMin);

% Pool all events together
[isAfterStimPooled, relEventTimesMinPooled, eventPropertiesPooled] = ...
    argfun(@(x) apply_over_cells(@vertcat, x), ...
            isAfterStim, relEventTimesMin, eventProperties);

% Extract events before stim
[relEventTimesMinBefore, eventPropertiesBefore] = ...
    argfun(@(x) x(~isAfterStimPooled), ...
            relEventTimesMinPooled, eventPropertiesPooled);

% Extract events after stim
[relEventTimesMinAfter, eventPropertiesAfter] = ...
    argfun(@(x) x(isAfterStimPooled), ...
            relEventTimesMinPooled, eventPropertiesPooled);

%% Compare event properties
% Create figure
if ~isempty(figName)
    fig = set_figure_properties('AlwaysNew', true);
end

% Hold on
wasHold = hold_on;

% Place groups together
twoGroups = {eventPropertiesBefore, eventPropertiesAfter};
xTickLabels = {'Before', 'After'};

% Plot violin plot
violins = plot_violin({eventPropertiesBefore, eventPropertiesAfter}, ...
                'XTickLabels', xTickLabels, 'YLabel', propertyLabel, ...
                otherArguments);

% Update y axis limits to include 0
switch propertyName
    case 'duration'
        yLimitsOrig = get(gca, 'YLim');
        yLimitsNew = yLimitsOrig;
        yLimitsNew(1) = 0;
        set(gca, 'YLim', yLimitsNew);
    otherwise
        % Do nothing
end

% TODO: plot_significance.m
% Test difference
statsStruct = test_difference(twoGroups, 'IsPaired', false);
symbol = statsStruct.symbol;
pValue = statsStruct.pValue;
testFunction = statsStruct.testFunction;

% Plot text for difference
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
