function handles = plot_relative_events (varargin)
%% Plots events (such as SWDs) relative to stim (such as gas pulses) (unfinished)
% Usage: handles = plot_relative_events (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_relative_events('Directory', '/media/shareX/2019octoberR01/Figures/Figure1c')
%       plot_relative_events('RelativeTimeWindow', [-20, 20]);
%       plot_relative_events('PlotType', 'psth', 'Edges', -20:2:20);
%       plot_relative_events('RelativeTimeWindow', [-20, 20], 'PlotType', 'chevron', 'PlotErrorBars', true);
%       plot_relative_events('RelativeTimeWindow', [-15, 15]);
%       plot_relative_events('PlotType', 'psth', 'Edges', -15:2.5:15, 'StimIndices', 'odd');
%       plot_relative_events('PlotType', 'psth', 'Edges', -15:2.5:15, 'StimIndices', 'even');
%       plot_relative_events('RelativeTimeWindow', [-15, 15], 'PlotType', 'chevron');
%       plot_relative_events('RelativeTimeWindow', [-15, 15], 'PlotType', 'chevron', 'StimIndices', 'odd');
%       plot_relative_events('RelativeTimeWindow', [-15, 15], 'PlotType', 'chevron', 'StimIndices', 'even');
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
%                   default == 'SWDs'
%                   - 'StimTableSuffix': Suffix for the stim table
%                   must be a character vector or a string scalar 
%                   default == 'pulses'
%                   - 'CorrTableSuffix': Suffix for the correlation table
%                   must be a character vector or a string scalar 
%                   default == 'resp_table'
%                   - 'CorrMeasures': Variable names in the correlation table
%                                       to correlate
%                   must be a character vector or a string scalar 
%                   default == {'respRates'; 'respAmps'}
%                   - 'Directory': directory to look for event table
%                                   and stim table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'RelativeTimeWindowMin': relative time window
%                   must be a 2-element numeric vector
%                   default == interStimInterval * 0.5 * [-1, 1]
%                   - 'BeforeWindowMin': time window before stim
%                   must be empty or a 2-element numeric vector
%                   default == set in count_events.m
%                   - 'AfterWindowMin': time window after stim
%                   must be empty or a 2-element numeric vector
%                   default == set in count_events.m
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
%                   - Any other parameter-value pair for plot_raster() 
%                           or plot_chevron() or plot_psth() 
%
% Requires:
%       cd/argfun.m
%       cd/apply_iteratively.m
%       cd/compute_relative_event_times.m
%       cd/count_events.m
%       cd/create_label_from_sequence.m
%       cd/create_subplots.m
%       cd/extract_subvectors.m
%       cd/extract_fileparts.m
%       cd/extract_vars.m
%       cd/force_column_cell.m
%       cd/read_matching_sheets.m
%       cd/plot_chevron.m
%       cd/plot_psth.m
%       cd/plot_raster.m
%       cd/plot_vertical_shade.m
%
% convert_units
% create_shifted_vectors
% find_matching_files
% plot_grouped_scatter
% plot_correlation_coefficient
% save_all_figtypes
% set_figure_properties
% update_figure_for_corel
%
% Used by:
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-10 Created by Adam Lu
% 2019-09-11 Added 'PlotType' as an optional argument
% 2019-09-15 Added 'StimTableSuffix' and 'EventTableSuffix' 
%               as optional arguments
% 2019-09-25 Finished the raster plot code
% 2019-09-30 Now uses read_matching_sheets.m
% 2019-10-04 Added 'StimIndices' as an optional arguments
% 2019-10-06 Added 'YLimits' as an optional argument
% 2019-10-06 Now plots SWD count ratio instead of percentage
% 2019-10-06 Added 'YLimitsLog2Ratio' as an optional argument
% 2019-10-08 Now plot log2 data, then change y tick labels to reflect original values
% 2019-10-08 Now conduct t-test on log2 data
% 2020-07-23 Added stimulation window to raster plot
% 2020-08-18 Added 'BeforeWindowMin' & 'AfterWindowMin' as optional arguments
% 2020-08-18 Now writes row names for Chevron tables
% 2020-08-19 Added 'corr' to valid plot types

%% Hard-coded parameters
SEC_PER_MIN = 60;
validPlotTypes = {'raster', 'psth', 'chevron', 'corr'};
plotLog2Ratio = true;

% TODO: Make optional arguments
stimStartLineColor = [0.5, 0.5, 0.5];
stimStartLineWidth = 1;
pathBase = '';
sheetType = 'csv';
figSuffix = '';
labels = {};
timeStr = 'timeInstants';

%% Default values for optional arguments
plotTypeDefault = 'raster';
stimIndicesDefault = [];        % take all stims by default
eventTableSuffixDefault = 'SWDs';
stimTableSuffixDefault = 'pulses';
corrTableSuffixDefault = 'resp_table';
corrMeasuresDefault = {'respRates'; 'respAmps'};
directoryDefault = '';          % set later
relativeTimeWindowMinDefault = [];
beforeWindowMinDefault = [];
afterWindowMinDefault = [];
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
addParameter(iP, 'CorrTableSuffix', corrTableSuffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'CorrMeasures', corrMeasuresDefault, ...
    @(x) validateattributes(x, {'cell', 'char', 'string'}, {'2d'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RelativeTimeWindowMin', relativeTimeWindowMinDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'BeforeWindowMin', beforeWindowMinDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'AfterWindowMin', afterWindowMinDefault, ...
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
corrTableSuffix = iP.Results.CorrTableSuffix;
corrMeasures = iP.Results.CorrMeasures;
directory = iP.Results.Directory;
relTimeWindowMin = iP.Results.RelativeTimeWindowMin;
beforeWindowMin = iP.Results.BeforeWindowMin;
afterWindowMin = iP.Results.AfterWindowMin;
avgStimDurationMin = iP.Results.StimDurationMin;
yLimits = iP.Results.YLimits;
yLimitsLog2Ratio = iP.Results.YLimitsLog2Ratio;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_raster() 
%                           or plot_chevron() or plot_psth() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output
handles = struct;

% Set default directory
if isempty(directory)
    directory = pwd;
end

% Set default figure suffix
if isempty(figSuffix)
    if isempty(stimIndices)
        figSuffixAddition = 'all_stims';
    elseif isnumeric(stimIndices)
        figSuffixAddition = ['stim', create_label_from_sequence(stimIndices)];
    elseif ischar(stimIndices)
        figSuffixAddition = ['all_', stimIndices, '_stims'];
    else
        error('stimNumber unrecognized!');
    end

    if ~isempty(beforeWindowMin)
        figSuffixAddition = sprintf('%s_bef%gto%g', figSuffixAddition, ...
                                    beforeWindowMin(1), beforeWindowMin(2));
    end

    if ~isempty(afterWindowMin)
        figSuffixAddition = sprintf('%s_aft%gto%g', figSuffixAddition, ...
                                    afterWindowMin(1), afterWindowMin(2));
    end

    figSuffix = [plotType, '_', figSuffixAddition];
end

% Set default figure name
if isempty(figName)
    % Extract the directory base
    dirBase = extract_fileparts(directory, 'dirbase');

    if ~isempty(relTimeWindowMin)
        % Create a sequence for create_label_from_sequence
        relTimeWindowSeq = linspace(relTimeWindowMin(1), ...
                                    relTimeWindowMin(2), 2);

        % Create a figure name
        figName = fullfile(directory, [dirBase, '_', ...
                            create_label_from_sequence(relTimeWindowSeq), ...
                            '_', figSuffix]);
    else
        % Create a figure name
        figName = fullfile(directory, [dirBase, '_', figSuffix]);
    end
end

% Set default figure title
if isempty(figTitle)
    if isempty(stimIndices)
        figTitle = 'SWD count around all stims';
    elseif isnumeric(stimIndices)
        figTitle = ['SWD count around stims ', ...
                    create_label_from_sequence(stimIndices)];
    elseif ischar(stimIndices)
        figTitle = ['SWD count around ', stimIndices, ' stims'];
    end
end

% Decide on default y limits
if isempty(yLimits)
    switch plotType
        case 'chevron'
            yLimits = [0, Inf];
        otherwise
            % Keep empty
    end
end

%% Get relative event times
% Load matching stimulus and event tables
[stimTables, swdTables, distinctParts] = ...
    read_matching_sheets(stimTableSuffix, eventTableSuffix, ...
                        'Directory', directory);

% Set default labels for each raster
if isempty(labels)
    labels = distinctParts;
end

% Extract all start times in seconds
[stimStartTimesSec, swdStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Extract stim durations in seconds
stimDurationsSec = cellfun(@(y) y.duration, stimTables, 'UniformOutput', false);

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
[stimStartTimesMin, stimDurationsMin, swdStartTimesMin] = ...
    argfun(@(x) cellfun(@(y) y / SEC_PER_MIN, x, 'UniformOutput', false), ...
            stimStartTimesSec, stimDurationsSec, swdStartTimesSec);

% Compute default stimulation duration in minutes
if isempty(avgStimDurationMin)
    % Find the minimum and maximum stimulation duration in minutes
    minStimDurationMin = apply_iteratively(@min, stimDurationsMin);
    maxStimDurationMin = apply_iteratively(@max, stimDurationsMin);

    % If they don't agree within 1%, plot stimulus duration as 0
    if (maxStimDurationMin - minStimDurationMin) / minStimDurationMin > 0.05
        fprintf(['Maximum stimulus duration %g and ' , ...
                    'minimum stimulus duration %g ', ...
                    'are more than 5%% apart, ', ...
                    'so stimulus duration will be plotted as 0!\n'], ...
                    maxStimDurationMin, minStimDurationMin);
        avgStimDurationMin = 0;
    else
        avgStimDurationMin = mean([maxStimDurationMin, minStimDurationMin]);
    end
end

%% Compute relative event times for each set
switch plotType
    case {'raster', 'chevron', 'corr'}
        % Compute the relative event times
        %   Note: this should return a cell matrix:
        %           Each column is a file
        %           Each row is a stim
        [relEventTimes, relTimeWindowMin] = ...
            compute_relative_event_times(swdStartTimesMin, stimStartTimesMin, ...
                                    'RelativeTimeWindow', relTimeWindowMin, ...
                                    'ForceMatrixOutput', true);
    case 'psth'
        % Relative event times computed in plot_psth.m
    otherwise
        error('plotType unrecognized!');
end

%% Plot event times
switch plotType
case 'raster'
    %% Plot rasters by stim #
    % Count the appropriate number of subplots
    nSubplots = size(relEventTimes, 1);

    % Count the number of files
    nFiles = size(relEventTimes, 2);

    if nSubplots == 1 && nFiles > 1
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
    
    % Create subplots
    [fig, ax] = create_subplots(1, nSubplots);

    % Plot the rasters
    for iAx = 1:numel(ax)
        % Use this subplot
        subplot(ax(iAx));

        % Create a figure title
        figTitle = ['Events around stim #', num2str(iAx)];

        % Plot raster
        rasters = plot_raster(relEventTimes(iAx, :), 'Labels', labels, ...
                                'XLabel', 'Time (min)', ...
                                'YLabel', 'File #', ...
                                'FigTitle', figTitle, ...
                                'XLimits', relTimeWindowMin, ...
                                'YLimits', yLimits, ...
                                'LegendLocation', legendLocation, ...
                                otherArguments);

        % Plot stim start line
        stimStartLine = plot_vertical_line(0, 'Color', stimStartLineColor, ...
                                        'LineWidth', stimStartLineWidth);
                                    
        % Plot stimulation window
        if ~isempty(avgStimDurationMin)
            stimWindow = [0, avgStimDurationMin];
            stimWindowShade = plot_vertical_shade(stimWindow);
        end
    end

    % Save figure
    save_all_figtypes(fig, figName, figTypes);

    % Return in handles
    handles.fig = fig;
    handles.ax = ax;
    handles.rasters = rasters;
    handles.stimStartLine = stimStartLine;
    handles.stimWindowShade = stimWindowShade;
case 'psth'
    %% Plot the peri-stimulus time histogram
    handles = plot_psth('EventTimes', swdStartTimesMin, ...
                        'StimTimes', stimStartTimesMin, ...
                        'GroupingLabels', labels, ...
                        'YLimits', yLimits, ...
                        'XLabel', 'Time (min)', ...
                        'RelativeTimeWindow', relTimeWindowMin, ...
                        'StimDuration', avgStimDurationMin, ...
                        'FigTitle', figTitle, ...
                        'FigName', figName, 'FigTypes', figTypes, ...
                        otherArguments);
case {'chevron', 'corr'}
    %% Plot two Chevron plots
    % Decide on p tick labels
    pTickLabels = {'Before', 'After'};

    % Modify the figure title for the log2 ratio plot
    figTitleRatio = replace(figTitle, 'SWD count', 'SWD count ratio');

    % Decide on file names
    figPathBase = extract_fileparts(figName, 'pathbase');
    sheetPath = [figPathBase, '.csv'];
    sheetPathLog2Ratio = [figPathBase, '_log2ratio', '.csv'];

    % Transpose so that each column is a stim
    relEventTimesTrans = transpose(relEventTimes);
    
    % Compute the number of events before and after, 
    %       summing across stims for each file
    %   Note: relEventTimes must be a cell array of numeric vectors
    [nEventsBeforeEachStim, nEventsAfterEachStim] = ...
        count_events(relEventTimesTrans, 'BeforeWindow', beforeWindowMin, ...
                    'AfterWindow', afterWindowMin);
    [nEventsBefore, nEventsAfter] = ...
        argfun(@(x) sum(x, 2), nEventsBeforeEachStim, nEventsAfterEachStim);

    switch plotType
        case 'chevron'
            % Compute log2 ratio data
            nEventsBeforeLog2Ratio = log2(nEventsBefore ./ nEventsBefore);
            nEventsAfterLog2Ratio = log2(nEventsAfter ./ nEventsBefore);

            % Save the data in tables
            chevronTable = table(nEventsBefore, nEventsAfter, ...
                                'RowNames', labels);
            writetable(chevronTable, sheetPath, 'WriteRowNames', true);

            % Save log2 ratio data in a table
            log2ratioChevronTable = table(nEventsBeforeLog2Ratio, ...
                                    nEventsAfterLog2Ratio, 'RowNames', labels);
            writetable(log2ratioChevronTable, sheetPathLog2Ratio, ...
                        'WriteRowNames', true);

            % TODO: Use plot_small_chevrons.m
            % Create subplots
            if plotLog2Ratio
                [fig, ax] = create_subplots(1, 2, 'FigExpansion', [1, 1]);
            else
                [fig, ax] = create_subplots(1, 1);
            end

            % Plot Chevron plot and save figure
            plot_chevron(chevronTable, 'FigTitle', figTitle, ...
                    'ReadoutLabel', 'SWD Count', 'PTickLabels', pTickLabels, ...
                    'ReadoutLimits', yLimits, 'LegendLocation', 'northeast', ...
                    'AxesHandle', ax(1), 'FigExpansion', [], otherArguments);

            % Plot log2 ratio Chevron plot and save figure
            if plotLog2Ratio
                plot_chevron(log2ratioChevronTable, ...
                    'FigTitle', figTitleRatio, ...
                    'IsLog2Data', true, ...
                    'ReadoutLabel', 'SWD Count Ratio', ...
                    'PTickLabels', pTickLabels, ...
                    'ReadoutLimits', yLimitsLog2Ratio, ...
                    'LegendLocation', 'northeast', ...
                    'AxesHandle', ax(2), 'FigExpansion', [], otherArguments);
            end

            % Save figure
            save_all_figtypes(fig, figName, figTypes);

            % Return in handles
            handles.fig = fig;
        case 'corr'
            % Find matching correlation table paths
            [~, corrPaths] = ...
                find_matching_files(distinctParts, 'Directory', directory, ...
                                    'Suffix', corrTableSuffix, ...
                                    'Extension', 'csv', 'ReturnEmpty', true);

            % Read the tables
            corrTables = cellfun(@readtable, corrPaths, 'UniformOutput', false);

            % Force corrMeasures as a column cell array
            corrMeasures = force_column_cell(corrMeasures);

            % Plot for each correlation measure
            handles = cellfun(@(x) plot_correlation_measure(timeStr, x, ...
                                    corrTables, stimStartTimesMin, ...
                                    beforeWindowMin, afterWindowMin, ...
                                    nEventsBefore, nEventsAfter, ...
                                    labels, figName), ...
                            corrMeasures);
        otherwise     
    end
otherwise
    error('plotType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = plot_correlation_measure (timeStr, corrMeasure, ...
                            corrTables, stimStartTimesMin, ...
                            beforeWindowMin, afterWindowMin, ...
                            nEventsBefore, nEventsAfter, labels, figNameOrig)
%% Plots the correlation of change in SWD count against change in a measure

%% Hard-coded parameters
figWidth = 3.4;
figHeight = 3;

% Replace in figure name
figName = replace(figNameOrig, 'corr', ['corr_', corrMeasure]);

% Extract path base
figPathBase = extract_fileparts(figName, 'pathbase');

% Extract the time instants from the tables
timeVecsSec = cellfun(@(x) extract_vars(x, timeStr), ...
                        corrTables, 'UniformOutput', false);

% Convert time to minutes
timeVecsMin = convert_units(timeVecsSec, 's', 'min');

% Extract the measure values from the tables
measureVecs = cellfun(@(x) extract_vars(x, corrMeasure), ...
                        corrTables, 'UniformOutput', false);

% Extract a time window for each stimulus time
[beforeWindows, afterWindows] = ...
    argfun(@(x) cellfun(@(y) create_shifted_vectors(x, y), ...
                        stimStartTimesMin, 'UniformOutput', false), ...
            beforeWindowMin, afterWindowMin);

% Compute the average of the measure values before and after stimuli
[avgMeasureBefore, avgMeasureAfter] = ...
    argfun(TODO , ...
            beforeWindows, afterWindows);

%% Compute percent changes
measurePercChange = 100 * (avgMeasureAfter - avgMeasureBefore) ./ ...
                        avgMeasureBefore;
nEventsPercChange = 100 * (nEventsAfter - nEventsBefore) ./ nEventsBefore;

%% Plot correlation
% Create figure
fig = set_figure_properties('AlwaysNew', true);

xLabel = sprintf('\\Delta %s (%%)', corrMeasure);
yLabel = '\Delta Event Count (%)';

% Plot scatter plot
handles = plot_grouped_scatter(measurePercChange, nEventsPercChange, labels, ...
                    'XLabel', xLabel, 'YLabel', yLabel, ...
                    'Markersize', 10, 'MarkerLineWidth', 3);

% Plot correlation coefficient
handles.corrText = plot_correlation_coefficient;

% Save figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Update for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'ScatterMarkerSize', 3, 'ScatterLineWidth', 1, ...
                        'RemoveLegends', true, 'RemoveText', true);

% Save figure
save_all_figtypes(fig, figPathBase, {'epsc', 'png'});

handles.fig = fig;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

relEventTimes = extract_in_order(relEventTimesCellCell);

% Extract the relevant event times
if stimIndices
    % Restrict to just the first event times
    %   Note: this should return a cell array of numeric vectors
    relEventTimes = cellfun(@(x) extract_element_by_index(x, 1, ...
                                                            false), ...
                        relEventTimesCellCell, 'UniformOutput', false);
else
    relEventTimes = extract_in_order(relEventTimesCellCell);
end

if stimIndices
    % Create a figure title
    figTitle = ['Events around stim #', num2str(1)];

    % Plot raster
    handles = plot_raster(relEventTimes, 'Labels', labels, ...
                            'XLabel', 'Time (min)', ...
                            'FigTitle', figTitle, ...
                            'XLimits', relTimeWindowMin, otherArguments);

    % Plot stim start line
    plot_vertical_line(0, 'LineWidth', 2, 'Color', 'k');
else
end

% Count the appropriate number of subplots
nSubplots = numel(relEventTimes);

% Force as column vectors
[nEventsBefore, nEventsAfter] = ...
    argfun(@force_column_vector, nEventsBefore, nEventsAfter);

% Generate the data for the Chevron plot
chevronData = transpose([nEventsBefore, nEventsAfter]);

if isempty(labels)
    labels = replace(distinctParts, '_', '\_');
end

stimStartTimesSec = cellfun(@(x) x(stimIndices), stimStartTimesSec, ...
                            'UniformOutput', false);
stimDurationsSec = cellfun(@(x) x(stimIndices), stimDurationsSec, ...
                            'UniformOutput', false);

% Extract the relevant event times
% TODO: May not be necessary
if stimIndices
    relEventTimes = relEventTimes(1, :);
end

figNameNormalized = strcat(figPathBase, '_normalized');

fig(1) = set_figure_properties('AlwaysNew', true);
save_all_figtypes(fig(1), figName, figTypes);
fig(2) = set_figure_properties('AlwaysNew', true);
save_all_figtypes(fig(2), figNameNormalized, figTypes);

% Check if event times exist
if isempty(relEventTimes)
    disp('There are no relative event times to plot!');
    return
end

figTitleNormalize = ['% ', figTitle];

% Decide on y limits
% TODO: Adjust to be integer?
if isempty(yLimitsLog2Ratio)
    yLimitsLog2Ratio = [-1, Inf];
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
