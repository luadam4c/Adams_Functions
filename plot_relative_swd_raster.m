function handles = plot_relative_swd_raster (varargin)
%% Plots a relative SWD raster from all gas_pulses.csv files and SWDs.csv files in a directory (unfinished)
% Usage: handles = plot_relative_swd_raster (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_relative_swd_raster('Directory', '/media/shareX/2019octoberR01/Figures/Figure1c')
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'FirstOnly': whether to take only the first window 
%                                   from each file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': directory to look for SWD table files
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
%                   - Any other parameter-value pair for plot_raster()
%
% Requires:
%       cd/compute_relative_event_times.m
%       cd/plot_raster.m
%
%   TODO: Update docs
%       cd/all_files.m
%       cd/create_label_from_sequence.m
%       cd/extract_distinct_fileparts.m
%       cd/extract_elements.m
%       cd/extract_fileparts.m
%
% Used by:
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-10 Created by Adam Lu
% TODO: Use load_matching_sheets.m
% TODO: Combine with plot_swd_raster.m
% 

%% Hard-coded parameters
SEC_PER_MIN = 60;

%   TODO: Create load_gas_pulses.m;
% Note: Must be consistent with parse_iox.m and parse_gas_trace.m
stimTableSuffix = '_gas_pulses';

% Note: Must be consistent with all_swd_sheets.m
swdTableSuffix = '_SWDs';

% TODO: Make optional arguments
pathBase = '';
sheetType = 'csv';
figSuffix = '';
labels = {};

%% Default values for optional arguments
firstOnlyDefault = false;       % take all windows by default
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
addParameter(iP, 'FirstOnly', firstOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
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
firstOnly = iP.Results.FirstOnly;
directory = iP.Results.Directory;
relTimeWindowMin = iP.Results.RelativeTimeWindowMin;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_raster() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default figure suffix
if isempty(figSuffix)
    if firstOnly
        figSuffix = 'first_windows';
    else
        figSuffix = 'all_windows';
    end
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
    if firstOnly
        figTitle = '# of SWDs around first stimulus starts';
    else
        figTitle = '# of SWDs around all stimulus starts';
    end
end

%% Do the job
% Get all stim pulse files
[~, stimPaths] = ...
    all_files('Directory', directory, 'Keyword', pathBase, ...
                'Suffix', stimTableSuffix, 'Extension', sheetType);

% Extract distinct prefixes
distinctPrefixes = extract_distinct_fileparts(stimPaths);

% Look for corresponding SWD spreadsheet files
[~, swdPaths] = ...
    cellfun(@(x) all_files('Prefix', x, 'Directory', directory, ...
                    'Suffix', swdTableSuffix, 'Extension', sheetType), ...
            distinctPrefixes, 'UniformOutput', false);

% Read all tables
[stimTables, swdTables] = ...
    argfun(@(x) cellfun(@readtable, x, 'UniformOutput', false), ...
            stimPaths, swdPaths);

% Set default labels for each raster
if isempty(labels)
    labels = strrep(distinctPrefixes, '_', '\_');
end

% Extract all start times in seconds
[stimStartTimesSec, swdStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Extract stim durations in seconds
[stimStartTimesSec, swdStartTimesSec] = ...
    argfun(@(x) cellfun(@(y) y.startTime, x, 'UniformOutput', false), ...
            stimTables, swdTables);

% Restrict to the first events only if requested
if firstOnly
    % Note: TODO: Use extract_elements.m with 'UniformOutput', false
    stimStartTimesSec = cellfun(@(x) x(1), stimStartTimesSec, ...
                                'UniformOutput', false);
end

% Convert to minutes
[stimStartTimesMin, swdStartTimesMin] = ...
    argfun(@(x) cellfun(@(y) y / SEC_PER_MIN, x, 'UniformOutput', false), ...
            stimStartTimesSec, swdStartTimesSec);

% Compute the relative event times
%   Note: this should return a cell array of cell arrays
relEventTimesCellCell = ...
    compute_relative_event_times(swdStartTimesMin, stimStartTimesMin, ...
                                    'RelativeTimeWindow', relTimeWindowMin);

% Restrict to just the first event times
if firstOnly
    %   Note: this should return a cell array of numeric vectors
    relEventTimes = cellfun(@extract_first_element, relEventTimesCellCell, ...
                            'UniformOutput', false);
else
    error('Not implemented yet!');
end

%% Plot the raster
if firstOnly
    handles = plot_raster(relEventTimes, 'Labels', labels, ...
                            'XLabel', 'Time (min)', ...
                            'XLimits', relTimeWindowMin, otherArguments);
    plot_vertical_line(0, 'LineWidth', 2, 'Color', 'k');
    save_all_figtypes(gcf, '/media/shareX/2019octoberR01/Figures/Figure1c/Figure1c', {'png', 'epsc2'})
else
    error('Not implemented yet!');
end

%% Plot a Chevron plot
% TODO: Move this to its own function
if firstOnly
    % Compute the number of events before and after
    %   TODO: Always cell arrays?
    nEventsBefore = cellfun(@(x) numel(x(x < 0)), relEventTimes);
    nEventsAfter = cellfun(@(x) numel(x(x >= 0)), relEventTimes);

    % Force as column vectors
    %   TODO: may not be necessary
    [nEventsBefore, nEventsAfter] = ...
        argfun(@force_column_vector, nEventsBefore, nEventsAfter);

    % Generate the data for the Chevron plot
    chevronData = transpose([nEventsBefore, nEventsAfter]);

    % TODO: plot_chevron.m
    [lowBefore, lowAfter] = ...
        argfun(@(x) compute_stats(x, 'lower95'), nEventsBefore, nEventsAfter);
    [highBefore, highAfter] = ...
        argfun(@(x) compute_stats(x, 'upper95'), nEventsBefore, nEventsAfter);

    [meanBefore, meanAfter] = ...
        argfun(@(x) compute_stats(x, 'mean'), nEventsBefore, nEventsAfter);
    pValue = [1, 2];

    % Plot a tuning curve
    figure1e = set_figure_properties('AlwaysNew', true, 'FigExpansion', [0.7, 0.7]);
    plot_tuning_curve(pValue, chevronData, 'PLimits', [0.5, 2.5], ...
                        'RunTTest', true, 'RunRankTest', true, ...
                        'Marker', 'o', 'MarkerFaceColor', [0, 0, 0], ...
                        'MarkerSize', 8, 'ColorMap', [0, 0, 0], ...
                        'FigHandle', figure1e, ...
                        'LegendLocation', 'suppress');
    hold on 
    plot(pValue, [meanBefore, meanAfter], 'r-o', ...
        'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot_error_bar(pValue, [lowBefore, lowAfter], [highBefore, highAfter], ...
        'Color', 'r', 'LineWidth', 2);
    save_all_figtypes(figure1e, '/media/shareX/2019octoberR01/Figures/Figure1e/Figure1e', {'png', 'epsc2'})
else
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function first = extract_first_element (vec)
%% Take the first element or return empty
% TODO: Merge with extract_elements

if numel(vec) >= 1
    if iscell(vec)
        first = vec{1};
    else
        first = vec(1);
    end
else
    first = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
