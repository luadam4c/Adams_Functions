function [bars, fig, outlines] = plot_grouped_histogram (varargin)
%% Plots a grouped histogram
% Usage: [bars, fig, outlines] = plot_grouped_histogram (stats (opt), grouping (opt), varargin)
% Explanation:
%       Plots a grouped histogram, placing bars from different groups
%           side-by-side by default
%       Notes: 
%           1. The bar() function is actually used for the main histogram
%           2. To compute just the histogram counts, 
%               use compute_grouped_histcounts.m
% Example(s):
%       randVec = randn(100, 1);
%       stats1 = [randVec, randVec + 1, randVec - 1];
%       plot_grouped_histogram(stats1)
%       plot_grouped_histogram(stats1, 'Style', 'stacked')
%       plot_grouped_histogram(stats1, 'Style', 'overlapped')
%       stats2 = [randVec; randVec + 1; randVec - 1];
%       grouping2 = [repmat({'Mark'}, 100, 1); repmat({'Peter'}, 100, 1); repmat({'Katie'}, 100, 1)];
%       plot_grouped_histogram(stats2, grouping2)
%       plot_grouped_histogram(stats2, grouping2, 'Style', 'stacked')
%       plot_grouped_histogram(stats2, grouping2, 'Style', 'overlapped')
% Outputs:
%       bars        - the histogram(s) returned as Bar object(s)
%                   specified as a Bar object handle array
%       fig         - figure handle for the created figure
%                   specified as a Figure object handle
%       outlines    - outlines of the histogram
%                   specified as a Stair object handle array
%
% Arguments:
%       stats       - (opt) data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       grouping    - (opt) group assignment for each data point
%                   must be an array of one the following types:
%                       'cell', 'string', numeric', 'logical', 
%                           'datetime', 'duration'
%                   default == the column number for a 2D array
%       varargin    - 'Counts': bin counts, with each group 
%                                   being a different column
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_grouped_histcounts(stats)
%                   - 'Edges': bin edges
%                   must be a vector of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%                   default == returned by compute_grouped_histcounts(stats)
%                   - 'Style': histogram style
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'side-by-side'  - same as 'grouped' bar graph
%                       'stacked'       - same as 'stacked' bar graph
%                       'overlapped'    - overlapped histograms
%                   default == 'side-by-side'
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == strcat('Value (', xUnits, ')')
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Count'
%                   - 'GroupingLabels': labels for the groupings, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == {'Group #1', 'Group #2', ...}
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nGroups == 1 
%                               'northeast' if nGroups is 2~9
%                               'eastoutside' if nGroups is 10+
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == strcat('Distribution of ', xLabel)
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'OutFolder': directory to save figure, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the bar() function
%
% Requires:
%       cd/compute_centers_from_edges.m
%       cd/compute_grouped_histcounts.m
%       cd/construct_fullpath.m
%       cd/create_default_grouping.m
%       cd/create_error_for_nargin.m
%       cd/islegendlocation.m
%       cd/ispositiveintegerscalar.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_histogram.m
%       cd/plot_swd_histogram.m
%       cd/ZG_fit_IEI_distributions.m

% TODO:
%       /media/adamX/Paula_IEIs/paula_iei4.m
%       /home/Matlab/Marks_Functions/paula/Oct2017/freqsPostJustinPartTwo.m
%
% 2017-12-11 Created by Adam Lu
% 2018-05-18 Added outFolder as a parameter
% 2018-05-25 Now doesn't plot if stats is empty
% 2018-12-27 Made all except one arguments optional and added many more
% 2018-12-28 Added usage of struct2arglist.m
% 2019-01-09 Now allows grouping to be a cell array of character vectors
% 2019-01-09 Added 'overlapped' as a style
% 2019-01-11 Improved default bin edges
% 2019-01-13 Now returns counts and edges as outputs
% 2019-01-15 Added 'Counts' as an optional parameter
% 2019-01-15 Moved code to compute_grouped_histcounts.m

%% Hard-coded parameters
validStyles = {'side-by-side', 'stacked', 'overlapped'};
barWidth = 1;                   % set this to 1 to make bars touch each other
maxInFigure = 8;                % maximum number of groups to keep the legend
                                %   inside the figure

%% Default values for optional arguments
statsDefault = [];              % set later
groupingDefault = [];           % set later
countsDefault = [];             % set later
edgesDefault = [];              % set later
styleDefault = 'side-by-side';  % plot bars side-by-side by default
xLimitsDefault = 'suppress';    % set later
yLimitsDefault = 'suppress';    % set later
xUnitsDefault = 'unit';         % the default x-axis units
xLabelDefault = '';             % set later
yLabelDefault = 'Count';        % set later
groupingLabelsDefault = '';     % set later
legendLocationDefault = 'auto'; % set later
figTitleDefault = '';           % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
outFolderDefault = '';          % default directory to save figure
figNameDefault = '';            % don't save figure by default
figTypesDefault = 'png';        % save as png file by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'stats', statsDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addOptional(iP, 'grouping', groupingDefault, ...
    @(x) validateattributes(x, {'cell', 'string', 'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Counts', countsDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Edges', edgesDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'Style', styleDefault, ...
    @(x) any(validatestring(x, validStyles)));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XUnits', xUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'GroupingLabels', groupingLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
stats  = iP.Results.stats;
grouping = iP.Results.grouping;
counts = iP.Results.Counts;
edges = iP.Results.Edges;
style = validatestring(iP.Results.Style, validStyles);
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xUnits = iP.Results.XUnits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
groupingLabels = iP.Results.GroupingLabels;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the bar() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Return if there is nothing to plot
if isempty(stats) && (isempty(counts) || isempty(edges))
    fprintf('Both counts and edges must be provided if stats is empty!\n');
    bars = gobjects(0);
    fig = gobjects(0);
    return
end

% Return warning
if ~isempty(stats) && ~isempty(counts)
    fprintf('Warning: Custom counts are provided, stats will be ignored!\n');
end

% Decide on the grouping vector and possibly labels
[grouping, groupingLabels] = ...
    create_default_grouping('Stats', stats, 'Counts', counts, ...
                        'Grouping', grouping, 'GroupingLabels', groupingLabels);

% Get all unique group values
groupValues = unique(grouping);

% Count the number of groups
nGroups = numel(groupValues);

% Compute the bin counts and bin edges
if isempty(counts)
    [counts, edges, binCenters] = ...
        compute_grouped_histcounts(stats, grouping, 'Edges', edges);
else
    binCenters = compute_centers_from_edges(edges);
end

%% Preparation for the plot
% Set the bar graph style
switch style
    case {'side-by-side', 'overlapped'}
        barStyle = 'grouped';
    case 'stacked'
        barStyle = 'stacked';
    otherwise
        error('style is unrecognized!!');
end

% Set the default x-axis labels
if isempty(xLabel)
    xLabel = strcat('Value (', xUnits, ')');
end

% Set the default figure title
if isempty(figTitle)
    figTitle = ['Distribution of ', xLabel];
end

% If the figure name is not a full path, create full path
if ~isempty(figName)
    figName = construct_fullpath(figName, 'Directory', outFolder);
end

% Decide on the figure to plot on
if ~isempty(figHandle)
    fig = figure(figHandle);
elseif ~isempty(figNumber)
    fig = figure(figNumber);
elseif ~isempty(figName)
    fig = figure;
else
    fig = gcf;
end

% Set legend location based on number of groups
if strcmpi(legendLocation, 'auto')
    if nGroups > 1 && nGroups <= maxInFigure
        legendLocation = 'northeast';
    elseif nGroups > maxInFigure
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

% Compute appropriate edge and face alpha
if strcmpi(style, 'overlapped')
    faceAlpha = 1/nGroups;
    edgeAlpha = 1/nGroups;
end

%% Plot and save histogram
if strcmpi(style, 'overlapped')
    % Store hold status
    wasHold = ishold;

    % Hold on
    hold on;

    % Plot bars
    bars = arrayfun(@(x, y) bar(binCenters, counts(:, x), ...
                        barWidth, barStyle, 'FaceAlpha', faceAlpha, ...
                        'EdgeAlpha', edgeAlpha, otherArguments{:}), ...
                    transpose(1:nGroups), groupingLabels);
else
    % Plot bars
    bars = bar(binCenters, counts, barWidth, barStyle, otherArguments{:});
end

% Count the number of groups
nGroups = numel(bars);

% Set the legend labels for each Bar object
for iBar = 1:nGroups
    set(bars(iBar), 'DisplayName', groupingLabels{iBar});
end

% Set x axis limits
if ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~isempty(yLimits) && ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Generate a legend if there is more than one group
if ~strcmpi(legendLocation, 'suppress')
    legend('location', legendLocation, 'AutoUpdate', 'off');
end

% Create staircase outlines if style is 'overlapped'
if strcmpi(style, 'overlapped')
    % Retrieve the colors used for the histograms
    colorsUsed = arrayfun(@(x) get(x, 'FaceColor'), bars, 'UniformOutput', false);
    colorsUsed = vertcat(colorsUsed{:});

    % Create outlines
    outlines = arrayfun(@(x, y) stairs(edges, [counts(:, x); 0], ...
                                'LineWidth', 2, 'Color', colorsUsed(x, :)), ...
                    transpose(1:nGroups), groupingLabels);
else
    outlines = gobjects(nGroups, 1);
end

% Generate an x-axis label
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end

% Generate a y-axis label
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end

% Generate a title
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

% Hold off if it was originally so
if ~wasHold
    hold off
end

% Save the figure
if ~isempty(figName)
    % Save the figure in all file types requested
    save_all_figtypes(fig, figName, figTypes);

    % Close figure
    close(fig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

histg(stats, grouping);

bars = figure('Visible', 'on');

saveas(bars, figName, 'png');
close(bars);

plot_grouped_histogram(figName, stats, grouping, groupingLabels, ...
                        xLabel, xUnits, figTitle, varargin)

addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric', 'categorical', ...
        'datetime', 'duration'}, {'vector', 'numel', 2}));

if ~isempty(xUnits)
    xlabel([xLabel, ' (', xUnits, ')']);
end

% Does not work for datetime arrays
binCenters = (edges(1:end-1) + edges(2:end)) / 2;

groupingLabels = create_labels_from_numbers(groupValues, ...
                                            'Prefix', 'Group #');

legend(groupingLabels, 'location', legendLocation, ...
        'Interpreter', 'none', 'AutoUpdate','off');

title(figTitle, 'Interpreter', 'none');

bars = bar(binCenters, counts, barWidth, barStyle);

% The following is slower:
counts = zeros(nBins, nGroups);
parfor iGroup = 1:nGroups
    % Get the value of this group
    groupValueThis = groupValues(iGroup);

    % Extract the stats for this group
    statsThisGroup = stats(grouping == groupValueThis);

    % Compute the bin counts for this group
    counts(:, iGroup) = compute_bins(statsThisGroup, 'Edges', edges);
end

% Count the number of bins
nBins = numel(edges) - 1;

if ~isempty(stats)
end

[~, edges] = compute_bins(stats, 'Edges', edges);

legend(groupingLabels, 'location', legendLocation, 'AutoUpdate','off');

% If the figure name is not a full path, create full path
if ~isempty(figName) && ~any(strfind(figName, filesep))
    % Set dependent argument defaults
    if isempty(outFolder)
        % Default output directory is present working directory
        outFolder = pwd;
    end

    % Create full figure file name
    figName = fullfile(outFolder, figName);
end

% Clear figure
clf(fig);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
