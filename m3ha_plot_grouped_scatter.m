function handles = m3ha_plot_grouped_scatter (statsPath, varargin)
%% Plots grouped scatter plots from a statistics table returned by m3ha_compute_statistics.m
% Usage: handles = m3ha_plot_grouped_scatter (statsPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - structures with fields:
%                       fig
%                       scatters
%                   specified as a structure array
%
% Arguments:
%       statsPath  - path to a .mat file containing the variables:
%                       statsTable - statistics table 
%                                       returned by m3ha_compute_statistics.m
%                           must have these columns:
%                               measureTitle
%                               measureStr
%                               allValues
%                       pharmLabels
%                       conditionLabel
%                   must be a string scalar or a character vector
%       varargin    - 'RowsToPlot': rows to plot
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'OutFolder': the directory where plots will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigWidth': figure width
%                   must be a positive scalar
%                   default == 5 cm
%                   - 'FigHeight': figure height
%                   must be a positive scalar
%                   default == 3.4 cm
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for plot_grouped_scatter.m
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/convert_units.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/isfigtype.m
%       cd/ispositiveintegervector.m
%       cd/m3ha_decide_on_ylimits.m
%       cd/plot_correlation_coefficient.m
%       cd/plot_grouped_scatter.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m
%
% Used by:
%       cd/m3ha_plot_figure08.m

% File History:
% 2020-08-02 Modified from m3ha_plot_violin.m

%% Hard-coded parameters

%% Default values for optional arguments
rowsToPlotDefault = 'all';
outFolderDefault = '';          % set later
figWidthDefault = 3.4;
figHeightDefault = 3;
figTypesDefault = {'png', 'epsc'};

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

% Add required inputs to the Input Parser
addRequired(iP, 'statsPath', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RowsToPlot', rowsToPlotDefault, ...
    @(x) assert(ispositiveintegervector(x) || ischar(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['RowsToPlot must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigWidth', figWidthDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigWidth must be a empty or a positive scalar!'));
addParameter(iP, 'FigHeight', figHeightDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigHeight must be a empty or a positive scalar!'));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, statsPath, varargin{:});
rowsToPlot = iP.Results.RowsToPlot;
outFolder = iP.Results.OutFolder;
figWidth = iP.Results.FigWidth;
figHeight = iP.Results.FigHeight;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_grouped_scatter() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output
handles = struct;

% Set default output directory
if isempty(outFolder)
    outFolder = extract_fileparts(statsPath, 'directory');
end

% Load stats table
disp('Loading statistics for grouped scatter plots ...');
if isfile(statsPath)
    load(statsPath, 'statsTable', 'pharmLabels', 'conditionLabel');
else
    fprintf('%s does not exist!\n', statsPath);
    return;
end

% Restrict to measures to plot
if ~(ischar(rowsToPlot) && strcmp(rowsToPlot, 'all'))
    statsTable = statsTable(rowsToPlot, :);
end

% Extract variables
allMeasureTitles = statsTable.measureTitle;
allMeasureStrs = statsTable.measureStr;
allValues = statsTable.allValues;

% Count the number of variables
nMeasures = numel(allMeasureStrs);

% Create all combinations
allIndexPairs = force_column_cell(transpose(nchoosek(1:nMeasures, 2)));

% Create all pairs of variables
[allTitlePairs, allMeasurePairs, allValuePairs] = ...
    argfun(@(a) cellfun(@(x) {a{x(1)}; a{x(2)}}, ...
                        allIndexPairs, 'UniformOutput', false), ...
            allMeasureTitles, allMeasureStrs, allValues);

%% Do the job
% Plot all grouped scatter plots
disp('Plotting grouped scatter plots ...');
handles = ...
    cellfun(@(a, b, c) m3ha_plot_grouped_scatter_helper(a, b, c, ...
                            pharmLabels, figHeight, figWidth, figTypes, ...
                            outFolder, conditionLabel, otherArguments), ...
            allValuePairs, allMeasurePairs, allTitlePairs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = ...
            m3ha_plot_grouped_scatter_helper (valuePair, strPair, titlePair, ...
                            pharmLabels, figHeight, figWidth, figTypes, ...
                            outFolder, conditionLabel, otherArguments)

%% Hard-coded parameters

%% Preparation
% Create figure base
figBase = strcat(strjoin(fliplr(strPair), '_vs_'), '_', conditionLabel);

% Create full path base
figPathBase = fullfile(outFolder, figBase);

% Decide on the color map
colorMap = decide_on_colormap([], numel(pharmLabels));

%% Plot
% Create figure
fig = set_figure_properties('AlwaysNew', true);

% Convert onset times from ms to seconds
for i = 1:numel(titlePair)
    if contains(titlePair{i}, 'onset')
        % Update values
        valuePair{i} = convert_units(valuePair{i}, 'ms', 's');

        % Update title
        titlePair{i} = replace(titlePair{i}, 'ms', 's');
    end
end

% Plot groups as a grouped scatter plot
scatters = plot_grouped_scatter(valuePair{1}, valuePair{2}, ...
                'GroupingLabels', pharmLabels, 'ColorMap', colorMap, ...
                'XLabel', titlePair{1}, 'YLabel', titlePair{2}, ...
                otherArguments);

% Plot the correlation coefficient
plot_correlation_coefficient;

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Set x axis limits
xLimits = m3ha_decide_on_ylimits(titlePair{1});
if ~isempty(xLimits)
    xlim(xLimits);
end

% Set y axis limits
yLimits = m3ha_decide_on_ylimits(titlePair{2});
if ~isempty(yLimits)
    ylim(yLimits);
end

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'PlotMarkerSize', 3, 'RemoveLegend', true, ...
                        'RemoveText', true);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

% Save in handles
handles.fig = fig;
handles.scatters = scatters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
