function handles = m3ha_plot_bar3 (statsPath, varargin)
%% Plots 3-dimensional bar plots from a statistics table returned by m3ha_compute_statistics.m
% Usage: handles = m3ha_plot_bar3 (statsPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - structures with fields:
%                       fig
%                       bars
%                       errorBarVert
%                       errorBarHorz
%                   specified as a structure array
%
% Arguments:
%       statsPath  - path to a .mat file containing the variables:
%                       statsTable - statistics table 
%                                       returned by m3ha_compute_statistics.m
%                       pharmLabels
%                       gIncrLabels
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
%                   default == 6 cm
%                   - 'FigHeight': figure height
%                   must be a positive scalar
%                   default == 6 cm
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for bar3()
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/extract_fileparts.m
%       cd/isfigtype.m
%       cd/ispositiveintegervector.m
%       cd/struct2arglist.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m
%
% Used by:
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-12-27 Moved from m3ha_plot_figure02.m
% 2019-12-28 Added 'FigTypes' as an optional argument
% 2019-12-28 Added 'FigHeight' and 'FigWidth' as optional arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
rowsToPlotDefault = 'all';
outFolderDefault = '';          % set later
figWidthDefault = 6;
figHeightDefault = 6;
figTypesDefault = {'png', 'epsc2'};

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

% Keep unmatched arguments for the bar3() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Initialize output
handles = struct;

% Set default output directory
if isempty(outFolder)
    outFolder = extract_fileparts(statsPath, 'directory');
end

% Load stats table
disp('Loading statistics for 3D bar plots ...');
if isfile(statsPath)
    load(statsPath, 'statsTable', 'pharmLabels', ...
        'gIncrLabels', 'conditionLabel');
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
allMeanValues = statsTable.meanValue;
allUpper95Values = statsTable.upper95Value;

% Create figure bases
allFigBases = combine_strings({allMeasureStrs, conditionLabel});

% Create full path bases
allFigPathBases = fullfile(outFolder, allFigBases);

%% Do the job
% Plot all 3D bar plots
disp('Plotting 3D bar plots ...');
handles = cellfun(@(a, b, c, d) m3ha_plot_bar3_helper(a, b, c, ...
                                pharmLabels, gIncrLabels, ...
                                d, figHeight, figWidth, ...
                                figTypes, otherArguments), ...
                allMeanValues, allUpper95Values, ...
                allMeasureTitles, allFigPathBases);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_bar3_helper(meanValue, upper95Value, ...
                        measureTitle, pharmLabels, gIncrLabels, ...
                        figPathBase, figHeight, figWidth, figTypes, ...
                        otherArguments)

% Create figure for conductance traces
fig = set_figure_properties('AlwaysNew', true);

% Flip the g incr axis
[meanValue, upper95Value, gIncrLabels] = ...
    argfun(@fliplr, meanValue, upper95Value, gIncrLabels);

% Set x and y tick labels
xTickLabels = pharmLabels;
yTickLabels = gIncrLabels;

% TODO: Add the following to plot_bar.m?
% Hard-coded parameters
relativeBarWidth = 0.2;
xTickAngle = 320;
barSeparation = 1;

% Decide on the color map
cm = decide_on_colormap([], 4);

% Set the color map
colormap(cm);

% Prepare for bar3
meanValueTransposed = transpose(meanValue);
upper95ValueTransposed = transpose(upper95Value);

% Plot the means as bars
bars = bar3(meanValueTransposed, relativeBarWidth, 'detached', ...
            otherArguments{:});

% Plot error bars
% TODO: Incorporate into plot_error_bar.m?

% Set the relative error bar width to be the same as the bars themselves
%   Note: error bar width must not exceed the bar width, 
%           otherwise the edges would be cut off
relativeErrorBarWidth = relativeBarWidth;

% Compute the actual error bar width
errorBarWidth = relativeErrorBarWidth * barSeparation;

% Compute the x and y values corresponding to each data point
[xValues, yValues] = meshgrid(1:numel(xTickLabels), 1:numel(yTickLabels));

% Compute the left and right positions of the horizontal parts of the error bars
xPosBarLeft = xValues - errorBarWidth / 2;
xPosBarRight = xValues + errorBarWidth / 2;

% Plot the vertical part of the error bars
errorBarVert = ...
    arrayfun(@(a, b, c, d, e, f) line([a, b], [c, d], [e, f], 'Color', 'k'), ...
            xValues, xValues, yValues, yValues, ...
            meanValueTransposed, upper95ValueTransposed);

% Plot the horizontal part of the error bars
errorBarHorz = ...
    arrayfun(@(a, b, c, d, e, f) line([a, b], [c, d], [e, f], 'Color', 'k'), ...
            xPosBarLeft, xPosBarRight, yValues, yValues, ...
            upper95ValueTransposed, upper95ValueTransposed);

% Plot z axis label
zlabel(measureTitle);

% Set x tick labels
set(gca, 'XTickLabel', xTickLabels);

% Set x tick angle
xtickangle(xTickAngle);

% Set y tick labels
set(gca, 'YTickLabel', yTickLabels);

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

% Save in handles
handles.fig = fig;
handles.bars = bars;
handles.errorBarVert = errorBarVert;
handles.errorBarHorz = errorBarHorz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%