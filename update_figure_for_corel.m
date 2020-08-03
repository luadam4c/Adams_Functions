function figHandle = update_figure_for_corel (varargin)
%% Update figure to be journal-friendly (ready for CorelDraw)
% Usage: figHandle = update_figure_for_corel (figHandle (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = update_figure_for_corel(fig);
%
% Outputs:
%       figHandle   - handle to updated figure
%                   specified as a Figure object handle
%
% Arguments:
%       figHandle   - (opt) figure handle to update
%                   must be a Figure object handle
%       varargin    - 'AlignSubplots': whether to align subplots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'BoxOn': whether to add axes outlines
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveTicks': whether to remove all ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveXTicks': whether to remove all x ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveYTicks': whether to remove all y ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveRulers': whether to remove all rulers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveXRulers': whether to remove all x rulers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveYRulers': whether to remove all y rulers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveLabels': whether to remove x and y labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveXLabels': whether to remove x labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveYLabels': whether to remove y labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveTitles': whether to remove titles
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveLegends': whether to remove all legends
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveTexts': whether to remove all texts
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveCircles': whether to remove plots that are circles
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'XTickLocs': locations of X ticks
%                   must be 'suppress' or a numeric vector
%                   default == 'suppress'
%                   - 'YTickLocs': locations of Y ticks
%                   must be 'suppress' or a numeric vector
%                   default == 'suppress'
%                   - 'LabelsFontSize': labels and titles font size in points
%                   must be a positive integer scalar
%                   default == 8
%                   - 'AxisFontSize': axis font size in points
%                   must be a positive integer scalar
%                   default == 7
%                   - 'TextFontSize': text font size in points
%                   must be a positive integer scalar
%                   default == 7
%                   - 'RulerLineWidth': axis ruler line width in points
%                   must be a positive integer scalar
%                   default == 1
%                   - 'PlotLineWidth': line width of Line objects that 
%                                       have more than two data points
%                   must be empty or a positive scalar
%                   default == [] (don't change)
%                   - 'PlotMarkerSize': marker size of Line objects that 
%                                       have more than two data points
%                   must be empty or a positive scalar
%                   default == [] (don't change)
%                   - 'ScatterMarkerSize': marker size of Scatter objects
%                   must be empty or a positive scalar
%                   default == [] (don't change)
%                   - Any other parameter-value pair for set_figure_properties()
%
% Requires:
%       cd/align_subplots.m
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/extract_fields.m
%       cd/force_column_cell.m
%       cd/is_field.m
%       cd/match_positions.m
%       cd/remove_non_axes.m
%       cd/set_figure_properties.m
%       cd/set_visible_off.m
%
% Used by:
%       cd/m3ha_compute_and_plot_IT.m
%       cd/m3ha_oscillations_analyze.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_plot_figure07.m
%       cd/m3ha_plot_grouped_scatter.m
%       cd/m3ha_plot_violin.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_measures.m
%       cd/plot_tuning_curve.m
%       cd/plot_traces_spike2_mat.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-19 Created by Adam Lu
% 2019-09-20 Added 'RemoveTicks' as an optional argument
% 2019-10-04 Added 'RemoveLegends' as an optional argument
% 2019-10-05 Add textFontSize
% 2019-10-06 Changed default labelsFontSize from 8 to 7
% 2019-11-30 Changed default labelsFontSize from 7 to 8
% 2019-11-30 Changed default axisFontSize from 6 to 7
% 2019-11-30 Changed default textFontSize from 6 to 7
% 2019-12-01 Added 'RemoveXLabels' as an optional argument
% 2019-12-02 Fixed bug when there are multiple labels of the same type
% 2019-12-04 Added 'RemoveTexts' as an optional argument
% 2019-12-29 Added 'AlignSubplots' as an optional argument
% 2020-02-06 Added 'BoxOn' as an optional argument
% 2019-04-22 Added 'RemoveCircles' as an optional argument

%% Hard-coded parameters
BLACK = [0, 0, 0];

% TODO: Make optional parameters
units = 'inches';
tickLengthsInches = [0.025, 0.025];
annotationLineWidth = 1; % TODO

%% Default values for optional arguments
figHandleDefault = [];
alignSubplotsDefault = false;   % don't align by default
boxOnDefault = false;           % no box by default
removeTicksDefault = false;     % don't remove by default
removeXTicksDefault = false;    % don't remove by default
removeYTicksDefault = false;    % don't remove by default
removeRulersDefault = false;    % don't remove by default
removeXRulersDefault = false;   % don't remove by default
removeYRulersDefault = false;   % don't remove by default
removeLabelsDefault = false;    % don't remove by default
removeXLabelsDefault = false;   % don't remove by default
removeYLabelsDefault = false;   % don't remove by default
removeTitlesDefault = false;    % don't remove by default
removeLegendsDefault = false;   % don't remove by default
removeTextsDefault = false;     % don't remove by default
removeCirclesDefault = false;   % don't remove by default
xTickLocsDefault = 'suppress';  % don't change by default
yTickLocsDefault = 'suppress';  % don't change by default
labelsFontSizeDefault = 8;
axisFontSizeDefault = 7;
textFontSizeDefault = 7;
rulerLineWidthDefault = 1;
plotLineWidthDefault = [];
plotMarkerSizeDefault = [];
scatterMarkerSizeDefault = [];

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
addOptional(iP, 'figHandle', figHandleDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AlignSubplots', alignSubplotsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BoxOn', boxOnDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveTicks', removeTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveXTicks', removeXTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveYTicks', removeYTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveRulers', removeRulersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveXRulers', removeXRulersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveYRulers', removeYRulersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveLabels', removeLabelsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveXLabels', removeXLabelsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveYLabels', removeYLabelsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveTitles', removeTitlesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveLegends', removeLegendsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveTexts', removeTextsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveCircles', removeCirclesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XTickLocs', xTickLocsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || ...
                    isnumericvector(x) || iscellnumericvector(x), ...
        ['XTickLocs must be ''suppress'' or a numeric vector', ...
            'or a cell array of numeric vectors!']));
addParameter(iP, 'YTickLocs', yTickLocsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || ...
                    isnumericvector(x) || iscellnumericvector(x), ...
        ['YTickLocs must be ''suppress'' or a numeric vector', ...
            'or a cell array of numeric vectors!']));
addParameter(iP, 'LabelsFontSize', labelsFontSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'AxisFontSize', axisFontSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TextFontSize', textFontSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'RulerLineWidth', rulerLineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PlotLineWidth', plotLineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'PlotMarkerSize', plotMarkerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ScatterMarkerSize', scatterMarkerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
figHandle = iP.Results.figHandle;
alignSubplots = iP.Results.AlignSubplots;
boxOn = iP.Results.BoxOn;
removeTicks = iP.Results.RemoveTicks;
removeXTicks = iP.Results.RemoveXTicks;
removeYTicks = iP.Results.RemoveYTicks;
removeRulers = iP.Results.RemoveRulers;
removeXRulers = iP.Results.RemoveXRulers;
removeYRulers = iP.Results.RemoveYRulers;
removeLabels = iP.Results.RemoveLabels;
removeXLabels = iP.Results.RemoveXLabels;
removeYLabels = iP.Results.RemoveYLabels;
removeTitles = iP.Results.RemoveTitles;
removeLegends = iP.Results.RemoveLegends;
removeTexts = iP.Results.RemoveTexts;
removeCircles = iP.Results.RemoveCircles;
xTickLocs = iP.Results.XTickLocs;
yTickLocs = iP.Results.YTickLocs;
labelsFontSize = iP.Results.LabelsFontSize;
axisFontSize = iP.Results.AxisFontSize;
textFontSize = iP.Results.TextFontSize;
rulerLineWidth = iP.Results.RulerLineWidth;
plotLineWidth = iP.Results.PlotLineWidth;
plotMarkerSize = iP.Results.PlotMarkerSize;
scatterMarkerSize = iP.Results.ScatterMarkerSize;

% Keep unmatched arguments for set_figure_properties()
otherArguments = iP.Unmatched;

%% Preparation
% Compute font size multipliers
titleFontSizeMultiplier = labelsFontSize / axisFontSize;
labelFontSizeMultiplier = labelsFontSize / axisFontSize;

% Check if the figure handle is valid
if ~isempty(figHandle) && ~isvalid(figHandle)
    error('figHandle is not valid!');
end

% Update flags
if removeTicks
    removeXTicks = true;
    removeYTicks = true;
end
if removeRulers
    removeXRulers = true;
    removeYRulers = true;
end
if removeLabels
    removeXLabels = true;
    removeYLabels = true;
end

%% Set figure properties
% Might change sizes
%   Note: Changing the renderer to 'painters' ensure graphics 
%           are saved as vectors
figHandle = set_figure_properties('FigHandle', figHandle, ...
                        'Renderer', 'painters', otherArguments);

% Update figure position units
unitsOrig = get(figHandle, 'Units');
if ~strcmp(unitsOrig, units)
    set(figHandle, 'Units', units);
end

%% Set axes properties
% Find all subplots in the figure
ax = findall(figHandle, 'type', 'axes');

% Remove axes for suplabel or suptitle
ax = remove_non_axes(ax);

% Sort the subplots in the figure
ax = sort_subplots(ax);

% Count the number of axes
nAx = numel(ax);

% Make axes the same size if requested
if alignSubplots
    align_subplots(figHandle);
end

% Remove boxes
for iAx = 1:nAx
    box(ax(iAx), 'off');
end

% Make all ticks go outward
set(ax, 'TickDir', 'out', 'TickDirMode', 'manual');

% Deal with x axis ticks
if removeXTicks
    set(ax, 'XTick', []);
else
    % Change the x tick values
    if ~ischar(xTickLocs) || ~strcmpi(xTickLocs, 'suppress')
        if iscell(xTickLocs)
            xTickLocs = force_column_cell(xTickLocs);
            cellfun(@(a, b) update_ticks(a, 'x', b), num2cell(ax), xTickLocs);
        else
            arrayfun(@(a) update_ticks(a, 'x', xTickLocs), ax);
        end
    end
end

% Deal with y axis ticks
if removeYTicks
    set(ax, 'YTick', []);
else
    % Change the y tick values
    if ~ischar(yTickLocs) || ~strcmpi(yTickLocs, 'suppress')
        if iscell(yTickLocs)
            yTickLocs = force_column_cell(yTickLocs);
            cellfun(@(a, b) update_ticks(a, 'y', b), num2cell(ax), yTickLocs);
        else
            arrayfun(@(a) update_ticks(a, 'y', yTickLocs), ax);
        end
    end
end

% Remove x axis if requested
if removeXRulers
    xAxises = get(ax, 'XAxis');
    set_visible_off(xAxises);
end

% Remove y axis if requested
if removeYRulers
    yAxises = get(ax, 'YAxis');
    set_visible_off(yAxises);
end

% Remove x labels if requested
if removeXLabels
    xLabels = get(ax, 'XLabel');
    set_string_empty(xLabels);
end

% Remove y labels if requested
if removeYLabels
    yLabels = get(ax, 'YLabel');
    set_string_empty(yLabels);
end

% Remove titles if requested
if removeTitles
    titles = get(ax, 'Title');
    set_string_empty(titles);
end

% Remove legends if requested
if removeLegends
    lgds = findobj(gcf, 'Type', 'Legend');
    delete(lgds);
end

% Remove texts if requested
if removeTexts
    texts = findobj(gcf, 'Type', 'Text');
    delete(texts);
end

% Add axes outline if requested
if boxOn
    arrayfun(@(x) box(x, 'On'), ax);
end

% Set font
set(ax, 'FontName', 'Arial');
set(ax, 'FontSize', axisFontSize);
set(ax, 'TitleFontSizeMultiplier', titleFontSizeMultiplier);
set(ax, 'TitleFontWeight', 'normal');
set(ax, 'LabelFontSizeMultiplier', labelFontSizeMultiplier);
set(ax, 'XColor', BLACK, 'YColor', BLACK, 'ZColor', BLACK);

% Change the fontsize of texts
texts = findobj(figHandle, 'Type', 'Text');
if ~isempty(texts)
    set(texts, 'Fontsize', textFontSize);
end

% Set ruler linewidths
set(ax, 'LineWidth', rulerLineWidth);

% Update plot line widths
if ~isempty(plotLineWidth)
    lines = findobj(figHandle, 'Type', 'Line');
    plots = lines(arrayfun(@(x) is_plot(x), lines));
    set(plots, 'LineWidth', plotLineWidth);
end

% Remove markers
if removeCircles
    lines = findobj(figHandle, 'Type', 'Line');
    markers = extract_fields(lines, 'Marker', 'UniformOutput', false);
    lineStyles = extract_fields(lines, 'LineStyle', 'UniformOutput', false);
    toRemove = strcmp(markers, 'o') & strcmp(lineStyles, 'none');
    delete(lines(toRemove))
end

% Update marker sizes
if ~isempty(plotMarkerSize)
    lines = findobj(figHandle, 'Type', 'Line');
    plots = lines(arrayfun(@(x) is_plot(x), lines));
    set(plots, 'MarkerSize', plotMarkerSize);
end

% Update Scatter object marker sizes
if ~isempty(scatterMarkerSize)
    scatters = findobj(figHandle, 'Type', 'Scatter');
    set(scatters, 'SizeData', scatterMarkerSize^2);
end

% Update all face alphas to 1
%   Note: If not 1, CorelDraw can't import the figure as vector graphics
dots = findobj(figHandle, 'Type', 'Scatter');
set(dots, 'MarkerFaceAlpha', 1);
patches = findobj(figHandle, 'Type', 'Patch');
set(patches, 'FaceAlpha', 1, 'EdgeAlpha', 1);

% Update annotation line widths
% TODO: How to distinguish?

% Set tick lengths
figPosition = get(figHandle, 'Position');
for iAx = 1:nAx
    % Get the length of the longest axis
    axPosition = get(ax(iAx), 'Position');
    axisLengthUnits = max(axPosition(3:4) .* figPosition(3:4));

    % Convert to units relative to longest axis
    tickLengthsRel = tickLengthsInches / axisLengthUnits;

    % Set new tick lengths
    set(ax(iAx), 'TickLength', tickLengthsRel);
end

%% Restore things
% Restore figure position units
unitsNow = get(figHandle, 'Units');
if ~strcmp(unitsNow, unitsOrig)
    set(figHandle, 'Units', unitsOrig);
end

%% Change color coding
% Find all objects in the figure
objects = findall(figHandle);
arrayfun(@convert_color_to_rgb, objects);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isPlot = is_plot(lineObject)
% TODO: Pull out as its own function

% Get x, y and z data
xData = lineObject.XData;
yData = lineObject.YData;
zData = lineObject.ZData;

% The line is a plot if there are more than two data points
isPlot = numel(xData) > 2 || numel(yData) > 2 || numel(zData) > 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_string_empty (textObject)
% TODO: Pull out as its own function

if iscell(textObject)
    cellfun(@set_string_empty, textObject);
else
    textObject.String = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_ticks (ax, axisType, tickLocs)
%% Updated the tick locations for an axes
%   Note: this function won't remove ticks if tickLocs is empty
% TODO: Pull out as its own function

if isempty(tickLocs)
    return
end

% Decide on which ticks to update
switch axisType
    case 'x'
        tickStr = 'XTick';
        tickLabelStr = 'XTickLabel';
        tickLabelModeStr = 'XTickLabelMode';
    case 'y'
        tickStr = 'YTick';        
        tickLabelStr = 'YTickLabel';
        tickLabelModeStr = 'YTickLabelMode';
    otherwise
        error('axisType unrecognized!');
end

% Update ticks and tick labels
if strcmp(ax.(tickLabelModeStr), 'manual')
    tickLocsOrig = get(ax, tickStr);
    tickLabelsOrig = get(ax, tickLabelStr);
    yTickLabels = match_positions(tickLabelsOrig, tickLocsOrig, ...
                                    tickLocs);
    set(ax, tickStr, tickLocs);
    set(ax, tickLabelStr, yTickLabels);
else
    set(ax, tickStr, tickLocs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ax = sort_subplots (ax)
%% Sort subplots of a figure in the same order as the subplot function
% TODO: Pull out as its own function

if isempty(ax)
    return
end

% Extract outer positions
outerPositions = arrayfun(@(x) x.OuterPosition, ax, 'UniformOutput', false);

% Extract the left and bottom positions
[leftPositions, bottomPositions] = ...
    argfun(@(a) extract_elements(outerPositions, 'specific', 'Index', a), ...
            1, 2);

% Put together into a table
subplotTable = table(ax, outerPositions, leftPositions, bottomPositions);

% First sort the bottom positions by descending order, 
%   then sort the left positions by ascending order
subplotTable = sortrows(subplotTable, {'bottomPositions', 'leftPositions'}, ...
                        {'descend', 'ascend'});

% Extract the axes
ax = subplotTable.ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_color_to_rgb(object)
% TODO: Convert all grayscale to rgb?

if is_field(object, 'Color')
    % object.Color
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for iAx = 1:nAx
    ax(iAx).XAxis.LineWidth = 1;
    ax(iAx).YAxis.LineWidth = 1;
end

% Set other axes properties
if ~isempty(otherArguments)
    set(ax, otherArguments{:});    
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
