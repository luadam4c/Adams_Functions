function figHandle = update_figure_for_corel (figHandle, varargin)
%% Update figure to be journal-friendly (ready for CorelDraw)
% Usage: figHandle = update_figure_for_corel (figHandle, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       fig = update_figure_for_corel(fig);
%
% Outputs:
%       figHandle   - figure handle updated
%                   specified as a Figure object handle
%
% Arguments:
%       figHandle   - figure handle to update
%                   must be a Figure object handle
%       varargin    - 'RemoveTicks': whether to remove all ticks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RemoveLegends': whether to remove all legends
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for set_figure_properties()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_traces_spike2_mat.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-19 Created by Adam Lu
% 2019-09-20 Added 'RemoveTicks' as an optional argument
% 2019-10-04 Added 'RemoveLegends' as an optional argument
% 2019-10-05 Add textFontSize
% 2019-10-06 Changed default labelsFontSize from 8 to 7

%% Hard-coded parameters
% TODO: Make optional parameters
labelsFontSize = 7;
axisFontSize = 6; 
textFontSize = 6;
plotLineWidth = 0.025; %0.5;
annotationLineWidth = 1;
rulerLineWidth = 1;
units = 'inches';
tickLengthsUnits = [0.025, 0.025];

%% Default values for optional arguments
removeTicksDefault = false;  % set later
removeLegendsDefault = false;  % set later

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
addRequired(iP, 'figHandle');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RemoveTicks', removeTicksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RemoveLegends', removeLegendsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, figHandle, varargin{:});
removeTicks = iP.Results.RemoveTicks;
removeLegends = iP.Results.RemoveLegends;

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


%% Set figure properties
% Might change sizes
figHandle = set_figure_properties('FigHandle', figHandle, otherArguments);

% Update figure position units
unitsOrig = get(figHandle, 'Units');
if ~strcmp(unitsOrig, units)
    set(figHandle, 'Units', units);
end

%% Set axes properties
% Find all axes in the figure
ax = findall(figHandle, 'type', 'axes');

% Count the number of axes
nAx = numel(ax);

% Remove boxes
for iAx = 1:nAx
    box(ax(iAx), 'off');
end

% Make all ticks go outward
set(ax, 'TickDir', 'out', 'TickDirMode', 'manual');

% Remove x and y axis ticks
if removeTicks
    set(ax, 'XTick', []);
    set(ax, 'YTick', []);
end

% Remove legends if requested
if removeLegends
    lgds = findobj(gcf, 'Type', 'Legend');
    delete(lgds);
end

% Set font
set(ax, 'FontName', 'Arial');
set(ax, 'FontSize', axisFontSize);
set(ax, 'TitleFontSizeMultiplier', titleFontSizeMultiplier);
set(ax, 'TitleFontWeight', 'normal');
set(ax, 'LabelFontSizeMultiplier', labelFontSizeMultiplier);

% Change the fontsize of texts
texts = findobj(figHandle, 'Type', 'Text');
if ~isempty(texts)
    set(texts, 'Fontsize', textFontSize);
end

% Set ruler linewidths
set(ax, 'LineWidth', rulerLineWidth);

% Update plot line widths
lines = findobj(figHandle, 'Type', 'Line');
plots = lines(arrayfun(@(x) is_plot(x), lines));

set(plots, 'LineWidth', plotLineWidth);
% Update annotation line widths
% TODO: How to distinguish?

% Set tick lengths
figPosition = get(figHandle, 'Position');
for iAx = 1:nAx
    % Get the length of the longest axis
    axPosition = get(ax(iAx), 'Position');
    axisLengthUnits = max(axPosition(3:4) .* figPosition(3:4));

    % Convert to units relative to longest axis
    tickLengthsRel = tickLengthsUnits / axisLengthUnits;

    % Set new tick lengths
    set(ax(iAx), 'TickLength', tickLengthsRel);
end

%% Restore things
% Restore figure position units
unitsNow = get(figHandle, 'Units');
if ~strcmp(unitsNow, unitsOrig)
    set(figHandle, 'Units', unitsOrig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isPlot = is_plot(lineObject)

% Get x, y and z data
xData = lineObject.XData;
yData = lineObject.YData;
zData = lineObject.ZData;

% The line is a plot if there are more than two data points
isPlot = numel(xData) > 2 || numel(yData) > 2 || numel(zData) > 2;

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