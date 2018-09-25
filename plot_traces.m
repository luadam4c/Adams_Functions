function h = plot_traces(tVec, data, varargin)
%% Plots traces all in one place, overlapped or in parallel
% Usage: h = plot_traces(tVec, data, varargin)
% Outputs:
%       h           - figure handle for the created figure
%                   must be a figure handle
%
% Arguments:
%       tVec        - time vector for plotting
%                   must be a numeric vector
%       data        - data array (each column is a data vector)
%                   must be a numeric 2-D array
%       varargin    - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(tVec), max(tVec)]
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Time'
%                   - 'YLabel': label for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Data'
%                   - 'TraceLabels': labels for the traces, 
%                               suppress by setting value to 'suppress'
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Trace #1', 'Trace #2', ...}
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/plot_traces_abf.m

% File History:
% 2018-09-18 Moved from plot_traces_abf.m
% 2018-09-25 Implemented the input parser
% TODO: Combine with plot_signals.m by defining a parameter PlotMode

%% Default values for optional arguments
xlimitsDefault = [];
ylimitsDefault = [];
xLabelDefault = 'Time';
yLabelDefault = 'Data';
traceLabelsDefault = '';                % set later
figTitleDefault = '';                   % set later
figNumberDefault = [];
figNameDefault = '';
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'tVec', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'XLimits', xlimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'YLimits', ylimitsDefault, ...
    @(x) ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TraceLabels', traceLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, tVec, data, varargin{:});
xlimits = iP.Results.XLimits;
ylimits = iP.Results.YLimits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
traceLabels = iP.Results.TraceLabels;
figTitle = iP.Results.FigTitle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Preparation
% Extract number of traces
nTraces = size(data, 2);

% Set trace labels
if isempty(traceLabels)
    traceLabels = cell(1, nTraces);
    for iTrace = 1:nTraces
        traceLabels{c} = ['Trace #', num2str(iTrace)];
    end
end

% Check if traceLabels has the correct length
if iscell(traceLabels) && numel(traceLabels) ~= nTraces
    error('traceLabels has %d elements instead of %d!!', ...
            numel(traceLabels), nTraces);
end

% Set the default figure title
if isempty(figTitle)
    if ~isempty(figName)
        figTitle = ['Traces for ', figName];
    else
        figTitle = [yLabel, ' over ', xLabel];
    end
end

% Compute minimum and maximum Y values
minY = min(min(data));
maxY = max(max(data));
rangeY = maxY - minY;

% Decide on the figure to plot on
if ~isempty(figName)
    % Create an invisible figure and clear it
    if ~isempty(figNumber)
        h = figure(figNumber);
        set(h, 'Visible', 'Off');
    else
        h = figure('Visible', 'Off');
    end
    clf(h);
else
    % Get the current figure
    h = gcf;
end

%% Plot
% Plot all traces
hold on;
for iTrace = 1:nTraces
    plot(tVec, data(:, iTrace), ...
        'DisplayName', traceLabels{iTrace});
end

% Determine the appropriate time axis limits
if ~isempty(xlimits)
    if ~strcmpi(xlimits, 'suppress')
        xlim(xlimits);
    end
else
    xlim([min(tVec), max(tVec)]);
end


% Determine the appropriate y axis limits
if ~isempty(ylimits)
    if ~strcmpi(ylimits, 'suppress')
        % Use x limits
        ylim(ylimits);
    end
else
    if rangeY ~= 0
        ylim([minY - 0.2 * rangeY, maxY + 0.2 * rangeY]);
    end
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
    title(figTitle, 'Interpreter', 'none');
end

% Generate a legend if there is more than one trace
if nTraces > 1 && nTraces < 10
    legend('location', 'northeast');
elseif nTraces >= 10
    legend('location', 'eastoutside');
end

%% Save
% Save figure
if ~isempty(figName)
    save_all_figtypes(h, figName, figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

function h = plot_traces(tVec, data, xlimits, xLabel, yLabel, ...
                            traceLabels, figTitle, figName, figNum)
%       xlimits     - x-axis limits
%       xLabel      - x-axis label
%       yLabel      - y-axis label
%       traceLabels - legend labels for each trace
%       figTitle    - figure title
%       figName     - figure name
%       figNum      - figure number

% Hold off and close figure
hold off;
close(h);

saveas(h, figName, 'png');

h = figure(figNum);
set(h, 'Visible', 'Off');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
