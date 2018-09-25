function h = plot_tuning_curve(pValues, readout, varargin)
%% Plot a 1-dimensional tuning curve
% Usage: h = plot_tuning_curve(pValues, readout, varargin)
% Outputs:
%       h           - figure handle for the created figure
%                   must be a figure handle
% Arguments:
%       pValues     - column vector of parameter values
%                   must be a numeric vector
%       readout     - a readout matrix where each column is a readout vector
%                   must be a numeric 2-D array
%       varargin    - 'ColsToPlot': columns of the readout matrix to plot
%                   must be a numeric vector
%                   default == 1:size(readout, 2);
%                   - 'PisLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == [false, false];
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'ColumnLabels': labels for the readout columns, 
%                               suppress by setting value to {'suppress'}
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Column #1', 'Column #2', ...}
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'YLimits': limits of y axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'SingleColor': color when colsToPlot == 1
%                   must be a 3-element vector
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
%       /media/adamX/RTCl/tuning_curves.m
%
% 2017-04-17 Moved from tuning_curves.m
% 2017-04-17 Simplified code
% 2017-04-17 Set default arguments
% 2017-04-17 Color map is now based on number of columns to plot
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-09-25 Made almost all arguments parameter-value pairs
%

%% Default values for optional arguments
pislogDefault = [false, false];
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
columnLabelsDefault = '';               % set later
xlimitsDefault = [];
ylimitsDefault = [];
singleColorDefault = [0, 0, 1];
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
addRequired(iP, 'pValues', ...              % vector of parameter values
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColsToPlot', colsToPlotDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'PisLog', pislogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColumnLabels', columnLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'XLimits', xlimitsDefault, ...
    @(x) ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', ylimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, ...
                            {'increasing', 'vector', 'numel', 2}));
addParameter(iP, 'SingleColor', singleColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
colsToPlot = iP.Results.ColsToPlot;
pIsLog = iP.Results.PisLog;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
columnLabels = iP.Results.ColumnLabels;
xlimits = iP.Results.XLimits;
ylimits = iP.Results.YLimits;
singlecolor = iP.Results.SingleColor;
figName = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

%% Prepare for tuning curve
% Extract number of columns
nCols = size(readout, 2);

% Set column labels
if isempty(columnLabels)
    columnLabels = cell(1, nCols);
    for c = 1:nCols
        columnLabels{c} = ['Column #', num2str(c)];
    end
end

% Compute the number of parameter values that 
%   don't give infinite values
nNonInf = sum(~isinf(readout), 1);

% Count the number of columns to plot
nColsToPlot = length(colsToPlot);

% Define the color map to use
cm = colormap(jet(nColsToPlot));

% Decide on the figure to plot on
if ~isempty(figName)
    % Create and clear figure if figName provided
    h = figure(floor(rand()*10^4)+1);
    clf(h);
else
    % Get the current figure
    h = gcf;
end

%% Plot tuning curve
% Plot readout values against parameter values
for c = 1:nColsToPlot
    % Get the column to plot
    col = colsToPlot(c);

    % Plot curve
    if pIsLog
        % Note: can't have hold on before semilogx
        p = semilogx(pValues, readout(:, col)); hold on;
    else
        p = plot(pValues, readout(:, col)); hold on;
    end
    
    % Set color
    if nColsToPlot > 1
        set(p, 'Color', cm(c, :))
    elseif nColsToPlot == 1
        set(p, 'Color', singlecolor);
    end

    % Set display name
    if ~strcmpi(columnLabels, 'suppress')
        set(p, 'DisplayName', strrep(columnLabels{col}, '_', '\_'));
    end

    % If there is only one value for this column, mark with a circle
    if nNonInf(col) == 1
        set(p, 'Marker', 'o');
    end
end

% Show legend only if readout has more than one columns
if nCols > 1
    legend('Location', 'eastoutside');
end

% Restrict x axis if xlimits provided; 
%   otherwise expand the x axis by a little bit
if ~isempty(xlimits)
    if ~strcmpi(xlimits, 'suppress')
        % Use x limits
        xlim(xlimits);
    end
else
    xlim([pValues(1) - (pValues(2) - pValues(1)), ...
        pValues(end) + (pValues(end) - pValues(end-1))]);
end

% Restrict y axis if ylimits provided
if ~isempty(ylimits)
    ylim(ylimits);
end

% Set title and axes labels
if ~strcmpi(pLabel, 'suppress')
    xlabel(pLabel);
end
if ~strcmpi(readoutLabel, 'suppress')
    ylabel(readoutLabel);
end
if ~strcmpi(pLabel, 'suppress') && ~strcmpi(readoutLabel, 'suppress')
    title(strrep([readoutLabel, ' vs. ', pLabel], '_', '\_'));
end

%% Post-plotting
% Save figure if figName provided
if ~isempty(figName)
    save_all_figtypes(h, figName, figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Usage: plot_tuning_curve(pValues, readout, colsToPlot, pIsLog, pLabel, readoutLabel, columnLabels, xlimits, ylimits, figName, varargin)

if ~isequal(columnLabels, {'suppress'})

if isequal(xlimits, -1)

if ~isequal(pLabel, 'suppress')
if ~isequal(readoutLabel, 'suppress')
if ~isequal(pLabel, 'suppress') && ~isequal(readoutLabel, 'suppress')

%}
