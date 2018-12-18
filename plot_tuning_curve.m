function h = plot_tuning_curve (pValues, readout, varargin)
%% Plot a 1-dimensional tuning curve
% Usage: h = plot_tuning_curve (pValues, readout, varargin)
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:
%       pValues     - column vector of parameter values
%                   must be a numeric vector
%       readout     - a readout matrix where each column is a readout vector
%                   must be a numeric 2-D array
%       varargin    - 'ColsToPlot': columns of the readout matrix to plot
%                   must be a numeric vector
%                   default == 1:size(readout, 2);
%                   - 'LineSpec': line specification
%                   must be a character array
%                   default == '-'
%                   - 'PisLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == [false, false];
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'YLimits': limits of y axis
%                   must be a 2-element increasing numeric vector
%                   default == []
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
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
%                   - 'SingleColor': color when colsToPlot == 1
%                   must be a 3-element vector
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nTraces == 1 
%                               'northeast' if nTraces is 2~9
%                               'eastoutside' if nTraces is 10+
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
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isfigtype.m
%       cd/islegendlocation.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/plot_struct.m
%       /media/adamX/RTCl/tuning_curves.m
%
% 2017-04-17 Moved from tuning_curves.m
% 2017-04-17 Simplified code
% 2017-04-17 Set default arguments
% 2017-04-17 Color map is now based on number of columns to plot
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-09-25 Made almost all arguments parameter-value pairs
% 2018-12-15 Added 'LineSpec' as a parameter-value pair argument
% 2018-12-18 Now uses iP.KeepUnmatched
%

%% Hard-coded parameters
pTickAngle = 60;                % x tick angle in degrees
lineWidth = 2;
                    
%% Default values for optional arguments
colsToPlotDefault = [];         % set later
lineSpecDefault = '-';
pislogDefault = [false, false];
xlimitsDefault = [];
ylimitsDefault = [];
pTicksDefault = [];
pTickLabelsDefault = {};
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
columnLabelsDefault = '';       % set later
singleColorDefault = [0, 0, 1];
legendLocationDefault = 'auto'; % set later
figTitleDefault = '';           % set later
figNumberDefault = [];          % invisible figure by default
figNameDefault = '';
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
addRequired(iP, 'pValues', ...              % vector of parameter values
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColsToPlot', colsToPlotDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PisLog', pislogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XLimits', xlimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', ylimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) iscellstr(x) || isstring(x));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColumnLabels', columnLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'SingleColor', singleColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) isempty(x) || isnumeric(x) && isscalar(x) && x > 0);
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
colsToPlot = iP.Results.ColsToPlot;
lineSpec = iP.Results.LineSpec;
pIsLog = iP.Results.PisLog;
xlimits = iP.Results.XLimits;
ylimits = iP.Results.YLimits;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
columnLabels = iP.Results.ColumnLabels;
singlecolor = iP.Results.SingleColor;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
        numel(pTicks) ~= numel(pTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    h = [];
    return
end

%% Prepare for tuning curve
% Count number of entries
nEntries = length(pValues);

% Count number of columns
nCols = size(readout, 2);

% Set default columns to plot
if isempty(colsToPlot)
    colsToPlot = 1:size(readout, 2);
end

% Set column labels
if isempty(columnLabels)
    columnLabels = cell(1, nCols);
    for c = 1:nCols
        columnLabels{c} = ['Column #', num2str(c)];
    end
end

% Set legend location based on number of traces
if strcmpi(legendLocation, 'auto')
    if nCols > 1 && nCols < 10
        legendLocation = 'northeast';
    elseif nCols >= 10
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

% Set the default figure title
if isempty(figTitle)
    if ~strcmpi(readoutLabel, 'suppress') && ~strcmpi(pLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. ', pLabel], '_', '\_');
    elseif ~strcmpi(readoutLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. parameter'], '_', '\_');
    else
        figTitle = 'Readout vs. parameter';
    end
end

% Compute the number of parameter values that 
%   don't give infinite values
nNonInf = sum(~isinf(readout), 1);

% Count the number of columns to plot
nColsToPlot = length(colsToPlot);

% Define the color map to use
cm = colormap(jet(nColsToPlot));

if ~isempty(figName)
    % Create an invisible figure and clear it
    if ~isempty(figNumber)
        h = figure(figNumber);
        set(h, 'Visible', 'Off');
    else
        h = figure(floor(rand()*10^4)+1);
        set(h, 'Visible', 'Off');
    end
    clf(h);
else
    % Get the current figure
    h = gcf;
end

%% Plot tuning curve
% Hold on if more than one column
if nColsToPlot > 1
    hold on
end

% Plot readout values against parameter values
for c = 1:nColsToPlot
    % Get the column to plot
    col = colsToPlot(c);

    % Plot curve
    if pIsLog
        % Note: can't have hold on before semilogx
        p = semilogx(pValues, readout(:, col), lineSpec, ...
                        'LineWidth', lineWidth, otherArguments);
    else
        p = plot(pValues, readout(:, col), lineSpec, ...
                        'LineWidth', lineWidth, otherArguments);
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

% Hold off if more than one column
if nColsToPlot > 1
    hold off
end

% Generate a legend if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
    legend('location', legendLocation);
end

% Restrict x axis if xlimits provided; 
%   otherwise expand the x axis by a little bit
if ~isempty(xlimits)
    if ~strcmpi(xlimits, 'suppress')
        % Use x limits
        xlim(xlimits);
    end
else
    if nEntries > 1
        xlim([pValues(1) - (pValues(2) - pValues(1)), ...
            pValues(end) + (pValues(end) - pValues(end-1))]);
    end
end

% Restrict y axis if ylimits provided
if ~isempty(ylimits)
    ylim(ylimits);
end

% Set title and axes labels
if ~isempty(pTicks)
    set(gca, 'XTick', pTicks);
    % xticks(pTicks);
end
if ~isempty(pTickLabels)
    set(gca, 'XTickLabel', pTickLabels);
    % xticklabels(pTicks);
end
if ~strcmpi(pLabel, 'suppress')
    xlabel(pLabel);
end
if ~strcmpi(readoutLabel, 'suppress')
    ylabel(readoutLabel);
end
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

% Rotate p tick labels if too long
% TODO
xtickangle(pTickAngle);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
