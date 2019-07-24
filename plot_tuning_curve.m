function [fig, lines, boundaries] = plot_tuning_curve (pValues, readout, varargin)
%% Plot a 1-dimensional tuning curve
% Usage: [fig, lines, boundaries] = plot_tuning_curve (pValues, readout, varargin)
% Outputs:
%       fig         - figure handle for the created figure
%                   specified as a figure object handle
%       lines       - lines for the tuning curves
%                   specified as a line object handle array
%       boundaries  - boundary lines
%                   specified as a line object handle array
% Arguments:
%       pValues     - column vector of parameter values
%                   must be a numeric vector
%       readout     - a readout matrix where each column is a readout vector
%                   must be a numeric 2-D array
%       varargin    - 'PhaseVectors': phase information for each readout vector
%                   must be a numeric matrix or a cell array of numeric vectors
%                   default == {}
%                   - 'RemoveOutliers': whether to remove outliers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ColsToPlot': columns of the readout matrix to plot
%                   must be a numeric vector
%                   default == 1:size(readout, 2);
%                   - 'LineSpec': line specification
%                   must be a character array
%                   default == '-'
%                   - 'PisLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PLimits': limits of parameter axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'ReadoutLimits': limits of readout axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == []
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'PTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
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
%                   - 'PBoundaries': parameter boundary values
%                   must be a numeric vector
%                   default == []
%                   - 'RBoundaries': readout boundary values
%                   must be a numeric vector
%                   default == []
%                   - 'IndSelected': selected indices to mark differently
%                   must be a numeric vector
%                   default == []
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
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
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       ~/Downloaded_Functions/rgb.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_fighandle.m
%       cd/force_column_vector.m
%       cd/isfigtype.m
%       cd/islegendlocation.m
%       cd/plot_horizontal_line.m
%       cd/plot_vertical_line.m
%       cd/remove_outliers.m
%       cd/save_all_figtypes.m
%       cd/unique_custom.m
%
% Used by:
%       cd/plot_struct.m
%       /media/adamX/RTCl/tuning_curves.m

% 2017-04-17 Moved from tuning_curves.m
% 2017-04-17 Simplified code
% 2017-04-17 Set default arguments
% 2017-04-17 Color map is now based on number of columns to plot
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-09-25 Made almost all arguments parameter-value pairs
% 2018-12-15 Added 'LineSpec' as a parameter-value pair argument
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Changed lineSpec default to o and singleColorDefault to SkyBlue
% 2019-03-14 Added 'RemoveOutliers' as an optional argument
% 2019-03-25 Added 'PhaseVectors' as an optional argument
% 2019-03-25 Now expands the y limits by a little by default
% 2019-05-10 Now uses decide_on_fighandle.m
% 2019-06-10 Added 'PBoundaries' and 'RBoundaries' as optional arguments
% TODO: Make 'ColorMap' and optional argument
% TODO: Return handles to plots
%

%% Hard-coded parameters

%% Default values for optional arguments
phaseVectorsDefault = {};           % no phase vectors by default
removeOutliersDefault = false;      % don't remove outliers by default
colsToPlotDefault = [];             % set later
lineSpecDefault = '-';
lineWidthDefault = 2;
pislogDefault = false;
pLimitsDefault = [];
readoutLimitsDefault = [];
pTicksDefault = [];
pTickLabelsDefault = {};
pTickAngleDefault = [];             % set later
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
columnLabelsDefault = '';           % set later
singleColorDefault = rgb('SkyBlue');
legendLocationDefault = 'auto';     % set later
pBoundariesDefault = [];
rBoundariesDefault = [];
indSelectedDefault = [];
figTitleDefault = '';               % set later
figHandleDefault = [];              % no existing figure by default
figNumberDefault = [];              % no figure number by default
figNameDefault = '';                % don't save figure by default
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
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, {'vector'}));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PhaseVectors', phaseVectorsDefault, ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['PhaseVectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'RemoveOutliers', removeOutliersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ColsToPlot', colsToPlotDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'PisLog', pislogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PLimits', pLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ReadoutLimits', readoutLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PTickAngle', pTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
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
addParameter(iP, 'PBoundaries', pBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RBoundaries', rBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'IndSelected', indSelectedDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
phaseVectors = iP.Results.PhaseVectors;
removeOutliers = iP.Results.RemoveOutliers;
colsToPlot = iP.Results.ColsToPlot;
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
pIsLog = iP.Results.PisLog;
pLimits = iP.Results.PLimits;
readoutLimits = iP.Results.ReadoutLimits;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pTickAngle = iP.Results.PTickAngle;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
columnLabels = iP.Results.ColumnLabels;
singlecolor = iP.Results.SingleColor;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
pBoundaries = iP.Results.PBoundaries;
rBoundaries = iP.Results.RBoundaries;
indSelected = iP.Results.IndSelected;
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
        numel(pTicks) ~= numel(pTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    fig = [];
    return
end

%% Prepare for tuning curve
% Count number of entries
nEntries = length(pValues);

% Count number of columns
nCols = size(readout, 2);

% Deal with phase vectors
if ~isempty(phaseVectors)
    % Force as a cell array of column vectors
    phaseVectors = force_column_vector(phaseVectors, 'ForceCellOutput', true);

    % Make sure the number of phase vectors matchs the number of readout columns
    phaseVectors = match_row_count(phaseVectors, nCols);

    % Get the unique phases for each readout column
    uniquePhases = cellfun(@(x) unique_custom(x, 'stable', ...
                                                'IgnoreNaN', true), ...
                            phaseVectors, 'UniformOutput', false);

    % Count the number of phases for each readout column
    nPhases = count_samples(uniquePhases);

    % Compute the maximum number of phases
    maxNPhases = max(nPhases);

    % Create phase labels
    phaseLabels = create_labels_from_numbers(1:maxNPhases, 'Prefix', 'Phase #');
else
    uniquePhases = {};
    nPhases = [];
    maxNPhases = 1;   
    phaseLabels = {};
end

% Remove outliers if requested
if removeOutliers
    readout = remove_outliers(readout, 'OutlierMethod', 'isoutlier', ...
                                'ReplaceWithNans', true);
end

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

% Count the number of boundaries
nPBoundaries = size(pBoundaries, 2);
nRBoundaries = size(rBoundaries, 2);
nBoundaries = nPBoundaries + nRBoundaries;

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
if isempty(phaseVectors)
    cm = colormap(jet(nColsToPlot));
else
    % Create a color map
    cm = colormap(hsv(maxNPhases));
end

% Decide on the figure to plot on
fig = decide_on_fighandle('FigHandle', figHandle, 'FigNumber', figNumber);

% Set the default parameter tick angle
if isempty(pTickAngle)
    if ~isempty(pTickLabels)
        % TODO: Rotate p tick labels only if too long
        pTickAngle = 60;
    else
        pTickAngle = 0;
    end
end

%% Plot tuning curve
% Hold on if more than one column
if nColsToPlot > 1
    hold on
end

% Plot readout values against parameter values
lines = gobjects(nColsToPlot, maxNPhases);
for c = 1:nColsToPlot
    % Get the column to plot
    col = colsToPlot(c);

    % Plot the tuning curve for this column
    if isempty(phaseVectors)
        lines(c, 1) = plot_one_line(pIsLog, pValues, readout(:, col), lineSpec, ...
                            lineWidth, otherArguments);
    else
        hold on;
        
        % Get the current phase vector
        phaseVectorThis = phaseVectors{c};
        uniquePhasesThis = uniquePhases{c};
        nPhasesThis = nPhases(c);

        % Get the distinct phase indices for this readout vector
        phaseIndices = arrayfun(@(x) find(phaseVectorThis == x), ...
                                uniquePhasesThis, 'UniformOutput', false);

        % Get the last index for this readout vector
        lastIndex = numel(phaseVectorThis);

        % Add the next index
        phaseIndices = cellfun(@(x) add_next_index(x, lastIndex), ...
                                phaseIndices, 'UniformOutput', false);

        % Plot readout vector for all phases
        lines(c, 1:nPhasesThis) = ...
            cellfun(@(x) plot_one_line(pIsLog, pValues(x), readout(x, col), ...
                        lineSpec, lineWidth, otherArguments), phaseIndices);
    end

    % Set color
    if isempty(phaseVectors)
        % Set color by columns
        if nColsToPlot > 1
            set(lines(c, 1), 'Color', cm(c, :))
        elseif nColsToPlot == 1
            set(lines(c, 1), 'Color', singlecolor);
        end
    else
        % Set color by phase
        for iPhase = 1:nPhasesThis
            set(lines(c, iPhase), 'Color', cm(iPhase, :));
        end
    end

    % Set display name
    if isempty(phaseVectors)
        if ~strcmpi(columnLabels, 'suppress')
            set(lines(c, 1), 'DisplayName', ...
                replace(columnLabels{col}, '_', '\_'));
        end
    else
        for iPhase = 1:nPhasesThis
            set(lines(c, iPhase), 'DisplayName', ...
                replace(phaseLabels{iPhase}, '_', '\_'));
        end
    end

    % If there is only one value for this column, mark with a circle
    % TODO: for ~isempty(phaseVectors)?
    if nNonInf(col) == 1
        set(lines(c, 1), 'Marker', 'o');
    end
end

% Hold off if more than one column
if nColsToPlot > 1
    hold off
end

% Generate a legend if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
    if isempty(phaseVectors)
        legend('location', legendLocation);
    else
        legend(lines(1, :), 'location', legendLocation);
    end
end

% Restrict x axis if pLimits provided; 
%   otherwise expand the x axis by a little bit
if ~isempty(pLimits)
    if ~strcmpi(pLimits, 'suppress')
        % Use x limits
        xlim(pLimits);
    end
else
    if nEntries > 1 && nEntries < 4
        xlim(compute_axis_limits(pValues, 'x', 'Coverage', 90));
    elseif nEntries >= 4
        xlim([pValues(1) - (pValues(2) - pValues(1)), ...
            pValues(end) + (pValues(end) - pValues(end-1))]);
    end
end

% Restrict y axis if readoutLimits provided
%   otherwise expand the y axis by a little bit
if ~isempty(readoutLimits)
    ylim(readoutLimits);
else
    ylim(compute_axis_limits(get(gca, 'YLim'), 'y', 'Coverage', 80));
end

% Set title and axes labels
if ~isempty(pTicks)
    xticks(pTicks);
end
if ~isempty(pTickLabels)
    xticklabels(pTicks);
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

% Set the angle for parameter ticks
xtickangle(pTickAngle);

% Plot parameter boundaries
if nPBoundaries > 0
    hold on
    pLines = plot_vertical_line(pBoundaries, 'LineWidth', 0.5, ...
                                'LineStyle', '--', 'Color', 'g');
else
    pBoundaries = gobjects;
end

% Plot readout boundaries
if nRBoundaries > 0
    hold on
    rLines = plot_horizontal_line(rBoundaries, 'LineWidth', 0.5, ...
                                'LineStyle', '--', 'Color', 'r');
else
    rBoundaries = gobjects;
end

% Plot selected values if any
if ~isempty(indSelected)
    % TODO
end

%% Post-plotting
% Return boundaries
boundaries = transpose(vertcat(pBoundaries, rBoundaries));

% Save figure if figName provided
if ~isempty(figName)
    save_all_figtypes(fig, figName, figtypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = plot_one_line(pIsLog, pValues, readout, lineSpec, ...
                            lineWidth, otherArguments)

if pIsLog
    % Note: can't have hold on before semilogx
    p = semilogx(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
else
    p = plot(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = add_next_index(indices, lastIndex)
%% Add the next index if not the last
% TODO: What if indices not contiguous?

if indices(end) < lastIndex
    indices = [indices; indices(end) + 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Usage: plot_tuning_curve(pValues, readout, colsToPlot, pIsLog, pLabel, ...
            readoutLabel, columnLabels, pLimits, readoutLimits, figName, varargin)

if ~isequal(columnLabels, {'suppress'})

if isequal(pLimits, -1)

if ~isequal(pLabel, 'suppress')
if ~isequal(readoutLabel, 'suppress')
if ~isequal(pLabel, 'suppress') && ~isequal(readoutLabel, 'suppress')

singleColorDefault = [0, 0, 1];
lineSpecDefault = '-';

set(fig, 'Visible', 'Off');
fig = figure(floor(rand()*10^4)+1);

if pIsLog
    % Note: can't have hold on before semilogx
    p = semilogx(pValues, readout(:, col), lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
else
    p = plot(pValues, readout(:, col), lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
end

if ~isempty(figName)
    % Create a figure
    if ~isempty(figNumber)
        % Create an invisible figure
        fig = figure(figNumber);
        set(fig, 'Visible', 'Off');
    else
        % Create a new figure
        fig = figure;
    end

    % Clear the figure
    clf(fig);
else
    % Get the current figure
    fig = gcf;
end

set(gca, 'XTick', pTicks);
set(gca, 'XTickLabel', pTickLabels);
pTickAngle = 60;                % x tick angle in degrees

% TODO FOR UNDERGRAD: Create unique_custom.m with the option 'IgnoreNan'
phaseVectorsNoNaN = cellfun(@(x) x(~isnan(x)), phaseVectors, ...
                            'UniformOutput', false);
uniquePhases = cellfun(@(x) unique(x, 'stable'), phaseVectorsNoNaN, ...
                        'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
