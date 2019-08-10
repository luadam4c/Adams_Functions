function [fig, lines, boundaries] = plot_tuning_curve (pValues, readout, varargin)
%% Plot 1-dimensional tuning curve(s), can include confidence intervals or test p values
% Usage: [fig, lines, boundaries] = plot_tuning_curve (pValues, readout, varargin)
% Explanation:
%       TODO
% Examples:
%       pValue = transpose(1:10);
%       readout1 = randi(numel(pValue), 10, 1);
%       upperCI1 = readout1 + randi(numel(pValue), 10, 1) / 10;
%       lowerCI1 = readout1 - randi(numel(pValue), 10, 1) / 10;
%       readout2 = randi(numel(pValue), 10, 1);
%       upperCI2 = readout2 + randi(numel(pValue), 10, 1) / 10;
%       lowerCI2 = readout2 - randi(numel(pValue), 10, 1) / 10;
%       readoutAll = [readout1, readout2];
%       upperCIAll = [upperCI1, upperCI2];
%       lowerCIAll = [lowerCI1, lowerCI2];
%       
%       plot_tuning_curve(pValue, readout1, 'UpperCI', upperCI1, 'LowerCI', lowerCI1);
%       plot_tuning_curve(pValue, readoutAll, 'UpperCI', upperCIAll, 'LowerCI', lowerCIAll, 'ColorMap', hsv(2));
%
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
%       varargin    - 'LowerCI': lower bounds of confidence intervals
%                   must be a numeric 2-D array
%                   default == []
%                   - 'UpperCI': upper bounds of confidence intervals
%                   must be a numeric 2-D array
%                   default == []
%                   - 'PhaseVectors': phase information for each readout vector
%                   must be a numeric matrix or a cell array of numeric vectors
%                   default == {}
%                   - 'RemoveOutliers': whether to remove outliers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunTTest': whether to run paired t-test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunRankTest': whether to run paired 
%                                       Wilcoxon signed-rank test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ColumnsToPlot': columns of the readout matrix to plot
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
%                   - 'PhaseLabels': phase labels if phase vectors are provided
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Phase #1', 'Phase #2', ...}
%                   - 'ColorMap' - color map used when nColumnsToPlot > 1
%                   must be a 2-D numeric array with 3 columns
%                   default == jet(nColumnsToPlot) or 
%                               rgb('SkyBlue') == [0.5273, 0.8047, 0.9180]
%                                   if nColumnsToPlot == 1 or
%                               hsv(maxNPhases) if phaseVectors is provided
%                   - 'ConfIntColor': color for confidence intervals
%                   must be a 3-element vector
%                   default == rgb('LightGray') == [0.8242, 0.8242, 0.8242]
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
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
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
%       cd/count_samples.m
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
%       ~/Marks_Functions/Adam/CLC2/markCLC2figures.m
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
% 2019-08-07 Now changes the pTickAngle only if the labels are too long
% 2019-08-07 Added 'PhaseLabels' as an optional argument
% 2019-08-07 Added 'ColorMap' as an optional argument
% 2019-08-07 Added 'ClearFig' as an optional argument
% 2019-08-07 Now accepts infinite values for readout limits
% 2019-08-07 Added 'RunTTest' as an optional argument
% 2019-08-07 Added 'RunRankTest' as an optional argument
% 2019-08-09 Updated confidence interval plots
% 2019-08-09 Combined SingleColor with ColorMap
% 2019-08-09 Fixed confidence interval plots for matrices
% 2019-08-09 Added 'IndSelected' as an optional argument
% TODO: Return handles to plots
%

%% Hard-coded constants
WHITE = [1, 1, 1];

%% Hard-coded parameters
sigLevel = 0.05;                    % significance level for tests
fadePercentage = 0.25;              % fade percentage for confidence interval colors
markerLineWidth = 3;

%% Default values for optional arguments
lowerCIDefault = [];
upperCIDefault = [];
phaseVectorsDefault = {};           % no phase vectors by default
removeOutliersDefault = false;      % don't remove outliers by default
runTTestDefault = false;            % don't run paired t-test by default
runRankTestDefault = false;         % don't run paired signed-rank test by default
columnsToPlotDefault = [];          % set later
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
phaseLabelsDefault = '';            % set later
colorMapDefault = [];               % set later
confIntColorDefault = [];           % light gray confidence intervals by default
legendLocationDefault = 'auto';     % set later
pBoundariesDefault = [];
rBoundariesDefault = [];
indSelectedDefault = [];
figTitleDefault = '';               % set later
figHandleDefault = [];              % no existing figure by default
figNumberDefault = [];              % no figure number by default
clearFigureDefault = false;         % don't clear figure by default
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
addParameter(iP, 'LowerCI', lowerCIDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'UpperCI', upperCIDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'PhaseVectors', phaseVectorsDefault, ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['PhaseVectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'RemoveOutliers', removeOutliersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunTTest', runTTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunRankTest', runRankTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ColumnsToPlot', columnsToPlotDefault, ...
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
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
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
addParameter(iP, 'PhaseLabels', phaseLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 3}));
addParameter(iP, 'ConfIntColor', confIntColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d', 'numel', 3}));
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
addParameter(iP, 'ClearFigure', clearFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
lowerCI = iP.Results.LowerCI;
upperCI = iP.Results.UpperCI;
phaseVectors = iP.Results.PhaseVectors;
removeOutliers = iP.Results.RemoveOutliers;
runTTest = iP.Results.RunTTest;
runRankTest = iP.Results.RunRankTest;
columnsToPlot = iP.Results.ColumnsToPlot;
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
phaseLabels = iP.Results.PhaseLabels;
colorMap = iP.Results.ColorMap;
confIntColor = iP.Results.ConfIntColor;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
pBoundaries = iP.Results.PBoundaries;
rBoundaries = iP.Results.RBoundaries;
indSelected = iP.Results.IndSelected;
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
clearFigure = iP.Results.ClearFigure;
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
    if isempty(phaseLabels)
        phaseLabels = create_labels_from_numbers(1:maxNPhases, ...
                                                'Prefix', 'Phase #');
    end
else
    uniquePhases = {};
    nPhases = [];
    maxNPhases = 1;
end

% Remove outliers if requested
if removeOutliers
    readout = remove_outliers(readout, 'OutlierMethod', 'isoutlier', ...
                                'ReplaceWithNans', true);
end

% Run paired t-tests if requested
if (runTTest || runRankTest) && size(readout, 1) > 1
    % Transpose the readout matrix
    readoutTransposed = transpose(readout);

    % Extract the before columns
    befores = readoutTransposed(:, 1:end-1);

    % Extract the after columns
    afters = readoutTransposed(:, 2:end);
end

% Run paired t-tests if requested
if runTTest
    % Run t-tests on each pair of columns
    [tTestResults, tTestPValues] = ttest(afters, befores);
else
    tTestResults = [];
    tTestPValues = [];
end

% Run paired Wilcoxon signed rank tests if requested
if runRankTest
    % Run t-tests on each pair of columns
    % TODO: Make a function vecfun.m
    [rankTestPValues, rankTestResults] = ...
        arrayfun(@(x) signrank(afters(:, x), befores(:, x)), 1:size(afters, 2));
else
    rankTestResults = [];
    rankTestPValues = [];
end

% Set default columns to plot
if isempty(columnsToPlot)
    columnsToPlot = 1:size(readout, 2);
end

% Set column labels
if isempty(columnLabels)
    columnLabels = cell(1, nCols);
    for iPlot = 1:nCols
        columnLabels{iPlot} = ['Column #', num2str(iPlot)];
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
nColumnsToPlot = length(columnsToPlot);

% Decide on the color map to use
if isempty(colorMap)
    if isempty(phaseVectors)
        if nColumnsToPlot == 1
            colorMap = rgb('SkyBlue');
        else
            colorMap = colormap(jet(nColumnsToPlot));
        end
    else
        colorMap = colormap(hsv(maxNPhases));
    end
end

% Decide on the confidence interval color map to use
if isempty(confIntColor)
    % Color of the confidence interval
    confIntColor = WHITE - (WHITE - colorMap) * fadePercentage;
end

% Decide on the figure to plot on
fig = decide_on_fighandle('FigHandle', figHandle, 'FigNumber', figNumber);

% Clear figure if requested
if clearFigure
    clf;
end

% Set the default parameter tick angle
if isempty(pTickAngle)
    if ~isempty(pTickLabels)
        maxCharPTickLabels = max(count_samples(pTickLabels, ...
                                    'TreatCellStrAsArray', false));
        if maxCharPTickLabels > 10
            pTickAngle = 60;
        elseif maxCharPTickLabels > 3
            pTickAngle = 45;
        else
            pTickAngle = 0;
        end
    else
        pTickAngle = 0;
    end
end

%% Plot tuning curve
% Hold on if more than one column
if nColumnsToPlot > 1
    hold on
end

% Plot readout values against parameter values
lines = gobjects(nColumnsToPlot, maxNPhases);
confInts = gobjects(nColumnsToPlot, maxNPhases);
for iPlot = 1:nColumnsToPlot
    % Get the column to plot
    col = columnsToPlot(iPlot);
    readoutThis = readout(:, col);

    % Plot the tuning curve for this column
    if isempty(phaseVectors)
        lines(iPlot, 1) = plot_one_line(pIsLog, pValues, readoutThis, ...
                                        lineSpec, lineWidth, otherArguments);
    else
        hold on;
        
        % Get the current phase vector
        phaseVectorThis = phaseVectors{iPlot};
        uniquePhasesThis = uniquePhases{iPlot};
        nPhasesThis = nPhases(iPlot);

        % Get the distinct phase indices for this readout vector
        phaseIndices = arrayfun(@(x) find(phaseVectorThis == x), ...
                                uniquePhasesThis, 'UniformOutput', false);

        % Get the last index for this readout vector
        lastIndex = numel(phaseVectorThis);

        % Add the next index
        phaseIndices = cellfun(@(x) add_next_index(x, lastIndex), ...
                                phaseIndices, 'UniformOutput', false);

        % Plot readout vector for all phases
        lines(iPlot, 1:nPhasesThis) = ...
            cellfun(@(x) plot_one_line(pIsLog, pValues(x), readout(x, col), ...
                        lineSpec, lineWidth, otherArguments), phaseIndices);
    end

    % If provided, plot a confidence interval for this column
    %   as a light-gray-shaded area
    if ~isempty(lowerCI) || ~isempty(upperCI)
        if ~isempty(lowerCI)
            lowerCIThis = lowerCI(:, col);
        else
            lowerCIThis = readoutThis;
        end
        if ~isempty(upperCI)
            upperCIThis = upperCI(:, col);
        else
            upperCIThis = readoutThis;
        end

        % Plot the confidence interval
        if ~isempty(phaseVectors) || pIsLog
            fprintf('Not Supported Yet!\n');
        else
            hold on;

            % Get the current Y limits
            yLimits = get(gca, 'YLim');

            % Compute the minimum y limits
            %TODO
            minY = apply_iteratively(@min, {yLimits, readoutLimits});

            % Remove tuning curve
            delete(lines(iPlot, 1));

            % The x and y values for the confidence intervals
            confIntXValues = [pValues; flipud(pValues)];
            confIntYValues = [upperCIThis; flipud(lowerCIThis)];

            % Fill the area between lowerCIThis and upperCIThis 
            %   with confIntColor
            confInts(iPlot, 1) = fill(confIntXValues, confIntYValues, ...
                                confIntColor(iPlot, :), 'LineStyle', 'none');

            % Plot tuning curve again
            lines(iPlot, 1) = plot_one_line(pIsLog, pValues, readoutThis, ...
                                        lineSpec, lineWidth, otherArguments);

            % Display tick marks and grid lines over graphics objects.
            set(gca, 'Layer', 'top');
        end
    end

    % Set color
    if isempty(phaseVectors)
        % Set color by columns
        set(lines(iPlot, 1), 'Color', colorMap(iPlot, :));
    else
        % Set color by phase
        for iPhase = 1:nPhasesThis
            set(lines(iPlot, iPhase), 'Color', colorMap(iPhase, :));
        end
    end

    % Set display name
    if isempty(phaseVectors) || isempty(phaseLabels)
        if ~strcmpi(columnLabels, 'suppress')
            set(lines(iPlot, 1), 'DisplayName', ...
                replace(columnLabels{col}, '_', '\_'));
        end
    else
        for iPhase = 1:nPhasesThis
            set(lines(iPlot, iPhase), 'DisplayName', ...
                replace(phaseLabels{iPhase}, '_', '\_'));
        end
    end

    % If there is only one value for this column, mark with a circle
    % TODO: for ~isempty(phaseVectors)?
    if nNonInf(col) == 1
        set(lines(iPlot, 1), 'Marker', 'o');
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
    yLimitsOrig = get(gca, 'YLim');
    if isinf(readoutLimits(1))
        readoutLimits(1) = yLimitsOrig(1);
    end
    if isinf(readoutLimits(2))
        readoutLimits(2) = yLimitsOrig(2);
    end
    ylim(readoutLimits);
else
    ylim(compute_axis_limits(get(gca, 'YLim'), 'y', 'Coverage', 80));
end

% Set title and axes labels
if ~isempty(pTicks)
    xticks(pTicks);
end
if ~isempty(pTickLabels)
    xticklabels(pTickLabels);
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
    % Selected x locations
    xLocsSelected = pValues(indSelected);

    % Selected y locations
    yLocsSelected = readout(indSelected, :);

    % Plot values
    hold on
    markers = plot(xLocsSelected, yLocsSelected, 'ro', ...
                    'LineWidth', markerLineWidth);
end

% TODO: Make function plot_text.m
% Plot t-test p values if any
if ~isempty(tTestPValues)
    hold on

    % Get the x locations
    xLocs = pValues(1:end-1) + (pValues(2) - pValues(1)) * 0.25;
    
    % Get current y axis limits
    yLimitsNow = get(gca, 'YLim');

    % Get y location
    yLoc = yLimitsNow(1) + (yLimitsNow(2) - yLimitsNow(1)) * 0.1;

    % Plot texts
    for iValue =  1:numel(tTestPValues)
        % Get the current values
        tTestPValueThis = tTestPValues(iValue);
        xLocThis = xLocs(iValue);

        % Create a p value string to 2 significant digits
        pString = ['p_t = ', num2str(tTestPValueThis, 2)];

        % Plot red if significant
        if tTestPValueThis < sigLevel
            text(xLocThis, yLoc, pString, 'Color', 'r');
        else
            text(xLocThis, yLoc, pString, 'Color', 'k');
        end
    end
end

% Plot rank test p values if any
if ~isempty(rankTestPValues)
    hold on

    % Get the x locations
    xLocs = pValues(1:end-1) + (pValues(2) - pValues(1)) * 0.25;
    
    % Get current y axis limits
    yLimitsNow = get(gca, 'YLim');

    % Get y location
    yLoc = yLimitsNow(1) + (yLimitsNow(2) - yLimitsNow(1)) * 0.05;

    % Plot texts
    for iValue =  1:numel(rankTestPValues)
        % Get the current values
        rankTestPValueThis = rankTestPValues(iValue);
        xLocThis = xLocs(iValue);

        % Create a p value string to 2 significant digits
        pString = ['p_r = ', num2str(rankTestPValueThis, 2)];

        % Plot red if significant
        if rankTestPValueThis < sigLevel
            text(xLocThis, yLoc, pString, 'Color', 'r');
        else
            text(xLocThis, yLoc, pString, 'Color', 'k');
        end
    end
end

% Hold off if more than one column
if nColumnsToPlot > 1
    hold off
end

%% Post-plotting
% Generate a legend for the lines only if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
    if isempty(phaseVectors) && nColumnsToPlot > 1
        legend(lines, 'location', legendLocation);
    else
        legend(lines(1, :), 'location', legendLocation);
    end
end

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

% Usage: plot_tuning_curve(pValues, readout, columnsToPlot, pIsLog, pLabel, ...
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

confInts = gobjects(nColumnsToPlot, 2);

% Make the area under upperCIThis light gray
confInts(iPlot, 1) = area(pValues, upperCIThis, minY, ...
                    'LineStyle', 'none', 'FaceColor', [0.9, 0.9, 0.9]);

% Make the area under lowerCIThis white
confInts(iPlot, 2) = area(pValues, lowerCIThis, minY, ...
                    'LineStyle', 'none', 'FaceColor', [1, 1, 1]);

%                   - 'SingleColor': color when nColumnsToPlot == 1
%                   must be a 3-element vector
%                   default == rgb('SkyBlue') == [0.5273, 0.8047, 0.9180]
singleColorDefault = rgb('SkyBlue');
addParameter(iP, 'SingleColor', singleColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
singlecolor = iP.Results.SingleColor;
if nColumnsToPlot > 1
    set(lines(iPlot, 1), 'Color', colorMap(iPlot, :));
elseif nColumnsToPlot == 1
    set(lines(iPlot, 1), 'Color', singlecolor);
end

confIntColorDefault = rgb('LightGray');  
                                    % light gray confidence intervals by default

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
