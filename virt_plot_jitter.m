function [results, handles] = virt_plot_jitter (dataTable, plotParams, varargin)
%% Plot whisk logarithmic decrements or correlation Z-scores as a grouped jitter plot
% Usage: [results, handles] = virt_plot_jitter (dataTable, plotParams, varargin)
% Explanation:
%       This function takes a table of whisk analysis data and generates a
%       grouped jitter plot of either logarithmic decrements or Fisher
%       Z-scores.
%       - Parametric: Plots Mean and 95% Confidence Interval.
%       - Nonparametric: Plots Median, 95% CI of Median (Notches), and IQR (Box).
%       It can create a new figure or update an existing one if a handles
%       structure is provided.
%
% Outputs:
%       results     - A structure containing computed statistics and data used for plotting.
%                   specified as a structure
%       handles     - A structure containing handles to the generated plot objects.
%                   specified as a structure
%
% Arguments:
%       dataTable   - A table containing whisk analysis data.
%                   must be a table
%       plotParams  - The plotting parameters structure (e.g., from P.Plotting).
%                   must be a structure
%       varargin    - 'GroupingColumn': The name of the column to group data by.
%                   must be a string scalar or a character vector
%                   default == 'seedNumber' or 'fileNumber'
%                   - 'DataColumn': The name of the column with the log decrement data.
%                   must be a string scalar or a character vector
%                   default == 'whiskLogDecrements'
%                   - 'DataMode': The type of data being plotted.
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'LogDecrement' - Data is logarithmic decrements.
%                       'FisherZScore' - Data is Fisher-transformed Z-scores.
%                   default == 'LogDecrement'
%                   - 'MaxOrders': The maximum number of decrement/correlation orders to analyze.
%                   must be a positive integer scalar
%                   default == Inf
%                   - 'Handles': A structure of graphics handles to update.
%                   must be a structure
%                   default == []
%                   - 'FigTitle': The title for the figure.
%                   must be a string scalar or a character vector
%                   default == set based on DataMode
%                   - 'FigName': The base file name for saving the figure.
%                   must be a string scalar or a character vector
%                   default == set based on DataMode
%                   - 'OutDir': The output directory for saving the figure.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigTypes': Figure type(s) for saving.
%                   default == {'png'}
%                   - 'ToSaveOutput': Whether to save the figure.
%                   must be a logical scalar
%                   default == true
%                   - 'ShowFigure': whether to show figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TestType': type of test to perform
%                   must be an unambiguous, case-insensitive match to one of:
%                       'auto'          - decide based on normality of all columns
%                       'parametric'    - force parametric test
%                       'nonparametric' - force nonparametric test
%                   default == 'auto'
%
% Requires:
%       cd/compute_combined_trace.m
%       cd/extract_fields.m
%       cd/force_matrix.m
%       cd/plot_grouped_jitter.m
%       cd/plot_test_result.m
%       cd/save_all_figtypes.m
%       cd/create_subplots.m
%       cd/test_difference.m
%       cd/test_normality.m
%       cd/vecfun.m
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m

% File History:
% 2025-10-02 Created by Gemini by pulling code from virt_analyze_sniff_whisk.m
% 2025-10-06 Modified by Gemini to accept an axes handle and return detailed plot handles.
% 2025-10-06 Modified by Gemini to include update logic from virt_plot_whisk_analysis.m
% 2025-10-06 Fixed by Gemini to handle more than one group.
% 2025-10-15 Renamed to virt_plot_jitter.m by Gemini.
% 2025-10-15 Added 'DataMode' to handle Fisher Z-scores by Gemini.
% 2025-10-17 Added 'ShowFigure' as an optional argument.
% 2026-01-09 Added 'TestType' optional argument.
% 2026-01-09 Modified by Gemini to plot notched boxplots for nonparametric data.
% 2026-01-09 Modified by Gemini to implement weighted global auto-detection for TestType.
% 2026-01-17 Now uses create_subplots.m instead of set_figure_properties.m
% 2026-01-17 Fixed showFigure issues by removing axes()
% 2026-01-19 Set colorJitter to black and colorSummary to red

%% Hard-coded parameters
yLocStarRel = 0.95;         % Relative y-location for significance stars
yLocPValueRel = 0.90;       % Relative y-location for p-value text
yLocTransformedLabelRel = 0.85; % Relative y-location for ratio/corr text
validDataModes = {'LogDecrement', 'FisherZScore'};
validTestTypes = {'auto', 'parametric', 'nonparametric'};
colorJitter = 'k';
colorSummary = 'r';

%% Default values for optional arguments
groupingColumnDefault = [];     % set later
dataColumnDefault = '';         % set later
dataModeDefault = 'LogDecrement';
maxOrdersDefault = Inf;         % Default to analyzing all available
handlesDefault = [];            % Default is to create a new plot
figTitleDefault = '';           % Default figure title is set based on DataMode
figNameDefault = '';            % Default file name is set based on DataMode
outDirDefault = pwd;            % Default output directory
figTypesDefault = {'png'};      % Default figure save format
toSaveOutputDefault = true;     % Default to save the output figure
showFigureDefault = true;       % Default to show the figure
testTypeDefault = 'auto';       % Default to auto-detect test type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addRequired(iP, 'dataTable', @istable);
addRequired(iP, 'plotParams', @isstruct);
addParameter(iP, 'GroupingColumn', groupingColumnDefault);
addParameter(iP, 'DataColumn', dataColumnDefault, @ischar);
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) any(validatestring(x, validDataModes)));
addParameter(iP, 'MaxOrders', maxOrdersDefault, @isnumeric);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);
addParameter(iP, 'ShowFigure', showFigureDefault, @islogical);
addParameter(iP, 'TestType', testTypeDefault, ...
    @(x) any(validatestring(x, validTestTypes)));

% Parse the inputs
parse(iP, dataTable, plotParams, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
dataColumn = iP.Results.DataColumn;
dataMode = validatestring(iP.Results.DataMode, validDataModes);
maxOrdersToAnalyze = iP.Results.MaxOrders;
handlesIn = iP.Results.Handles;
figTitleUser = iP.Results.FigTitle;
figNameUser = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;
showFigure = iP.Results.ShowFigure;
testType = validatestring(iP.Results.TestType, validTestTypes);

%% Preparation
% Initialize output structures
results = struct;
handles = struct;

% Decide on data column
if isempty(dataColumn)
    switch dataMode
    case 'LogDecrement'
        dataColumn = 'whiskLogDecrements';
    case 'FisherZScore'
        dataColumn = 'whiskAmpCorrZScoresBasal';
    end
end

% Configure plot settings based on DataMode
switch dataMode
case 'LogDecrement'
    yLabel = 'Log Decrement (ln(A_{n+1}/A_{n}))';
    xTickLabelFormat = 'ln(A%d/A%d)';
    transformFunc = @exp;
    transformLabel = 'Ratio';
    figTitleDefault = 'Whisk Logarithmic Decrements';
    figNameDefault = 'whisk_log_decrements_jitter';
case 'FisherZScore'
    yLabel = 'Fisher Z-Score';
    xTickLabelFormat = 'A%d-A%d';
    transformFunc = @tanh;
    transformLabel = 'Corr';
    figTitleDefault = 'Whisk Amplitude Correlation Z-Scores';
    figNameDefault = 'whisk_corr_zscores_jitter';
end

% Use user-provided titles/names if available, otherwise use defaults
figTitle = figTitleUser;
if isempty(figTitle)
    figTitle = figTitleDefault;
end
figName = figNameUser;
if isempty(figName)
    figName = figNameDefault;
end

% Exit early if the table is empty or the required data column doesn't exist
if isempty(dataTable) || ~ismember(dataColumn, dataTable.Properties.VariableNames)
    fprintf('No data in column "%s" to plot. Skipping jitter plot!\n', dataColumn);
    return;
end

% Get the column names
columnNames = dataTable.Properties.VariableNames;

% Decide on grouping column
if isempty(groupingColumn)
    if ismember('seedNumber', columnNames)
        groupingColumn = 'seedNumber';
    elseif ismember('fileNumber', columnNames)
        groupingColumn = 'fileNumber';
    elseif ismember('repetitionNumber', columnNames)
        groupingColumn = 'repetitionNumber';
    else
        % Create dummy column
        dataTable.grouping = ones(height(dataTable), 1);
        groupingColumn = 'grouping';
    end
end

% Extract the necessary columns from the input table
dataCell = dataTable.(dataColumn); % Get data (cell array of vectors)
groupingData = dataTable.(groupingColumn); % Get grouping data (e.g., seed number)

% Exit if there is no data to average
if isempty(dataCell) || all(cellfun(@isempty, dataCell))
    fprintf('No data found to plot.\n');
    return;
end

% Convert the cell array of log decrements or z scores into a matrix
%   Note: each column is a decrement/correlation order
allDataMatrix = force_matrix(dataCell, 'CombineMethod', 'leftAdjustPad')';

%% Restrict the matrix to the maximum number of decrements to analyze
% Get total number of decrement/correlation orders
nOrdersTotal = size(allDataMatrix, 2);

% Determine how many to plot
nOrdersToAnalyze = min(nOrdersTotal, maxOrdersToAnalyze);

% Trim matrix
allDataMatrix = allDataMatrix(:, 1:nOrdersToAnalyze);

%% Compute Statistics
% Calculate mean, 95% CI, median and IQR for each decrement order
[meanData, ~, lower95, upper95, medians, ...
        quartile25s, quartile75s, lower95med, upper95med] = ...
    compute_stats_for_cellnumeric(allDataMatrix');

% Determine the global test type if 'auto' is selected
if strcmpi(testType, 'auto')
    % Check normality for each column (decrement order)
    [isNormalCols, ~] = test_normality(allDataMatrix);
    
    % Define weights: more weight given to earlier columns (harmonic weights: 1, 1/2, 1/3...)
    nCols = numel(isNormalCols);
    weights = 1 ./ (1:nCols);
    
    % Calculate weighted fraction of normal columns
    % Ensure row vectors for dot product
    weightedFraction = sum(double(isNormalCols(:))' .* weights(:)') / sum(weights);
    
    % Set test type based on weighted majority
    if weightedFraction > 0.5
        testType = 'parametric';
    else
        testType = 'nonparametric';
    end
end

% Perform significance tests (e.g., t-test or ranksum) for each decrement order
stats = vecfun(@(x) test_difference(x, 'TestType', testType), allDataMatrix);
pValues = extract_fields(stats, 'pValue'); % Extract p-values
testFunctions = extract_fields(stats, 'testFunction'); % Extract name of statistical test used
symbols = extract_fields(stats, 'symbol'); % Extract significance symbols (*, **, etc.)
useParametric = extract_fields(stats, 'useParametric'); % Extract whether parametric test was used

% Determine which stats to plot based on the test type used
% For parametric: Mean
% For nonparametric: Median
plotData = NaN(size(meanData));
for iOrder = 1:numel(useParametric)
    if useParametric(iOrder)
        plotData(iOrder) = meanData(iOrder);
    else
        plotData(iOrder) = medians(iOrder);
    end
end

% Apply the inverse transformation (exp for log, tanh for Z-score)
avgTransformedValues = transformFunc(plotData);

%% Save Results
% Store computed statistics and formatted data in the results structure
switch dataMode
case 'LogDecrement'
    results.allLogDecrements = allDataMatrix;
    results.meanLogDecrements = meanData;
    results.lower95LogDecrements = lower95;
    results.upper95LogDecrements = upper95;
    results.medianLogDecrements = medians;
    results.firstQuartileLogDecrements = quartile25s;
    results.thirdQuartileLogDecrements = quartile75s;
    results.lower95MedianLogDecrements = lower95med;
    results.upper95MedianLogDecrements = upper95med;
    results.pValuesLogDecrements = pValues;
    results.avgWhiskAmpRatios = avgTransformedValues;
case 'FisherZScore'
    results.allCorrZScores = allDataMatrix;
    results.meanCorrZScores = meanData;
    results.lower95CorrZScores = lower95;
    results.upper95CorrZScores = upper95;
    results.medianCorrZScores = medians;
    results.firstQuartileCorrZScores = quartile25s;
    results.thirdQuartileCorrZScores = quartile75s;
    results.lower95MedianCorrZScores = lower95med;
    results.upper95MedianCorrZScores = upper95med;
    results.pValuesCorrZScores = pValues;
    results.avgWhiskAmpCorrs = avgTransformedValues;
end

%% Prepare data vectors for the plot_grouped_jitter function
% Flatten data matrix into a column vector
allDataVec = allDataMatrix(:);

% Create decrement order matrix
allOrdersVec = repmat(1:nOrdersToAnalyze, size(allDataMatrix, 1), 1);

% Flatten decrement order matrix into a column vector
allOrdersVec = allOrdersVec(:);

% Create grouping matrix
allGroupsVec = repmat(groupingData, 1, nOrdersToAnalyze);

% Flatten group matrix into a column vector
allGroupsVec = allGroupsVec(:);

%% Plotting
% Extract plotting parameters
jitterWidth = plotParams.jitterWidth;
markerSizeJitter = plotParams.markerSizeJitter;

% Determine x-coordinates for plots
xValues = (1:nOrdersToAnalyze)';

% Check if we are updating an existing plot or creating a new one
if ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig)
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn; % Use the passed-in handles struct
    axJitter = handles.axJitter;
    
    % Prepare data vectors for the plot_grouped_jitter function
    nAnalysisWindows = size(allDataMatrix, 1);
    ordersMatrix = repmat(1:nOrdersToAnalyze, nAnalysisWindows, 1);
    ordersVec = ordersMatrix(:);

    % Update x limits
    xlim(axJitter, [min(ordersVec) - 0.5, max(ordersVec) + 0.5]);

    % Compute the number of unique groups
    uniqueGroups = unique(groupingData);
    nGroups = numel(uniqueGroups);

    % Get the current jitter plot handles
    currentJitterHandles = handles.hJitter;

    % Update or plot jitter plots
    if nGroups == numel(currentJitterHandles)
        for iGroup = 1:nGroups
            groupValue = uniqueGroups(iGroup);
            groupMask = (allGroupsVec == groupValue);
            
            xDataGroup = ordersVec(groupMask);
            yDataGroup = allDataVec(groupMask);
            nPointsGroup = sum(groupMask);
            
            xDataJittered = xDataGroup + (jitterWidth * (rand(nPointsGroup, 1) - 0.5));
            
            set(currentJitterHandles(iGroup), 'XData', xDataJittered, 'YData', yDataGroup);
        end
    else
        % Fallback for safety if number of groups changes
        delete(currentJitterHandles);
        hOutJitter = plot_grouped_jitter(allDataVec, allGroupsVec, ordersVec, ...
            'AxesHandle', axJitter, 'UsePlotSpread', false, 'JitterWidth', jitterWidth, ...
            'PlotMeanValues', false, 'PlotErrorBars', false, 'RunTTest', false, ...
            'RunRankTest', false, 'MarkerSize', markerSizeJitter, ...
            'ColorMap', colorJitter, ...
            'LegendLocation', 'suppress', 'XTickLabels', get(axJitter, 'XTickLabel'));
        handles.hJitter = hOutJitter.distributions;
    end
    
    % --- Remove old stats objects to redraw them ---
    if isfield(handles, 'hErrorBars'), delete(handles.hErrorBars); end
    if isfield(handles, 'hMeans'), delete(handles.hMeans); end
    if isfield(handles, 'hBoxPlots'), delete(handles.hBoxPlots); end
    
    % --- Remove old annotations that need to be replotted ---
    delete(handles.pTextJitter);
    delete(handles.sigMarkerJitter);
    delete(handles.transformedLabels);

    % Get figure handle for saving
    fig = handles.fig; 
else
    % --- CREATE NEW PLOT ---
    % Create x-axis tick labels
    xTickLabels = arrayfun(@(x) sprintf(xTickLabelFormat, x+1, x), 1:nOrdersToAnalyze, 'UniformOutput', false);

    % Set up figure and axes
    [fig, axJitter] = create_subplots(1, 1, 'AlwaysNew', true, 'ShowFigure', showFigure);

    % Generate the base jitter plot without its own statistics
    hOutJitter = plot_grouped_jitter(allDataVec, allGroupsVec, allOrdersVec, ...
        'AxesHandle', axJitter, 'UsePlotSpread', false, 'JitterWidth', jitterWidth, ...
        'XTickLabels', xTickLabels, 'YLabel', yLabel, ...
        'LegendLocation', 'suppress', 'PlotMeanValues', false, 'PlotErrorBars', false, ...
        'RunTTest', false, 'RunRankTest', false, 'MarkerSize', markerSizeJitter, ...
        'ColorMap', colorJitter);

    hold on; % Hold the current axes to overlay statistical information

    % Add a horizontal line at y=0 to represent the null hypothesis (no change)
    hNull = yline(0, '--k', 'LineWidth', 1);

    % Store handles
    handles.fig = fig;
    handles.axJitter = axJitter;
    handles.hJitter = hOutJitter.distributions;
    handles.hNull = hNull;
    
    % Initialize empty handles for optional stats
    handles.hErrorBars = [];
    handles.hMeans = [];
    handles.hBoxPlots = [];
end

%% Draw Summary Statistics (Shared by New and Update)
% Identify which groups are parametric vs nonparametric
idxParametric = find(useParametric);
idxNonParametric = find(~useParametric);

% 1. Plot Parametric Stats: Mean +/- 95% Confidence Interval
if ~isempty(idxParametric)
    % Extract parametric data
    xPara = xValues(idxParametric);
    yPara = meanData(idxParametric);
    lErr = meanData(idxParametric) - lower95(idxParametric);
    uErr = upper95(idxParametric) - meanData(idxParametric);
    
    hErrorBars = errorbar(axJitter, xPara, yPara, lErr, uErr, '_', ...
             'Color', colorSummary, 'LineWidth', 2.5, 'CapSize', 20, 'Marker', 'none');
    hMeans = plot(axJitter, xPara, yPara, 'o', 'MarkerEdgeColor', colorSummary, ...
            'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.5);
         
    handles.hErrorBars = hErrorBars;
    handles.hMeans = hMeans;
end

% 2. Plot Nonparametric Stats: Notched Box Plot (Median, 95% CI of Median, IQR)
if ~isempty(idxNonParametric) && size(allDataMatrix, 1) > 1
    % Extract data for nonparametric groups
    % boxplot expects a matrix where each column is a group. 
    % We use 'Positions' to place them correctly on the x-axis.
    dataNonPara = allDataMatrix(:, idxNonParametric);
    xNonPara = xValues(idxNonParametric);
    
    % Draw Notched Box Plots
    % 'Symbol','' suppresses outliers (since jitter shows them)
    % 'Widths', 0.5 matches typical bar/errorbar width
    hBoxPlots = boxplot(axJitter, dataNonPara, 'Positions', xNonPara, ...
        'Notch', 'on', 'Symbol', '', 'Colors', colorSummary, 'Widths', 0.5);
    
    % Customize Box Plot appearance (make lines thicker)
    set(hBoxPlots, 'LineWidth', 2);
    
    % Note: boxplot sometimes resets XTickLabels, so we might need to restore them if this was a new plot
    % but standard practice suggests it should be fine if Positions are set.
    
    handles.hBoxPlots = hBoxPlots;
end

%% Update Annotations
% Annotate plot with statistical significance symbols and p-values
hOutTest = plot_test_result(pValues, 'TestFunction', testFunctions, ...
                'Symbol', symbols, 'XLocText', xValues, ...
                'YLocTextRel', yLocPValueRel, 'YLocStarRel', yLocStarRel, ...
                'AxesHandle', axJitter);

% Add text labels showing the transformed values (ratios or correlations)
yLimits = ylim(axJitter); % Get current y-axis limits to position the text
yPosText = yLimits(1) + yLocTransformedLabelRel * diff(yLimits); % Calculate y-position for the text
transformedLabels = arrayfun(@(x) sprintf('%s: %.2f', transformLabel, x), avgTransformedValues, 'UniformOutput', false); % Create labels
hTransformedLabels = text(axJitter, xValues, repmat(yPosText, size(xValues)), transformedLabels, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Add text to plot

handles.pTextJitter = hOutTest.pText;
handles.sigMarkerJitter = hOutTest.sigMarker;
handles.transformedLabels = hTransformedLabels;

% Finalize plot
if isempty(handlesIn)
    title(figTitle); % Set the figure title
    grid on; % Turn on the grid
end

% Save the figure to file if requested
if toSaveOutput
    figPath = fullfile(pathOutDir, figName); % Construct full path for saving
    save_all_figtypes(fig, figPath, figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [means, stderrs, lower95s, upper95s, medians, ...
            quartile25s, quartile75s, lower95med, upper95med] = ...
    compute_stats_for_cellnumeric(vecs)
% Helper function to compute basic stats on a cell array of numeric vectors

% Compute parametric stats
means = compute_combined_trace(vecs, 'mean');
stderrs = compute_combined_trace(vecs, 'stderr');
lower95s = compute_combined_trace(vecs, 'lower95');
upper95s = compute_combined_trace(vecs, 'upper95');

% Compute nonparametric stats
medians = compute_combined_trace(vecs, 'median');
quartile25s = compute_combined_trace(vecs, 'quartile25');
quartile75s = compute_combined_trace(vecs, 'quartile75');
lower95med = compute_combined_trace(vecs, 'lower95med');
upper95med = compute_combined_trace(vecs, 'upper95med');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%