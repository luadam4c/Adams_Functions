function [results, handles] = virt_plot_log_decrement_jitter (dataTable, plotParams, varargin)
%% Plot whisk logarithmic decrements as a grouped jitter plot
% Usage: [results, handles] = virt_plot_log_decrement_jitter (dataTable, plotParams, varargin)
% Explanation:
%       This function takes a table of whisk analysis data and generates a
%       grouped jitter plot of logarithmic decrements, with mean, 95% CI,
%       and statistical test results overlaid.
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
%                   default == 'repetitionNumber'
%                   - 'DataColumn': The name of the column with the log decrement data.
%                   must be a string scalar or a character vector
%                   default == 'whiskLogDecrements'
%                   - 'MaxDecrements': The maximum number of decrements to analyze.
%                   must be a positive integer scalar
%                   default == Inf
%                   - 'FigTitle': The title for the figure.
%                   must be a string scalar or a character vector
%                   default == 'Whisk Logarithmic Decrements'
%                   - 'FigName': The base file name for saving the figure.
%                   must be a string scalar or a character vector
%                   default == 'whisk_log_decrements_jitter'
%                   - 'OutDir': The output directory for saving the figure.
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigTypes': Figure type(s) for saving.
%                   default == {'png'}
%
% Requires:
%       cd/compute_combined_trace.m
%       cd/extract_fields.m
%       cd/force_matrix.m
%       cd/plot_grouped_jitter.m
%       cd/plot_test_result.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_difference.m
%       cd/vecfun.m
%
% Used by:
%       TODO: cd/virt_analyze_sniff_whisk.m
%       TODO: \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-02 - Created by Gemini.
%

%% Hard-coded parameters
yLocRelStar = 0.95;         % Relative y-location for significance stars
yLocRelPValue = 0.90;       % Relative y-location for p-value text
yLocRelRatioLabel = 0.85;   % Relative y-location for ratio text

%% Default values for optional arguments
groupingColumnDefault = 'repetitionNumber'; % Default column for grouping data points
dataColumnDefault = 'whiskLogDecrements';   % Default column containing the data
maxDecrementsDefault = Inf;                 % Default to analyzing all available decrements
figTitleDefault = 'Whisk Logarithmic Decrements'; % Default figure title
figNameDefault = 'whisk_log_decrements_jitter';   % Default file name for saving
outDirDefault = pwd;                        % Default output directory
figTypesDefault = {'png'};                  % Default figure save format

%% Set up Input Parser
iP = inputParser; % Create an input parser object
iP.FunctionName = mfilename; % Set the function name for error messages
addRequired(iP, 'dataTable', @istable); % The input data must be a table
addRequired(iP, 'plotParams', @isstruct); % Plotting parameters must be a structure
addParameter(iP, 'GroupingColumn', groupingColumnDefault, @ischar); % Add optional grouping column
addParameter(iP, 'DataColumn', dataColumnDefault, @ischar); % Add optional data column
addParameter(iP, 'MaxDecrements', maxDecrementsDefault, @isnumeric); % Add optional max decrements
addParameter(iP, 'FigTitle', figTitleDefault, @ischar); % Add optional figure title
addParameter(iP, 'FigName', figNameDefault, @ischar); % Add optional figure name
addParameter(iP, 'OutDir', outDirDefault, @ischar); % Add optional output directory
addParameter(iP, 'FigTypes', figTypesDefault); % Add optional figure types

% Parse the inputs
parse(iP, dataTable, plotParams, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
dataColumn = iP.Results.DataColumn;
maxDecrementsToAnalyze = iP.Results.MaxDecrements;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Initialize output structures
results = struct;
handles = struct;

% Exit early if the table is empty or the required data column doesn't exist
if isempty(dataTable) || ~ismember(dataColumn, dataTable.Properties.VariableNames)
    fprintf('No data in column "%s" to plot. Skipping jitter plot!\n', dataColumn); % Inform user
    return; % Stop execution
end

% Extract the necessary columns from the input table
logDecrementsCell = dataTable.(dataColumn); % Get log decrement data (cell array of vectors)
groupingData = dataTable.(groupingColumn); % Get grouping data (e.g., repetition number)

% Exit if there is no data to average
if isempty(logDecrementsCell) || all(cellfun(@isempty, logDecrementsCell))
    fprintf('No log decrement data found to plot.\n'); % Inform user
    return; % Stop execution
end

% Convert the cell array of log decrements into a matrix
allLogDecrementsMatrix = force_matrix(logDecrementsCell, 'CombineMethod', 'leftAdjustPad')';

% Restrict the matrix to the maximum number of decrements to analyze
nDecrementsTotal = size(allLogDecrementsMatrix, 2); % Get total number of decrement orders
nDecrementsToAnalyze = min(nDecrementsTotal, maxDecrementsToAnalyze); % Determine how many to plot
allLogDecrementsMatrix = allLogDecrementsMatrix(:, 1:nDecrementsToAnalyze); % Trim matrix

% Create matrices for decrement order and grouping, matching the data matrix size
nAnalysisWindows = size(allLogDecrementsMatrix, 1); % Number of rows (trials/cycles)
allDecrementOrdersMatrix = repmat(1:nDecrementsToAnalyze, nAnalysisWindows, 1); % Matrix of decrement orders (1, 2, 3...)
allGroupsMatrix = repmat(groupingData, 1, nDecrementsToAnalyze); % Matrix of group IDs

%% Compute Statistics
% Calculate mean and 95% confidence intervals for each decrement order
[meanLogDecrements, ~, lower95, upper95] = compute_stats_for_cellnumeric(allLogDecrementsMatrix');

% Perform significance tests (e.g., t-test) for each decrement order against zero
stats = vecfun(@test_difference, allLogDecrementsMatrix);
pValues = extract_fields(stats, 'pValue'); % Extract p-values
testFunctions = extract_fields(stats, 'testFunction'); % Extract name of statistical test used
symbols = extract_fields(stats, 'symbol'); % Extract significance symbols (*, **, etc.)

% Calculate the geometric mean of amplitude ratios from the mean log decrements
avgWhiskAmpRatios = exp(meanLogDecrements);

%% Save Results
% Store computed statistics and formatted data in the results structure
results.allLogDecrements = allLogDecrementsMatrix;
results.meanLogDecrements = meanLogDecrements;
results.lower95 = lower95;
results.upper95 = upper95;
results.pValues = pValues;
results.avgWhiskAmpRatios = avgWhiskAmpRatios;

%% Plotting
% Prepare data vectors for the plot_grouped_jitter function
allLogDecrementsVec = allLogDecrementsMatrix(:); % Flatten data matrix into a single vector
allDecrementOrdersVec = allDecrementOrdersMatrix(:); % Flatten order matrix
allGroupsVec = allGroupsMatrix(:); % Flatten group matrix

% Create x-axis tick labels (e.g., "ln(A2/A1)")
xTickLabels = arrayfun(@(x) sprintf('ln(A%d/A%d)', x+1, x), 1:nDecrementsToAnalyze, 'UniformOutput', false);

% Set up figure properties
handles.fig = set_figure_properties('AlwaysNew', true, 'ClearFigure', true);
figPath = fullfile(pathOutDir, figName); % Construct full path for saving the figure

% Generate the base jitter plot without its own statistics
plot_grouped_jitter(allLogDecrementsVec, allGroupsVec, allDecrementOrdersVec, ...
    'XTickLabels', xTickLabels, 'YLabel', 'Log Decrement (ln(A_{n+1}/A_{n}))', ...
    'LegendLocation', 'suppress', 'PlotMeanValues', false, 'PlotErrorBars', false, ...
    'RunTTest', false, 'RunRankTest', false, 'MarkerSize', plotParams.markerSizeJitter);

hold on; % Hold the current axes to overlay statistical information

% Plot the mean and 95% confidence interval error bars
xValues = (1:nDecrementsToAnalyze)'; % X-coordinates for means/error bars
errLower = meanLogDecrements - lower95; % Lower error bar length
errUpper = upper95 - meanLogDecrements; % Upper error bar length
errorbar(xValues, meanLogDecrements, errLower, errUpper, '_', ...
         'Color', 'k', 'LineWidth', 2.5, 'CapSize', 20, 'Marker', 'none');
plot(xValues, meanLogDecrements, 'o', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.5);

% Add a horizontal line at y=0 to represent the null hypothesis (no change)
yline(0, '--k', 'LineWidth', 1);

% Annotate plot with statistical significance symbols and p-values
plot_test_result(pValues, 'TestFunction', testFunctions, 'Symbol', symbols, ...
                 'XLocText', xValues, 'YLocTextRel', yLocRelPValue, 'YLocStarRel', yLocRelStar);

% Add text labels showing the geometric mean of amplitude ratios
yLimits = ylim; % Get current y-axis limits to position the text
yPosText = yLimits(1) + yLocRelRatioLabel * diff(yLimits); % Calculate y-position for the text
ratioLabels = arrayfun(@(x) sprintf('Ratio: %.2f', x), avgWhiskAmpRatios, 'UniformOutput', false); % Create labels
text(xValues, repmat(yPosText, size(xValues)), ratioLabels, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Add text to plot

% Finalize plot
title(figTitle); % Set the figure title
grid on; % Turn on the grid
hold off; % Release the axes hold

% Save the figure to file
save_all_figtypes(handles.fig, figPath, figTypes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [means, stderrs, lower95s, upper95s] = compute_stats_for_cellnumeric(vecs)
% Helper function to compute basic stats on a cell array of numeric vectors

means = compute_combined_trace(vecs, 'mean');
stderrs = compute_combined_trace(vecs, 'stderr');
lower95s = compute_combined_trace(vecs, 'lower95');
upper95s = compute_combined_trace(vecs, 'upper95');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%