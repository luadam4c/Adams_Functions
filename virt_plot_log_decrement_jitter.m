function [results, handles] = virt_plot_log_decrement_jitter (dataTable, plotParams, varargin)
%% Plot whisk logarithmic decrements as a grouped jitter plot
% Usage: [results, handles] = virt_plot_log_decrement_jitter (dataTable, plotParams, varargin)
% Explanation:
%       This function takes a table of whisk analysis data and generates a
%       grouped jitter plot of logarithmic decrements, with mean, 95% CI,
%       and statistical test results overlaid. It can create a new figure
%       or update an existing one if a handles structure is provided.
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
%                   - 'Handles': A structure of graphics handles to update.
%                   must be a structure
%                   default == []
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
%                   - 'ToSaveOutput': Whether to save the figure.
%                   must be a logical scalar
%                   default == true
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
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-02 Created by Gemini by pulling code from virt_analyze_sniff_whisk.m
% 2025-10-06 Modified by Gemini to accept an axes handle and return detailed plot handles.
% 2025-10-06 Modified by Gemini to include update logic from virt_plot_whisk_analysis.m
% 2025-10-06 Fixed by Gemini to handle more than one group.
%

%% Hard-coded parameters
yLocStarRel = 0.95;         % Relative y-location for significance stars
yLocPValueRel = 0.90;       % Relative y-location for p-value text
yLocRatioLabelRel = 0.85;   % Relative y-location for ratio text

%% Default values for optional arguments
groupingColumnDefault = 'repetitionNumber'; % Default column for grouping data points
dataColumnDefault = 'whiskLogDecrements';   % Default column containing the data
maxDecrementsDefault = Inf;                 % Default to analyzing all available decrements
handlesDefault = [];                        % Default is to create a new plot
figTitleDefault = 'Whisk Logarithmic Decrements'; % Default figure title
figNameDefault = 'whisk_log_decrements_jitter';   % Default file name for saving
outDirDefault = pwd;                        % Default output directory
figTypesDefault = {'png'};                  % Default figure save format
toSaveOutputDefault = true;                 % Default to save the output figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addRequired(iP, 'dataTable', @istable);
addRequired(iP, 'plotParams', @isstruct);
addParameter(iP, 'GroupingColumn', groupingColumnDefault, @ischar);
addParameter(iP, 'DataColumn', dataColumnDefault, @ischar);
addParameter(iP, 'MaxDecrements', maxDecrementsDefault, @isnumeric);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);

% Parse the inputs
parse(iP, dataTable, plotParams, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
dataColumn = iP.Results.DataColumn;
maxDecrementsToAnalyze = iP.Results.MaxDecrements;
handlesIn = iP.Results.Handles;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;

%% Preparation
% Initialize output structures
results = struct;
handles = struct;

% Exit early if the table is empty or the required data column doesn't exist
if isempty(dataTable) || ~ismember(dataColumn, dataTable.Properties.VariableNames)
    fprintf('No data in column "%s" to plot. Skipping jitter plot!\n', dataColumn);
    return;
end

% Extract the necessary columns from the input table
logDecrementsCell = dataTable.(dataColumn); % Get log decrement data (cell array of vectors)
groupingData = dataTable.(groupingColumn); % Get grouping data (e.g., repetition number)

% Exit if there is no data to average
if isempty(logDecrementsCell) || all(cellfun(@isempty, logDecrementsCell))
    fprintf('No log decrement data found to plot.\n');
    return;
end

% Convert the cell array of log decrements into a matrix
allLogDecrementsMatrix = force_matrix(logDecrementsCell, 'CombineMethod', 'leftAdjustPad')';

% Restrict the matrix to the maximum number of decrements to analyze
nDecrementsTotal = size(allLogDecrementsMatrix, 2); % Get total number of decrement orders
nDecrementsToAnalyze = min(nDecrementsTotal, maxDecrementsToAnalyze); % Determine how many to plot
allLogDecrementsMatrix = allLogDecrementsMatrix(:, 1:nDecrementsToAnalyze); % Trim matrix

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
% Extract plotting parameters
jitterWidth = plotParams.jitterWidth;
markerSizeJitter = plotParams.markerSizeJitter;

% Check if we are updating an existing plot or creating a new one
if ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig)
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn; % Use the passed-in handles struct
    axJitter = handles.axJitter;
    
    nAnalysisWindows = size(allLogDecrementsMatrix, 1);
    decrementOrdersMatrix = repmat(1:nDecrementsToAnalyze, nAnalysisWindows, 1);
    decrementOrdersVec = decrementOrdersMatrix(:);
    xValues = (1:nDecrementsToAnalyze)';
    errLower = meanLogDecrements - lower95;
    errUpper = upper95 - meanLogDecrements;

    % Update x limits
    xlim(axJitter, [min(decrementOrdersVec) - 0.5, max(decrementOrdersVec) + 0.5]);

    % Update jitter data for each group separately
    allLogDecrementsVec = allLogDecrementsMatrix(:);
    allGroupsVec = repmat(groupingData, 1, nDecrementsToAnalyze);
    allGroupsVec = allGroupsVec(:);

    uniqueGroups = unique(groupingData);
    nGroups = numel(uniqueGroups);
    currentJitterHandles = handles.hJitter;

    if nGroups == numel(currentJitterHandles)
        for iGroup = 1:nGroups
            groupValue = uniqueGroups(iGroup);
            groupMask = (allGroupsVec == groupValue);
            
            xDataGroup = decrementOrdersVec(groupMask);
            yDataGroup = allLogDecrementsVec(groupMask);
            nPointsGroup = sum(groupMask);
            
            xDataJittered = xDataGroup + (jitterWidth * (rand(nPointsGroup, 1) - 0.5));
            
            set(currentJitterHandles(iGroup), 'XData', xDataJittered, 'YData', yDataGroup);
        end
    else
        % Fallback for safety if number of groups changes
        delete(currentJitterHandles);
        hOutJitter = plot_grouped_jitter(allLogDecrementsVec, allGroupsVec, decrementOrdersVec, ...
            'AxesHandle', axJitter, 'UsePlotSpread', false, 'JitterWidth', jitterWidth, ...
            'PlotMeanValues', false, 'PlotErrorBars', false, 'RunTTest', false, ...
            'RunRankTest', false, 'MarkerSize', markerSizeJitter, ...
            'LegendLocation', 'suppress', 'XTickLabels', get(axJitter, 'XTickLabel'));
        handles.hJitter = hOutJitter.distributions;
    end

    % Update means and error bars
    handles.hErrorBars.XData = xValues';
    handles.hErrorBars.YData = meanLogDecrements';
    handles.hErrorBars.YNegativeDelta = errLower';
    handles.hErrorBars.YPositiveDelta = errUpper';
    handles.hMeans.XData = xValues';
    handles.hMeans.YData = meanLogDecrements';

    % Remove old annotations that need to be replotted
    delete(handles.pTextJitter);
    delete(handles.sigMarkerJitter);
    delete(handles.ratioLabels);

    % Update p value texts and significance markers
    hOut = plot_test_result(pValues, 'TestFunction', testFunctions, 'Symbol', symbols, ...        
                         'XLocText', xValues, 'YLocTextRel', yLocPValueRel, ...
                         'YLocStarRel', yLocStarRel, 'AxesHandle', axJitter);
    handles.pTextJitter = hOut.pText;
    handles.sigMarkerJitter = hOut.sigMarker;

    % Update ratio labels
    yLimits = ylim(axJitter);
    yPosText = yLimits(1) + yLocRatioLabelRel * diff(yLimits);
    ratioLabels = arrayfun(@(x) sprintf('Ratio: %.2f', x), avgWhiskAmpRatios, 'UniformOutput', false);
    handles.ratioLabels = text(axJitter, xValues, repmat(yPosText, size(xValues)), ratioLabels, ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

    fig = handles.fig; % Get figure handle for saving
else
    % --- CREATE NEW PLOT ---
    % Prepare data vectors for the plot_grouped_jitter function
    allLogDecrementsVec = allLogDecrementsMatrix(:); % Flatten data matrix into a single vector
    allDecrementOrdersVec = repmat(1:nDecrementsToAnalyze, size(allLogDecrementsMatrix, 1), 1);
    allDecrementOrdersVec = allDecrementOrdersVec(:); % Flatten order matrix
    allGroupsVec = repmat(groupingData, 1, nDecrementsToAnalyze); % Matrix of group IDs
    allGroupsVec = allGroupsVec(:); % Flatten group matrix

    % Create x-axis tick labels (e.g., "ln(A2/A1)")
    xTickLabels = arrayfun(@(x) sprintf('ln(A%d/A%d)', x+1, x), 1:nDecrementsToAnalyze, 'UniformOutput', false);

    % Set up figure and axes
    fig = set_figure_properties('AlwaysNew', true, 'ClearFigure', true);
    axJitter = gca;
    figPath = fullfile(pathOutDir, figName); % Construct full path for saving the figure

    % Generate the base jitter plot without its own statistics
    hOutJitter = plot_grouped_jitter(allLogDecrementsVec, allGroupsVec, allDecrementOrdersVec, ...
        'AxesHandle', axJitter, 'UsePlotSpread', false, 'JitterWidth', jitterWidth, ...
        'XTickLabels', xTickLabels, 'YLabel', 'Log Decrement (ln(A_{n+1}/A_{n}))', ...
        'LegendLocation', 'suppress', 'PlotMeanValues', false, 'PlotErrorBars', false, ...
        'RunTTest', false, 'RunRankTest', false, 'MarkerSize', markerSizeJitter);

    hold on; % Hold the current axes to overlay statistical information

    % Plot the mean and 95% confidence interval error bars
    xValues = (1:nDecrementsToAnalyze)'; % X-coordinates for means/error bars
    errLower = meanLogDecrements - lower95; % Lower error bar length
    errUpper = upper95 - meanLogDecrements; % Upper error bar length
    hErrorBars = errorbar(xValues, meanLogDecrements, errLower, errUpper, '_', ...
             'Color', 'k', 'LineWidth', 2.5, 'CapSize', 20, 'Marker', 'none');
    hMeans = plot(xValues, meanLogDecrements, 'o', 'MarkerEdgeColor', 'k', ...
         'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.5);

    % Add a horizontal line at y=0 to represent the null hypothesis (no change)
    hNull = yline(0, '--k', 'LineWidth', 1);

    % Annotate plot with statistical significance symbols and p-values
    hOutTest = plot_test_result(pValues, 'TestFunction', testFunctions, 'Symbol', symbols, ...
                     'XLocText', xValues, 'YLocTextRel', yLocPValueRel, 'YLocStarRel', yLocStarRel, 'AxesHandle', axJitter);

    % Add text labels showing the geometric mean of amplitude ratios
    yLimits = ylim; % Get current y-axis limits to position the text
    yPosText = yLimits(1) + yLocRatioLabelRel * diff(yLimits); % Calculate y-position for the text
    ratioLabels = arrayfun(@(x) sprintf('Ratio: %.2f', x), avgWhiskAmpRatios, 'UniformOutput', false); % Create labels
    hRatioLabels = text(xValues, repmat(yPosText, size(xValues)), ratioLabels, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % Add text to plot

    % Finalize plot
    title(figTitle); % Set the figure title
    grid on; % Turn on the grid
    hold off; % Release the axes hold

    % Store handles for output
    handles.fig = fig;
    handles.axJitter = axJitter;
    handles.hJitter = hOutJitter.distributions;
    handles.hErrorBars = hErrorBars;
    handles.hMeans = hMeans;
    handles.hNull = hNull;
    handles.pTextJitter = hOutTest.pText;
    handles.sigMarkerJitter = hOutTest.sigMarker;
    handles.ratioLabels = hRatioLabels;
end

% Save the figure to file if requested
if toSaveOutput
    figPath = fullfile(pathOutDir, figName); % Construct full path for saving
    save_all_figtypes(fig, figPath, figTypes);
end

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
