function [results, handles] = virt_plot_amplitude_correlation (dataTable, plotParams, varargin)
%% Correlate and plot successive whisk amplitudes
% Usage: [results, handles] = virt_plot_amplitude_correlation (dataTable, plotParams, varargin)
% Explanation:
%       This function correlates the amplitude of successive whisks
%       (e.g., A1 vs A2, A2 vs A3) and plots these correlations in
%       separate subplots, with points colored by a grouping variable.
%
% Outputs:
%       results     - A structure containing correlation coefficients and p-values.
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
%                   - 'DataColumn': The name of the column containing the amplitude data.
%                   must be a string scalar or a character vector
%                   default == 'whiskPeakAmplitudes'
%                   - 'NCorrelations': The number of successive correlations to plot.
%                   must be a positive integer scalar
%                   default == 4
%                   - 'FigTitle': The title for the figure.
%                   must be a string scalar or a character vector
%                   default == 'Successive Whisk Amplitude Correlations'
%                   - 'FigName': The base file name for saving the figure.
%                   must be a string scalar or a character vector
%                   default == 'whisk_amplitude_scatter'
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
%       /Shared/Code/Adams_Functions/create_subplots.m
%       /Shared/Code/Adams_Functions/force_matrix.m
%       /Shared/Code/Adams_Functions/plot_grouped_scatter.m
%       /Shared/Code/Adams_Functions/plot_regression_line.m
%       /Shared/Code/Adams_Functions/resize_subplots_for_labels.m
%       /Shared/Code/Adams_Functions/save_all_figtypes.m
%
% Used by:
%       TODO: cd/virt_analyze_sniff_whisk.m
%       TODO: \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-02 - Created by Gemini.
%

%% Hard-coded parameters
textLocBestFit = 'topleft';     % Location for the best-fit line equation text
textLocThrOrig = 'bottomright'; % Location for the through-origin line equation text

%% Default values for optional arguments
groupingColumnDefault = 'repetitionNumber';
dataColumnDefault = 'whiskPeakAmplitudes';
nCorrelationsDefault = 4;
figTitleDefault = 'Successive Whisk Amplitude Correlations';
figNameDefault = 'whisk_amplitude_scatter';
outDirDefault = pwd;
figTypesDefault = {'png'};
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
addParameter(iP, 'NCorrelations', nCorrelationsDefault, @isnumeric);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);

% Parse the inputs
parse(iP, dataTable, plotParams, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
dataColumn = iP.Results.DataColumn;
nCorrToAnalyze = iP.Results.NCorrelations;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;

%% Preparation
% Initialize output structures
results = struct;
handles.fig = gobjects; % Initialize figure handle as invalid graphics object

% Exit if the table is empty or the data column is missing
if isempty(dataTable) || ~ismember(dataColumn, dataTable.Properties.VariableNames)
    fprintf('No data in column "%s" to plot. Skipping scatter plot!\n', dataColumn);
    return;
end

% Extract data from the table
whiskAmplitudesCell = dataTable.(dataColumn); % Get amplitude data (cell array of vectors)
groupingData = dataTable.(groupingColumn); % Get grouping data

% Convert the cell array of amplitudes into a matrix
whiskAmplitudesMatrix = force_matrix(whiskAmplitudesCell, 'CombineMethod', 'leftAdjustPad')';

% Determine the number of correlations to compute
nWhiskPeakOrders = size(whiskAmplitudesMatrix, 2); % Max number of whisks in a cycle
nCorrToAnalyze = min(nCorrToAnalyze, nWhiskPeakOrders - 1); % Can't correlate more than N-1 pairs

% Exit if there isn't enough data for at least one correlation
if nCorrToAnalyze < 1
    disp('Not enough whisk data to correlate successive amplitudes.');
    return;
end

%% Plotting
% Create a figure with subplots for each correlation pair
[fig, ax] = create_subplots(nCorrToAnalyze, 'AlwaysNew', true, ...
                        'ClearFigure', true, 'FigExpansion', [1, 1]);
figPath = fullfile(pathOutDir, figName); % Construct full path for saving the figure

% Initialize storage for correlation results
corrCoeffs = nan(nCorrToAnalyze, 1);
pValues = nan(nCorrToAnalyze, 1);

% Loop through each successive pair of whisk amplitudes
for iCorr = 1:nCorrToAnalyze
    % Extract amplitude data for the current pair (e.g., A1 and A2)
    ampCurrent = whiskAmplitudesMatrix(:, iCorr);
    ampNext = whiskAmplitudesMatrix(:, iCorr + 1);

    % Plot the amplitudes against each other, grouped by the specified column
    plot_grouped_scatter(ampCurrent, ampNext, groupingData, ...
        'PlotEllipse', false, 'LinkXY', true, 'GridOn', true, ...
        'Color', plotParams.colorScatter, 'MarkerType', plotParams.markerTypeScatter, ...
        'MarkerSize', plotParams.markerSizeScatter, 'MarkerLineWidth', plotParams.markerLineWidthScatter, ...
        'XLabel', sprintf('Amplitude of Whisk #%d (deg)', iCorr), ...
        'YLabel', sprintf('Amplitude of Whisk #%d (deg)', iCorr + 1), ...
        'FigTitle', 'suppress', 'LegendLocation', 'suppress', 'AxesHandle', ax(iCorr));

    % Compute and plot the standard linear regression line
    [~, ~, ~, regResults] = plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
        'ThroughOrigin', false, 'Color', plotParams.colorBestFit, 'LineStyle', plotParams.lineStyleBestFit, ...
        'LineWidth', plotParams.lineWidthBestFit, 'AxesHandle', ax(iCorr), ...
        'TextLocation', textLocBestFit, 'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', true);
    
    % Store the correlation results
    corrCoeffs(iCorr) = regResults.corrCoeff;
    pValues(iCorr) = regResults.pValue;

    % Compute and plot a regression line forced through the origin
    plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
        'ThroughOrigin', true, 'Color', plotParams.colorThrOrig, 'LineStyle', plotParams.lineStyleThrOrig, ...
        'LineWidth', plotParams.lineWidthThrOrig, 'AxesHandle', ax(iCorr), ...
        'TextLocation', textLocThrOrig, 'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', false);
end

% Add an overall title to the figure
resize_subplots_for_labels('FigTitle', figTitle);

% Save the figure to file if requested
if toSaveOutput
    save_all_figtypes(fig, figPath, figTypes);
end

% Store handles and results for output
handles.fig = fig;
results.corrCoeffs = corrCoeffs;
results.pValues = pValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%