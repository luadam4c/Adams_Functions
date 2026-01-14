function [results, handles] = virt_plot_amplitude_correlation (dataTable, plotParams, varargin)
%% Correlate and plot successive whisk amplitudes
% Usage: [results, handles] = virt_plot_amplitude_correlation (dataTable, plotParams, varargin)
% Explanation:
%       This function correlates the amplitude of successive whisks
%       (e.g., A1 vs A2, A2 vs A3) and plots these correlations in
%       separate subplots, with points colored by a grouping variable.
%       It can create a new figure or update an existing one if a handles 
%       structure is provided.
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
%                   default == 'seedNumber' or 'fileNumber'
%                   - 'DataColumn': The name of the column containing the amplitude data.
%                   must be a string scalar or a character vector
%                   default == 'whiskPeakAmplitudes'
%                   - 'NCorrelations': The number of successive correlations to plot.
%                   must be a positive integer scalar
%                   default == 4
%                   - 'Handles': A structure of graphics handles to update.
%                   must be a structure
%                   default == []
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
%                   - 'ShowFigure': whether to show figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/create_subplots.m
%       cd/force_matrix.m
%       cd/plot_grouped_scatter.m
%       cd/plot_regression_line.m
%       cd/resize_subplots_for_labels.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-02 Created by Gemini by pulling code from virt_analyze_sniff_whisk.m
% 2025-10-06 Modified by Gemini to accept axes handles and return detailed plot handles.
% 2025-10-06 Modified by Gemini to include update logic from virt_plot_whisk_analysis.m
% 2025-10-06 Fixed by Gemini to handle more than one group.
% 2025-10-17 Added 'ShowFigure' as an optional argument.
% 2026-01-14 Fixed "object deleted" error by validating handles before update (Gemini)

%% Hard-coded parameters
textLocBestFitDefault = 'topleft';     % Location for the best-fit line equation text
textLocThrOrigDefault = 'bottomright'; % Location for the through-origin line equation text
axisCoveragePercScatter = 90;          % Coverage for compute_axis_limits

%% Default values for optional arguments
groupingColumnDefault = [];                 % set later
dataColumnDefault = 'whiskPeakAmplitudes';
nCorrelationsDefault = 4;
handlesDefault = [];
figTitleDefault = 'Successive Whisk Amplitude Correlations';
figNameDefault = 'whisk_amplitude_scatter';
outDirDefault = pwd;
figTypesDefault = {'png'};
toSaveOutputDefault = true;                 % Default to save the output figure
showFigureDefault = true;                   % Default to show the figure

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
addParameter(iP, 'NCorrelations', nCorrelationsDefault, @isnumeric);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);
addParameter(iP, 'ShowFigure', showFigureDefault, @islogical);

% Parse the inputs
parse(iP, dataTable, plotParams, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
dataColumn = iP.Results.DataColumn;
nCorrToAnalyze = iP.Results.NCorrelations;
handlesIn = iP.Results.Handles;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;
showFigure = iP.Results.ShowFigure;

%% Preparation
% Initialize output structures
results = struct;
handles.fig = gobjects; % Initialize figure handle as invalid graphics object

% Set text locations for regression lines, with defaults
if isfield(plotParams, 'textLocBestFit')
    textLocBestFit = plotParams.textLocBestFit;
else
    textLocBestFit = textLocBestFitDefault;
end
if isfield(plotParams, 'textLocThrOrig')
    textLocThrOrig = plotParams.textLocThrOrig;
else
    textLocThrOrig = textLocThrOrigDefault;
end

% Exit if the table is empty or the data column is missing
if isempty(dataTable) || ~ismember(dataColumn, dataTable.Properties.VariableNames)
    fprintf('No data in column "%s" to plot. Skipping scatter plot!\n', dataColumn);
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
% Check if we are updating an existing plot or creating a new one
if ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig) && ...
   isfield(handlesIn, 'hScatters') && numel(handlesIn.hScatters) == nCorrToAnalyze && ...
   isfield(handlesIn, 'axScatter') && all(isgraphics(handlesIn.axScatter))
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn;
    fig = handles.fig;

    % Delete old regression lines and text
    delete(handles.hBestFits);
    delete(handles.hTextBestFit);
    delete(handles.hThrOrigs);
    delete(handles.hTextThrOrig);

    % Run through all scatter plots to update them
    hBestFits = gobjects(nCorrToAnalyze, 1);
    hTextBestFit = gobjects(nCorrToAnalyze, 1);
    hThrOrigs = gobjects(nCorrToAnalyze, 1);
    hTextThrOrig = gobjects(nCorrToAnalyze, 1);
    for iCorr = 1:nCorrToAnalyze
        % Get the data for this pair
        ampCurrent = whiskAmplitudesMatrix(:, iCorr);
        ampNext = whiskAmplitudesMatrix(:, iCorr + 1);

        % Update scatter data for each group
        uniqueGroups = unique(groupingData);
        nGroups = numel(uniqueGroups);
        currentScatterHandles = handles.hScatters{iCorr};

        % Validate if the specific scatter objects are still valid
        areScatterHandlesValid = all(isgraphics(currentScatterHandles));

        % If handles are valid and the number of groups hasn't changed, update data
        if areScatterHandlesValid && nGroups == numel(currentScatterHandles)
            for iGroup = 1:nGroups
                groupValue = uniqueGroups(iGroup);
                groupMask = (groupingData == groupValue);
                set(currentScatterHandles(iGroup), ...
                    'XData', ampCurrent(groupMask), 'YData', ampNext(groupMask));
            end
        else
            % Fallback: Re-plot this specific subplot
            % Delete old handles if they exist but are invalid for reuse
            if areScatterHandlesValid
                delete(currentScatterHandles);
            end
            
            % Get axis handle (we already validated axScatter exists above)
            axScatter = handles.axScatter(iCorr);
            
            % Re-plot the grouped scatter
            hOutScatter = plot_grouped_scatter(ampCurrent, ampNext, groupingData, ...
                'PlotEllipse', false, 'LinkXY', true, 'GridOn', true, ...
                'Color', plotParams.colorScatter, 'MarkerType', plotParams.markerTypeScatter, ...
                'MarkerSize', plotParams.markerSizeScatter, 'MarkerLineWidth', plotParams.markerLineWidthScatter, ...
                'FigTitle', 'suppress', 'LegendLocation', 'suppress', 'AxesHandle', axScatter);
            handles.hScatters{iCorr} = hOutScatter.dots;
        end
        
        % Get the current subplot handle
        axScatter = handles.axScatter(iCorr);

        % Update axis limits
        axisLimits = compute_axis_limits({ampCurrent, ampNext}, 'Coverage', axisCoveragePercScatter);
        axScatter.XLim = axisLimits;
        axScatter.YLim = axisLimits;
        
        % Update regression lines
        if numel(ampCurrent) > 2
            [hBestFits(iCorr), hTextBestFit(iCorr)] = ...
                plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
                                    'ThroughOrigin', false, 'Color', plotParams.colorBestFit, ...
                                    'LineStyle', plotParams.lineStyleBestFit, ...
                                    'LineWidth', plotParams.lineWidthBestFit, ...
                                    'AxesHandle', axScatter, ...
                                    'TextLocation', textLocBestFit, ...
                                    'ShowEquation', true, 'ShowRSquared', true, ...
                                    'ShowCorrCoeff', true);
            [hThrOrigs(iCorr), hTextThrOrig(iCorr)] = ...
                plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
                                    'ThroughOrigin', true, 'Color', plotParams.colorThrOrig, ...
                                    'LineStyle', plotParams.lineStyleThrOrig, ...
                                    'LineWidth', plotParams.lineWidthThrOrig, ...
                                    'AxesHandle', axScatter, ...
                                    'TextLocation', textLocThrOrig, ...
                                    'ShowEquation', true, 'ShowRSquared', true, ...
                                    'ShowCorrCoeff', false);
        end
    end

    % Store new handles
    handles.hBestFits = hBestFits;
    handles.hTextBestFit = hTextBestFit;
    handles.hThrOrigs = hThrOrigs;
    handles.hTextThrOrig = hTextThrOrig;

else
    % --- CREATE NEW PLOT ---
    [fig, ax] = create_subplots(nCorrToAnalyze, 'AlwaysNew', true, ...
                            'ShowFigure', showFigure, 'FigExpansion', [1, 1]);
    
    % Initialize storage for plot handles
    % hScatters changed to a cell array to support multiple groups
    hScatters = cell(nCorrToAnalyze, 1);
    hBestFits = gobjects(nCorrToAnalyze, 1);
    hTextBestFit = gobjects(nCorrToAnalyze, 1);
    hThrOrigs = gobjects(nCorrToAnalyze, 1);
    hTextThrOrig = gobjects(nCorrToAnalyze, 1);

    % Loop through each successive pair of whisk amplitudes
    for iCorr = 1:nCorrToAnalyze
        % Extract amplitude data for the current pair (e.g., A1 and A2)
        ampCurrent = whiskAmplitudesMatrix(:, iCorr);
        ampNext = whiskAmplitudesMatrix(:, iCorr + 1);

        % Plot the amplitudes against each other, grouped by the specified column
        hOutScatter = plot_grouped_scatter(ampCurrent, ampNext, groupingData, ...
            'PlotEllipse', false, 'LinkXY', true, 'GridOn', true, ...
            'Color', plotParams.colorScatter, 'MarkerType', plotParams.markerTypeScatter, ...
            'MarkerSize', plotParams.markerSizeScatter, 'MarkerLineWidth', plotParams.markerLineWidthScatter, ...
            'XLabel', sprintf('Amplitude of Whisk #%d (deg)', iCorr), ...
            'YLabel', sprintf('Amplitude of Whisk #%d (deg)', iCorr + 1), ...
            'FigTitle', 'suppress', 'LegendLocation', 'suppress', 'AxesHandle', ax(iCorr));
        % Store the returned handle(s) in a cell array
        hScatters{iCorr} = hOutScatter.dots;

        % Compute and plot the standard linear regression line
        [hBestFits(iCorr), hTextBestFit(iCorr), ~, ~] = ...
            plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
            'ThroughOrigin', false, 'Color', plotParams.colorBestFit, 'LineStyle', plotParams.lineStyleBestFit, ...
            'LineWidth', plotParams.lineWidthBestFit, 'AxesHandle', ax(iCorr), ...
            'TextLocation', textLocBestFit, 'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', true);
        
        % Compute and plot a regression line forced through the origin
        [hThrOrigs(iCorr), hTextThrOrig(iCorr)] = ...
            plot_regression_line('XData', ampCurrent, 'YData', ampNext, ...
            'ThroughOrigin', true, 'Color', plotParams.colorThrOrig, 'LineStyle', plotParams.lineStyleThrOrig, ...
            'LineWidth', plotParams.lineWidthThrOrig, 'AxesHandle', ax(iCorr), ...
            'TextLocation', textLocThrOrig, 'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', false);
    end
    
    % Add an overall title to the figure
    resize_subplots_for_labels('FigTitle', figTitle);

    % Store handles for output
    handles.fig = fig;
    handles.axScatter = ax;
    handles.hScatters = hScatters;
    handles.hBestFits = hBestFits;
    handles.hTextBestFit = hTextBestFit;
    handles.hThrOrigs = hThrOrigs;
    handles.hTextThrOrig = hTextThrOrig;
end

% Compute correlation coefficients for results output
corrCoeffs = nan(nCorrToAnalyze, 1);
pValues = nan(nCorrToAnalyze, 1);
for iCorr = 1:nCorrToAnalyze
    ampCurrent = whiskAmplitudesMatrix(:, iCorr);
    ampNext = whiskAmplitudesMatrix(:, iCorr + 1);
    toKeep = ~isnan(ampCurrent) & ~isnan(ampNext);
    if sum(toKeep) > 2
        [R, P] = corrcoef(ampCurrent(toKeep), ampNext(toKeep));
        corrCoeffs(iCorr) = R(1,2);
        pValues(iCorr) = P(1,2);
    end
end
results.corrCoeffs = corrCoeffs;
results.pValues = pValues;


% Save the figure to file if requested
if toSaveOutput
    figPath = fullfile(pathOutDir, figName); % Construct full path for saving the figure
    save_all_figtypes(fig, figPath, figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%