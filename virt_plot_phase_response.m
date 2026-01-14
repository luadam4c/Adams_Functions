function [results, handles] = virt_plot_phase_response (dataTable, pPlot, varargin)
%% Plots a phase response curve from aggregated data
% Usage: [results, handles] = virt_plot_phase_response (dataTable, pPlot, varargin)
% Explanation:
%       Plots phaseReset vs. phaseChangeWhisk, grouped by a specified
%       column. 
%       - Plots a regression line with statistics.
%       - Optionally calculates and plots the mean/median change in phase 
%         in non-overlapping bins (default pi/4) and tests if significantly 
%         different from zero using test_difference.m.
%       - Optionally plots the distribution of phase resets and tests for 
%         uniformity.
%
% Outputs:
%       results     - A structure containing regression analysis results,
%                     binned statistics, and uniformity test results.
%       handles     - A structure containing handles to the plot objects.
%
% Arguments:
%       dataTable           - A table containing whisk analysis data.
%                   must be a table
%       pPlot       - The plotting parameters structure from virt_moore.m.
%                   must be a structure
%       varargin    - 'GroupingColumn': The name of the column to group data by.
%                   must be a string scalar or a character vector
%                   default == 'seedNumber' or 'fileNumber'
%                   - 'WhiskDir': Direction of whisk used for phase calculations.
%                   must be a string scalar or character vector
%                   default == 'retraction'
%                   - 'Handles': A structure of graphics handles to update.
%                   must be a structure
%                   default == []
%                   - 'FigTitle': The title for the figure.
%                   must be a string scalar or a character vector
%                   default == 'Whisk Phase Response Curve'
%                   - 'FigName': The base file name for saving the figure.
%                   must be a string scalar or a character vector
%                   default == 'phase_response_scatter'
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
%                   - 'BinWidth': Width of bins for binned analysis (radians)
%                   must be a positive scalar
%                   default == pi/4
%                   - 'ShowBinStats': Whether to show binned statistics and histogram
%                   must be logical
%                   default == true
%
% Requires:
%       cd/create_subplots.m
%       cd/plot_grouped_scatter.m
%       cd/plot_regression_line.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/test_difference.m
%       cd/test_normality.m
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
% 2025-10-17 Added 'ShowFigure' as an optional argument.
% 2025-11-26 Added binned phase analysis and uniformity distribution plot.

%% Default values for optional arguments
groupingColumnDefault = [];                 % set later
whiskDirDefault = 'retraction';
handlesDefault = [];
figTitleDefault = 'Whisk Phase Response Curve';
figNameDefault = 'phase_response_scatter';
outDirDefault = pwd;
figTypesDefault = {'png'};
toSaveOutputDefault = true;                 % Default to save the output figure
showFigureDefault = true;                   % Default to show the figure
binWidthDefault = pi/4;                     % Default bin width for analysis
showBinStatsDefault = true;                 % Default to show binned stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true; 

% Add required inputs
addRequired(iP, 'dataTable', @istable);
addRequired(iP, 'pPlot', @isstruct);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'GroupingColumn', groupingColumnDefault);
addParameter(iP, 'WhiskDir', whiskDirDefault, @ischar);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);
addParameter(iP, 'ShowFigure', showFigureDefault, @islogical);
addParameter(iP, 'BinWidth', binWidthDefault, @isnumeric);
addParameter(iP, 'ShowBinStats', showBinStatsDefault, @islogical);

% Read from the Input Parser
parse(iP, dataTable, pPlot, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
whiskDirForPhase = iP.Results.WhiskDir;
handlesIn = iP.Results.Handles;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;
showFigure = iP.Results.ShowFigure;
binWidth = iP.Results.BinWidth;
showBinStats = iP.Results.ShowBinStats;

%% Preparation
% Initialize outputs
results = struct;
handles.fig = gobjects;

% Check if phase response data is present
if isempty(dataTable) || ~ismember('phaseReset', dataTable.Properties.VariableNames)
    fprintf('No phase response data to plot. Skipping PRC plot!\n');
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

% Extract data and filter NaNs
phaseReset = dataTable.phaseReset;
phaseChangeWhisk = dataTable.phaseChangeWhisk;
groupingData = dataTable.(groupingColumn);
toKeep = ~isnan(phaseReset) & ~isnan(phaseChangeWhisk);
phaseReset = phaseReset(toKeep);
phaseChangeWhisk = phaseChangeWhisk(toKeep);
groupingData = groupingData(toKeep);

if isempty(phaseReset)
    fprintf('No valid phase response data found to plot.\n');
    return;
end

%% Setup Figure and Axes
if ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig)
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn;
    fig = handles.fig;
    axPRC = handles.axPRC;
    % Check if histogram axis exists, otherwise we might skip it
    if isfield(handles, 'axHist')
        axHist = handles.axHist;
    else
        axHist = gobjects; 
    end
else
    % --- CREATE NEW PLOT ---
    if showBinStats
        % Create 2 rows: Top for PRC, Bottom for Distribution
        % Pass [] for gridPositions to use default
        [fig, axs] = create_subplots(2, 1, [], ...
            'ShowFigure', showFigure);
        axPRC = axs(1);
        axHist = axs(2);
    else
        [fig, axs] = create_subplots(1, 1, [], ...
            'ShowFigure', showFigure);
        axPRC = axs(1);
        axHist = gobjects;
    end
    
    handles.fig = fig;
    handles.axPRC = axPRC;
    handles.axHist = axHist;
end

%% 1. Plot Main Scatter and Regression (Top Panel)
if ~isempty(handlesIn) && isfield(handlesIn, 'fig')
     % Update logic for scatter
    uniqueGroups = unique(groupingData);
    nGroups = numel(uniqueGroups);
    if isfield(handles, 'hPRCScatter') && nGroups == numel(handles.hPRCScatter)
        for iGroup = 1:nGroups
            groupValue = uniqueGroups(iGroup);
            groupMask = (groupingData == groupValue);
            set(handles.hPRCScatter(iGroup), ...
                'XData', phaseReset(groupMask), ...
                'YData', phaseChangeWhisk(groupMask));
        end
    else
         % Fallback: Replot if grouping changed
         cla(axPRC);
         hOutScatter = plot_grouped_scatter(phaseReset, phaseChangeWhisk, groupingData, ...
            'AxesHandle', axPRC, 'PlotEllipse', false, 'LinkXY', false, 'GridOn', true, ...
            'XLabel', ['Phase of breath onset in whisk ', whiskDirForPhase, ' cycle (radians)'], ...
            'YLabel', ['Phase change of following whisk ', whiskDirForPhase, ' (radians)'], ...
            'XLimits', [0, 2*pi], 'YLimits', [-2*pi, 2*pi], ...
            'LegendLocation', 'suppress');
         handles.hPRCScatter = hOutScatter.dots;
    end
    
    % Re-run regression plotting
    if isfield(handles, 'hPRCRegLine') && isvalid(handles.hPRCRegLine)
        delete(handles.hPRCRegLine); 
    end
    if isfield(handles, 'hPRCRegText') && isvalid(handles.hPRCRegText)
        delete(handles.hPRCRegText);
    end
    
    [handles.hPRCRegLine, handles.hPRCRegText, ~, regResults] = ...
            plot_regression_line('XData', phaseReset, 'YData', phaseChangeWhisk, ...
                             'AxesHandle', axPRC, 'ShowEquation', true, ...
                             'ShowRSquared', true, 'ShowCorrCoeff', true);
else
    % Create Scatter
    hOutScatter = plot_grouped_scatter(phaseReset, phaseChangeWhisk, groupingData, ...
        'AxesHandle', axPRC, 'PlotEllipse', false, 'LinkXY', false, 'GridOn', true, ...
        'XLabel', ['Phase of breath onset in whisk ', whiskDirForPhase, ' cycle (radians)'], ...
        'YLabel', ['Phase change of following whisk ', whiskDirForPhase, ' (radians)'], ...
        'XLimits', [0, 2*pi], 'YLimits', [-2*pi, 2*pi], ...
        'LegendLocation', 'suppress');
    handles.hPRCScatter = hOutScatter.dots;
    hold(axPRC, 'on');

    % Set up tick marks and labels
    xticks(axPRC, pi * (0:0.5:2));
    xticklabels(axPRC, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticks(axPRC, pi * (-2:1:2));
    yticklabels(axPRC, {'-2\pi', '-\pi', '0', '\pi', '2\pi'});

    % Plot reference lines
    line(axPRC, [0, 2*pi], [-2*pi, 0], 'Color', 'green', 'LineStyle', '-', 'LineWidth', 1);
    yline(axPRC, 0, '--g', 'LineWidth', 1);

    % Plot regression
    [handles.hPRCRegLine, handles.hPRCRegText, ~, regResults] = plot_regression_line('AxesHandle', axPRC, ...
            'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', true);
end

%% 2. Calculate and Plot Binned Stats (Top Panel)
if showBinStats
    % Define bins
    edges = 0:binWidth:2*pi;
    binCenters = edges(1:end-1) + binWidth/2;
    nBins = numel(binCenters);
    
    % Initialize storage
    binMeans = nan(nBins, 1);
    binErrs = nan(nBins, 1);
    binPvals = nan(nBins, 1);
    binNormality = false(nBins, 1);
    
    for i = 1:nBins
        % Find data in this bin
        inBin = phaseReset >= edges(i) & phaseReset < edges(i+1);
        ySub = phaseChangeWhisk(inBin);
        
        if sum(~isnan(ySub)) >= 2
            % Test difference from zero using provided function
            % Note: test_difference defaults to testing against 0 for one group
            % and handles normality checking internally.
            stats = test_difference(ySub, 'DisplayAnova', false);
            
            % Check if normal to decide what metric to plot
            if stats.isNormal_Group1
                % Normal: Use Mean and Standard Error
                binMeans(i) = stats.meanDiff; 
                binErrs(i) = stats.stderrDiff;
                binNormality(i) = true;
            else
                % Non-Normal: Use Median
                % Note: test_difference output 'meanDiff' usually contains the
                % mean regardless of normality, but we want the median if the
                % test was non-parametric (signrank).
                if strcmpi(stats.testFunction, 'signrank')
                    binMeans(i) = median(ySub, 'omitnan');
                else
                    binMeans(i) = stats.meanDiff;
                end
                
                % For visualization of error in non-parametric cases, we
                % approximate or use the standard error provided for consistency
                % of visual weight, but results struct contains exact p-value.
                binErrs(i) = stats.stderrDiff; 
            end
            binPvals(i) = stats.pValue;
        end
    end
    
    % Remove old error bars if updating
    if isfield(handles, 'hBinError') && isvalid(handles.hBinError)
        delete(handles.hBinError);
    end
    if isfield(handles, 'hBinSig') && isvalid(handles.hBinSig)
        delete(handles.hBinSig);
    end

    % Plot Error Bars
    hold(axPRC, 'on');
    handles.hBinError = errorbar(axPRC, binCenters, binMeans, binErrs, ...
        'Color', 'r', 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'CapSize', 8);
    
    % Annotate significant bins
    sigBins = binPvals < 0.05;
    if any(sigBins)
        % Place stars slightly above the error bar
        handles.hBinSig = plot(axPRC, binCenters(sigBins), ...
            binMeans(sigBins) + binErrs(sigBins) + 0.2, ...
            'r*', 'MarkerSize', 8);
    else
        handles.hBinSig = gobjects;
    end
    hold(axPRC, 'off');

    % Add Title
    title(axPRC, figTitle);
    
    % Store in results
    results.binned.centers = binCenters;
    results.binned.means = binMeans;
    results.binned.errs = binErrs;
    results.binned.pValues = binPvals;
    results.binned.isNormal = binNormality;
end

%% 3. Plot Distribution and Test Uniformity (Bottom Panel)
if showBinStats && isgraphics(axHist)
    % Calculate Histogram
    [counts, ~] = histcounts(phaseReset, edges);
    
    % Plot Bar Chart
    handles.hHist = bar(axHist, binCenters, counts, 1, 'FaceColor', [0.7 0.7 0.7]);
    
    % Formatting
    xlabel(axHist, 'Phase of breath onset (radians)');
    ylabel(axHist, 'Count');
    xticks(axHist, pi * (0:0.5:2));
    xticklabels(axHist, {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    xlim(axHist, [0, 2*pi]);
    grid(axHist, 'on');
    
    % Test for Uniformity (Chi-Square Goodness of Fit)
    % Observed = counts
    % Expected = Total / nBins (for uniform assumption)
    nTotal = sum(counts);
    
    if nTotal > 0
        expected = repmat(nTotal / nBins, 1, nBins);

        % Chi-Square Statistic: sum((O-E)^2 / E)
        % Handle cases where expected is 0 to avoid Inf
        validBins = expected > 0;
        chi2stat = sum((counts(validBins) - expected(validBins)).^2 ./ expected(validBins));
        df = sum(validBins) - 1;

        % P-value (Upper tail of Chi-square distribution)
        if df > 0
            pUniform = chi2cdf(chi2stat, df, 'upper');
        else
            pUniform = NaN;
        end
    else
        chi2stat = NaN;
        pUniform = NaN;
    end
    
    % Add Title with Stats
    if ~isnan(pUniform)
        if pUniform < 0.05
            unifStatus = 'Non-Uniform';
        else
            unifStatus = 'Uniform';
        end
        title(axHist, sprintf('Distribution of Phase Resets (%s, p=%.3f)', ...
            unifStatus, pUniform));
    else
        title(axHist, 'Distribution of Phase Resets');
    end
    
    % Store in results
    results.uniformity.counts = counts;
    results.uniformity.pValue = pUniform;
    results.uniformity.chi2stat = chi2stat;
end

%% Finalize Output
results.regression = regResults;

% Save the figure to file if requested
if toSaveOutput
    figPath = fullfile(pathOutDir, figName);
    save_all_figtypes(fig, figPath, figTypes);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%