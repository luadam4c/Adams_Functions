function [results, handles] = virt_plot_phase_response (dataTable, pPlot, varargin)
%% Plots a phase response curve from aggregated data
% Usage: [results, handles] = virt_plot_phase_response (dataTable, pPlot, varargin)
% Explanation:
%       Plots phaseReset vs. phaseChangeWhisk, grouped by a specified
%       column. Also plots a regression line with statistics. It can create
%       a new figure or update an existing one.
%
% Outputs:
%       results     - A structure containing regression analysis results.
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
%
% Requires:
%       cd/plot_grouped_scatter.m
%       cd/plot_regression_line.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
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
%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addRequired(iP, 'dataTable', @istable);
addRequired(iP, 'pPlot', @isstruct);
addParameter(iP, 'GroupingColumn', groupingColumnDefault);
addParameter(iP, 'WhiskDir', whiskDirDefault, @ischar);
addParameter(iP, 'Handles', handlesDefault);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);
addParameter(iP, 'ShowFigure', showFigureDefault, @islogical);

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

%% Plot scatter plot and regression
if ~isempty(handlesIn) && isfield(handlesIn, 'fig') && isgraphics(handlesIn.fig)
    % --- UPDATE EXISTING PLOT ---
    handles = handlesIn;
    fig = handles.fig;
    axPRC = handles.axPRC;
    
    % Update scatter data for each group
    uniqueGroups = unique(groupingData);
    nGroups = numel(uniqueGroups);
    currentScatterHandles = handles.hPRCScatter;

    if nGroups == numel(currentScatterHandles)
        for iGroup = 1:nGroups
            groupValue = uniqueGroups(iGroup);
            groupMask = (groupingData == groupValue);
            set(currentScatterHandles(iGroup), ...
                'XData', phaseReset(groupMask), ...
                'YData', phaseChangeWhisk(groupMask));
        end
    else
        % Fallback if number of groups changes
        delete(currentScatterHandles);
        hOutScatter = plot_grouped_scatter(phaseReset, phaseChangeWhisk, groupingData, ...
            'AxesHandle', axPRC, 'PlotEllipse', false, 'LinkXY', false, 'GridOn', true, ...
            'XLabel', get(get(axPRC, 'XLabel'), 'String'), 'YLabel', get(get(axPRC, 'YLabel'), 'String'), ...
            'LegendLocation', 'suppress');
        handles.hPRCScatter = hOutScatter.dots;
    end

    % Delete and replot regression line
    delete(handles.hPRCRegLine);
    delete(handles.hPRCRegText);
    if sum(toKeep) > 2
        [handles.hPRCRegLine, handles.hPRCRegText, ~, regResults] = ...
            plot_regression_line('XData', phaseReset, 'YData', phaseChangeWhisk, ...
                             'AxesHandle', axPRC, 'ShowEquation', true, ...
                             'ShowRSquared', true, 'ShowCorrCoeff', true);
    else
        handles.hPRCRegLine = gobjects;
        handles.hPRCRegText = gobjects;
        regResults = struct;
    end
else
    % --- CREATE NEW PLOT ---
    fig = set_figure_properties('AlwaysNew', true, 'ShowFigure', showFigure);
    ax = gca;

    % Generate the scatter plot
    hOutScatter = plot_grouped_scatter(phaseReset, phaseChangeWhisk, groupingData, ...
        'AxesHandle', ax, 'PlotEllipse', false, 'LinkXY', false, 'GridOn', true, ...
        'XLabel', ['Phase of breath onset in whisk ', whiskDirForPhase, ' cycle (radians)'], ...
        'YLabel', ['Phase change of following whisk ', whiskDirForPhase, ' (radians)'], ...
        'XLimits', [0, 2*pi], 'YLimits', [-2*pi, 2*pi], ...
        'FigTitle', figTitle, 'LegendLocation', 'suppress');
    hPRCScatter = hOutScatter.dots;
    hold(ax, 'on');

    % Set up tick marks and labels
    xticks(pi * (0:0.5:2));
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    yticks(pi * (-2:1:2));
    yticklabels({'-2\pi', '-\pi', '0', '\pi', '2\pi'});

    % Plot reference lines
    line([0, 2*pi], [-2*pi, 0], 'Color', 'green', 'LineStyle', '-', 'LineWidth', 1);
    yline(0, '--g', 'LineWidth', 1);

    % Plot regression and get results
    if numel(phaseReset) > 2
        [hPRCRegLine, hPRCRegText, ~, regResults] = plot_regression_line('AxesHandle', ax, ...
            'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', true);
    else
        hPRCRegLine = gobjects;
        hPRCRegText = gobjects;
        regResults = struct;
    end
    hold(ax, 'off');

    % Store handles for output
    handles.fig = fig;
    handles.axPRC = ax;
    handles.hPRCScatter = hPRCScatter;
    handles.hPRCRegLine = hPRCRegLine;
    handles.hPRCRegText = hPRCRegText;
end


% Save the figure to file if requested
if toSaveOutput
    figPath = fullfile(pathOutDir, figName);
    save_all_figtypes(fig, figPath, figTypes);
end

% Re-compute regression for results output if not done
if ~exist('regResults', 'var')
    if numel(phaseReset) > 2
        [~, ~, ~, regResults] = plot_regression_line('XData', phaseReset, 'YData', phaseChangeWhisk, 'ToPlot', false);
    else
        regResults = struct;
    end
end
results = regResults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%