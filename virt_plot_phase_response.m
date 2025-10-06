function [results, handles] = virt_plot_phase_response (T, pPlot, varargin)
%% Plots a phase response curve from aggregated data
% Usage: [results, handles] = virt_plot_phase_response (T, pPlot, varargin)
% Explanation:
%       Plots phaseReset vs. phaseChangeWhisk, grouped by a specified
%       column. Also plots a regression line with statistics.
%
% Outputs:
%       results     - A structure containing regression analysis results.
%       handles     - A structure containing handles to the plot objects.
%
% Arguments:
%       T           - A table containing whisk analysis data.
%                   must be a table
%       pPlot       - The plotting parameters structure from virt_moore.m.
%                   must be a structure
%       varargin    - 'GroupingColumn': The name of the column to group data by.
%                   must be a string scalar or a character vector
%                   default == 'fileNumber'
%                   - 'WhiskDir': Direction of whisk used for phase calculations.
%                   must be a string scalar or character vector
%                   default == 'retraction'
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
%
% Requires:
%       /Shared/Code/Adams_Functions/plot_grouped_scatter.m
%       /Shared/Code/Adams_Functions/plot_regression_line.m
%       /Shared/Code/Adams_Functions/save_all_figtypes.m
%       /Shared/Code/Adams_Functions/set_figure_properties.m
%
% Used by:
%       TODO: cd/virt_analyze_sniff_whisk.m
%       TODO: \Shared\Code\vIRt-Moore\virt_plot_whisk_analysis.m
%       \Shared\Code\vIRt-Moore\virt_moore_monte_carlo.m

% File History:
% 2025-10-02 Created by Gemini
%

%% Default values for optional arguments
groupingColumnDefault = 'fileNumber';
whiskDirDefault = 'retraction';
figTitleDefault = 'Whisk Phase Response Curve';
figNameDefault = 'phase_response_scatter';
outDirDefault = pwd;
figTypesDefault = {'png'};
toSaveOutputDefault = true;                 % Default to save the output figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addRequired(iP, 'T', @istable);
addRequired(iP, 'pPlot', @isstruct);
addParameter(iP, 'GroupingColumn', groupingColumnDefault, @ischar);
addParameter(iP, 'WhiskDir', whiskDirDefault, @ischar);
addParameter(iP, 'FigTitle', figTitleDefault, @ischar);
addParameter(iP, 'FigName', figNameDefault, @ischar);
addParameter(iP, 'OutDir', outDirDefault, @ischar);
addParameter(iP, 'FigTypes', figTypesDefault);
addParameter(iP, 'ToSaveOutput', toSaveOutputDefault, @islogical);

% Read from the Input Parser
parse(iP, T, pPlot, varargin{:});
groupingColumn = iP.Results.GroupingColumn;
whiskDirForPhase = iP.Results.WhiskDir;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
pathOutDir = iP.Results.OutDir;
figTypes = iP.Results.FigTypes;
toSaveOutput = iP.Results.ToSaveOutput;

%% Preparation
results = struct;
handles.fig = gobjects;

if isempty(T) || ~ismember('phaseReset', T.Properties.VariableNames)
    fprintf('No phase response data to plot. Skipping PRC plot!\n');
    return;
end

% Extract data and filter NaNs
phaseReset = T.phaseReset;
phaseChangeWhisk = T.phaseChangeWhisk;
groupingData = T.(groupingColumn);
toKeep = ~isnan(phaseReset) & ~isnan(phaseChangeWhisk);
phaseReset = phaseReset(toKeep);
phaseChangeWhisk = phaseChangeWhisk(toKeep);
groupingData = groupingData(toKeep);

if isempty(phaseReset)
    fprintf('No valid phase response data found to plot.\n');
    return;
end

%% Plot scatter plot and regression
fig = set_figure_properties('AlwaysNew', true, 'ClearFigure', true);
figPath = fullfile(pathOutDir, figName);

% Set up axes properties
ax = set_axes_properties;

plot_grouped_scatter(phaseReset, phaseChangeWhisk, groupingData, ...
    'PlotEllipse', false, 'LinkXY', false, 'GridOn', true, ...
    'XLabel', ['Phase of breath onset in whisk ', whiskDirForPhase, ' cycle (radians)'], ...
    'YLabel', ['Phase change of following whisk ', whiskDirForPhase, ' (radians)'], ...
    'XLimits', [0, 2*pi], 'YLimits', [-2*pi, 2*pi], ...
    'FigTitle', figTitle, 'LegendLocation', 'suppress');

xticks(pi * (0:0.5:2));
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
yticks(pi * (-2:1:2));
yticklabels({'-2\pi', '-\pi', '0', '\pi', '2\pi'});

line([0, 2*pi], [-2*pi, 0], 'Color', 'green', 'LineStyle', '-', 'LineWidth', 1);
yline(0, '--g', 'LineWidth', 1);

[~, ~, ~, regResults] = plot_regression_line('AxesHandle', ax, ...
    'ShowEquation', true, 'ShowRSquared', true, 'ShowCorrCoeff', true);

% Save the figure to file if requested
if toSaveOutput
    save_all_figtypes(fig, figPath, figTypes);
end

% Store handles and results for output
handles.fig = fig;
handles.ax = ax;
results = regResults;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%