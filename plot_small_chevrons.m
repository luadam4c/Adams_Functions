function handles = plot_small_chevrons (chevronSheetPath, varargin)
%% Plots two chevron tables
% Usage: handles = plot_small_chevrons (chevronSheetPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, paths] = all_files('Extension', 'csv');
%       for i = 1:numel(paths); plot_small_chevrons(paths{i}); end
%       for i = 1:numel(paths); plot_small_chevrons(paths{i}, 'ColorMap', 'k', 'LegendLocation', 'suppress', 'PlotErrorBars', true); end
%       for i = 1:numel(paths); plot_small_chevrons(paths{i}, 'ColorMap', 'k', 'LegendLocation', 'suppress', 'PlotErrorBars', true, 'RunTTest', false, 'RunRankTest', false); end
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       chevronSheetPath    - TODO: Description of chevronSheetPath
%                           must be a TODO
%       varargin    - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - 'YLimitsLog2Ratio': limits of y axis 
%                                           for the log2 ratio plot
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - Any other parameter-value pair for plot_chevron()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_subplots.m
%       cd/extract_fileparts.m
%       cd/plot_chevron.m
%       cd/save_all_figtypes.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-10-18 Adapted from plot_relative_events.m
% 

%% Hard-coded parameters

% TODO: Make optional arguments
% readoutLabel = 'SWD Count'
% readoutLabel = 'Oscillation Duration';
readoutLabel = 'Oscillation Period';
plotLog2Ratio = true; %false;
figTypes = {'png', 'epsc2'};
% pTickLabels = {'Before', 'After'};
pTickLabels = {'Baseline', 'GAT Block'};
varsToExclude = {'phaseNumber'};
log2Fun = @(x) log2(x);
figHeight = 2;
figWidth = 3.0;
alwaysNew = true;

%% Default values for optional arguments
yLimitsDefault = [];            % set later
yLimitsLog2RatioDefault = [];  % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'chevronSheetPath');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimitsLog2Ratio', yLimitsLog2RatioDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, chevronSheetPath, varargin{:});
yLimits = iP.Results.YLimits;
yLimitsLog2Ratio = iP.Results.YLimitsLog2Ratio;

% Keep unmatched arguments for the plot_chevron() function
otherArguments = iP.Unmatched;

%% Preparation
% Extract the file base
chevronPathBase = extract_fileparts(chevronSheetPath, 'pathbase');

% Create figure path base
figPathBase = chevronPathBase;

% Create the readout label
readoutLabelLog2 = [readoutLabel, ' ratio'];

%% Do the job
% Read the table
chevronTable = readtable(chevronSheetPath);

% Remove variables to exclude
chevronTable = removevars(chevronTable, varsToExclude);

% Compute log2 data
if height(chevronTable) < width(chevronTable)
    log2ChevronTable = varfun(log2Fun, chevronTable);
    log2ratioChevronTable = varfun(@(x) x/x(1), log2ChevronTable);
else
    log2ChevronTable = rowfun(log2Fun, chevronTable);
    log2ratioChevronTable = rowfun(@(x) x/x(1), log2ChevronTable);
end

% Create figure
fig = set_figure_properties('Units', 'inches', 'Height', figHeight, ...
                    'Width', figWidth, 'AlwaysNew', alwaysNew);

% Create subplots
if plotLog2Ratio
    [fig, ax] = create_subplots(1, 2, 'FigExpansion', [1, 1], ...
                                'FigHandle', fig);
else
    [fig, ax] = create_subplots(1, 1, 'FigHandle', fig);
end

% Plot Chevron plot and save figure
plot_chevron(chevronTable, 'FigTitle', 'suppress', ...
            'ReadoutLabel', readoutLabel, 'PTickLabels', pTickLabels, ...
            'ReadoutLimits', yLimits, 'LegendLocation', 'northeast', ...
            'AxesHandle', ax(1), 'FigExpansion', [], otherArguments);

% Plot log2 ratio Chevron plot and save figure
if plotLog2Ratio
    plot_chevron(log2ratioChevronTable, 'FigTitle', 'suppress', ...
                'IsLog2Data', true, ...
                'ReadoutLabel', readoutLabelLog2, ...
                'PTickLabels', pTickLabels, ...
                'ReadoutLimits', yLimitsLog2Ratio, ...
                'LegendLocation', 'northeast', ...
                'AxesHandle', ax(2), 'FigExpansion', [], otherArguments);
end

% Save figure
save_all_figtypes(fig, figPathBase, figTypes);

%% Output results
handles.fig = fig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%