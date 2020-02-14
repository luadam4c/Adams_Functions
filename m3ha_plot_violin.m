function handles = m3ha_plot_violin (statsPath, varargin)
%% Plots violin plots from a statistics table returned by m3ha_compute_statistics.m
% Usage: handles = m3ha_plot_violin (statsPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - structures with fields:
%                       fig
%                       violins
%                   specified as a structure array
%
% Arguments:
%       statsPath  - path to a .mat file containing the variables:
%                       statsTable - statistics table 
%                                       returned by m3ha_compute_statistics.m
%                           must have these columns:
%                               measureTitle
%                               measureStr
%                               allValues
%                       pharmLabels
%                       conditionLabel
%                   must be a string scalar or a character vector
%       varargin    - 'RowsToPlot': rows to plot
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'OutFolder': the directory where plots will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigWidth': figure width
%                   must be a positive scalar
%                   default == 5 cm
%                   - 'FigHeight': figure height
%                   must be a positive scalar
%                   default == 3.4 cm
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for plot_violin.m
%
% Requires:
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/isfigtype.m
%       cd/ispositiveintegervector.m
%       cd/plot_violin.m
%       cd/struct2arglist.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m
%       ~/Downloaded_Functions/Violinplot/violinplot.m
%
% Used by:
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_plot_figure07.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-12-27 Moved from m3ha_plot_figure02.m
% 2019-12-28 Added 'FigTypes' as an optional argument
% 2019-12-28 Added 'FigHeight' and 'FigWidth' as optional arguments
% 

%% Hard-coded parameters
violinRelativeBandWidth = 0.1;  % bandwidth relative to data range
medianColor = [0.6758, 1, 0.1836];     % rgb('GreenYellow') 
                                % color of median circle
medianSize = 6;                 % size of median circle in points

%% Default values for optional arguments
rowsToPlotDefault = 'all';
outFolderDefault = '';          % set later
figWidthDefault = 3.4; %5;
figHeightDefault = 3; %3.4;
figTypesDefault = {'png', 'epsc2'};

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
addRequired(iP, 'statsPath', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RowsToPlot', rowsToPlotDefault, ...
    @(x) assert(ispositiveintegervector(x) || ischar(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['RowsToPlot must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigWidth', figWidthDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigWidth must be a empty or a positive scalar!'));
addParameter(iP, 'FigHeight', figHeightDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                'FigHeight must be a empty or a positive scalar!'));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, statsPath, varargin{:});
rowsToPlot = iP.Results.RowsToPlot;
outFolder = iP.Results.OutFolder;
figWidth = iP.Results.FigWidth;
figHeight = iP.Results.FigHeight;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_violin() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output
handles = struct;

% Set default output directory
if isempty(outFolder)
    outFolder = extract_fileparts(statsPath, 'directory');
end

% Load stats table
disp('Loading statistics for violin plots ...');
if isfile(statsPath)
    load(statsPath, 'statsTable', 'pharmLabels', 'conditionLabel');
else
    fprintf('%s does not exist!\n', statsPath);
    return;
end

% Restrict to measures to plot
if ~(ischar(rowsToPlot) && strcmp(rowsToPlot, 'all'))
    statsTable = statsTable(rowsToPlot, :);
end

% Extract variables
allMeasureTitles = statsTable.measureTitle;
allMeasureStrs = statsTable.measureStr;
allValues = statsTable.allValues;

% Create figure bases
allFigBases = combine_strings({allMeasureStrs, conditionLabel});

% Create full path bases
allFigPathBases = fullfile(outFolder, allFigBases);

%% Do the job
% Plot all 2D violin plots
disp('Plotting 2D violin plots ...');
handles = cellfun(@(a, b, c) m3ha_plot_violin_helper(...
                            a, violinRelativeBandWidth, ...
                            medianColor, medianSize, b, ...
                            pharmLabels, c, ...
                            figHeight, figWidth, ...
                            figTypes, otherArguments), ...
                    allValues, allMeasureTitles, allFigPathBases);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_violin_helper (allValues, relativeBandWidth, ...
                                            medianColor, medianSize, ...
                                            measureTitle, pharmLabels, ...
                                            figPathBase, figHeight, figWidth, ...
                                            figTypes, otherArguments)

% Hard-coded parameters
MS_PER_S = 1000;
xTickAngle = 320;

% Compute statistics
% TODO

% Create figure for conductance traces
fig = set_figure_properties('AlwaysNew', true);

% Convert onset times from ms to seconds
%{
if contains(measureTitle, 'onset')
    % Update values
    allValues = cellfun(@(x) x ./ MS_PER_S, allValues, 'UniformOutput', false);

    % Update title
    measureTitle = replace(measureTitle, 'ms', 's');
end
%}

% Plot groups as a violin plot
violins = plot_violin(allValues, 'XTickLabels', pharmLabels, ...
                        'RelativeBandWidth', relativeBandWidth, ...
                        'MedianColor', medianColor, ...
                        'XTickAngle', xTickAngle, 'YLabel', measureTitle, ...
                        otherArguments);
% TODO: plot_jitter.m
% Plot the data points for each cell
% plotSpread(allValues);
% Set x tick labels
% xticklabels(pharmLabels);

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Set y axis limits based on measureTitle
switch measureTitle
    case 'LTS probability'
        ylim([0, 1]);
    case 'LTS onset time (ms)'
        ylim([0, 2000]);
    case 'Spikes Per LTS'
        ylim([0, 6.5]);
    case 'LTS maximum slope (V/s)'
        ylim([0, 5]);
    case 'LTS amplitude (mV)'
        ylim([-75, -45]);
        yticks(-75:10:-45);
    otherwise
        % Do nothing
end

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'ScatterMarkerSize', 3);

% Fix axes position
set(gca, 'Position', [0.2356, 0.1947, 0.6694, 0.7303]);

% Update median size
medianPlots = extract_fields(violins, 'MedianPlot', 'UniformOutput', true);
medianPlots = medianPlots(ishandle(medianPlots));
set(medianPlots, 'SizeData', medianSize^2);

% Save the figure
save_all_figtypes(fig, figPathBase, figTypes);

% Save in handles
handles.fig = fig;
handles.violins = violins;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%