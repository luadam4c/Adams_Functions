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
%                   - Any other parameter-value pair for violinplot.m
%
% Requires:
%       cd/addpath_custom.m
%       cd/apply_iteratively.m
%       cd/combine_strings.m
%       cd/create_error_for_nargin.m
%       cd/decide_on_colormap.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/locate_functionsdir.m
%       cd/isfigtype.m
%       cd/ispositiveintegervector.m
%       cd/struct2arglist.m
%       cd/save_all_figtypes.m
%       cd/set_figure_properties.m
%       cd/update_figure_for_corel.m
%       ~/Downloaded_Functions/Violinplot/violinplot.m
%
% Used by:
%       cd/m3ha_plot_figure02.m
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
figWidthDefault = 5;
figHeightDefault = 3.4;
figTypesDefault = {'png', 'epsc2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for violinplot()
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions', ...
                    'Violinplot'));
end

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

% Keep unmatched arguments for the violinplot() function
otherArguments = struct2arglist(iP.Unmatched);

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

function handles = m3ha_plot_violin_helper (allValues, violinRelativeBandWidth, ...
                            medianColor, medianSize, ...
                            measureTitle, pharmLabels, ...
                            figPathBase, figHeight, figWidth, ...
                            figTypes, otherArguments)

% Hard-coded parameters
MS_PER_S = 1000;
xTickAngle = 320;

% Create figure for conductance traces
fig = set_figure_properties('AlwaysNew', true);

% Count the number of groups
nGroups = numel(pharmLabels);

% Convert onset times from ms to seconds
%{
if contains(measureTitle, 'onset')
    % Update values
    allValues = cellfun(@(x) x ./ MS_PER_S, allValues, 'UniformOutput', false);

    % Update title
    measureTitle = replace(measureTitle, 'ms', 's');
end
%}

% Decide on the color map
cm = decide_on_colormap([], nGroups);

% Set the color map
%   TODO: Apply this in violinplot?
colormap(cm);

% TODO: plot_violin.m
% TODO: plot_jitter.m

% Force as a numeric array
allValues = force_matrix(allValues);

% If there is any column with all NaN values, remove it from the violin plot
if any(all(isnan(allValues), 1))
    % Test whether each column has no values
    noValues = all(isnan(allValues), 1);

    % Remove those columns
    disp('Groups with all NaN values are removed!');
    allValues = allValues(:, ~noValues);
    pharmLabels = pharmLabels(~noValues);
end

% Compute range of all values
rangeValues = apply_iteratively(@max, allValues) - apply_iteratively(@min, allValues);

% Compute the bandwidth for the kernel density estimates
bandWidth = violinRelativeBandWidth * rangeValues;

% Plot a violin plot
% violinplot(allValues, pharmLabels);
violins = violinplot(allValues, pharmLabels, 'BandWidth', bandWidth, ...
                        otherArguments{:});

% Plot the data points for each cell
% plotSpread(allValues);
% Set x tick labels
% xticklabels(pharmLabels);

% Modify x limits
xlim([0.5, nGroups + 0.5]);

% Modify x tick angle
xtickangle(xTickAngle);

% Set y label
ylabel(measureTitle);

% Find all median scatters and make the face color green and size bigger
medianScatters = findobj(gca, 'Type', 'Scatter', ...
                        'MarkerEdgeColor', [0.5, 0.5, 0.5]);
set(medianScatters, 'MarkerFaceColor', medianColor);

% Save the figure
save_all_figtypes(fig, [figPathBase, '_orig'], 'png');

% Update figure for CorelDraw
update_figure_for_corel(fig, 'Units', 'centimeters', ...
                        'Height', figHeight, 'Width', figWidth, ...
                        'ScatterMarkerSize', 3);

% Fix axes position
set(gca, 'Position', [0.2356, 0.1947, 0.6694, 0.7303]);

% Update median size
set(medianScatters, 'SizeData', medianSize^2);

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