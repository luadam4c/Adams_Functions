function handles = plot_tuning_curve (pValues, readout, varargin)
%% Plot 1-dimensional tuning curve(s), can include confidence intervals or test p values
% Usage: handles = plot_tuning_curve (pValues, readout, varargin)
% Explanation:
%       TODO
% 
% Examples:
%       pValues = transpose(1:10);
%       readout1 = randi(numel(pValues), size(pValues));
%       upperCI1 = readout1 + randi(numel(pValues), 10, 1) / 10;
%       lowerCI1 = readout1 - randi(numel(pValues), 10, 1) / 10;
%       readout2 = randi(numel(pValues), size(pValues));
%       upperCI2 = readout2 + randi(numel(pValues), 10, 1) / 10;
%       lowerCI2 = readout2 - randi(numel(pValues), 10, 1) / 10;
%       readoutAll = [readout1, readout2];
%       upperCIAll = [upperCI1, upperCI2];
%       lowerCIAll = [lowerCI1, lowerCI2];
%
%       pValues = transpose(20:-2:2);
%       
%       plot_tuning_curve(pValues, readout1, 'UpperCI', upperCI1, 'LowerCI', lowerCI1);
%       plot_tuning_curve(pValues, readoutAll, 'UpperCI', upperCIAll, 'LowerCI', lowerCIAll, 'ColorMap', hsv(2));
%
%       plot_tuning_curve((1:3)', [rand(1, 10); 2*rand(1, 10); 5*rand(1, 10)]);
%       plot_tuning_curve((1:3)', [rand(1, 10); 2*rand(1, 10); 5*rand(1, 10)], 'RunTTest', true, 'RunRankTest', true);
%       plot_tuning_curve((1:3)', [randn(1, 10); 5 + randn(1, 10); 4 + randn(1, 10)], 'RunTTest', true, 'RunRankTest', true);
%       plot_tuning_curve(1:5, 2:6);
%       plot_tuning_curve(1:20, 20:-1:1, 'ColorMap', @lines);
%       plot_tuning_curve(1:100, randi(100, 1, 100), 'ColorMap', @parula);
%
% Outputs:
%       handles     - handles structure with fields:
%                       fig         - figure handle
%                       ax          - axes handle
%                       curves      - handles to the plotted tuning curves
%                       confInts    - (optional) handles to confidence intervals
%                       boundaries  - (optional) handles to boundary lines
%                       selected    - (optional) handles to selected value markers
%                       averages    - (optional) handles to phase average lines
%                       avgWindows  - (optional) handles to average window bars
%                   specified as a scalar structure%
% Arguments:
%       pValues     - vector(s) of parameter values
%                   must be a numeric 2-D array
%       readout     - vector(s) of readout values 
%                       each column is a readout vector
%                   must be a numeric 2-D array
%       varargin    - 'LowerCI': lower bounds of confidence intervals
%                   must be a numeric 2-D array
%                   default == []
%                   - 'UpperCI': upper bounds of confidence intervals
%                   must be a numeric 2-D array
%                   default == []
%                   - 'RemoveOutliers': whether to remove outliers
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunTTest': whether to run paired t-test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunRankTest': whether to run paired 
%                                       Wilcoxon signed-rank test
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ColumnsToPlot': columns of the readout matrix to plot
%                   must be a numeric vector
%                   default == 1:size(readout, 2);
%                   - 'LineSpec': line specification
%                   must be a character array
%                   default == '-'
%                   - 'LineWidth': line width
%                   must be a positive scalar
%                   default == 2
%                   - 'PIsLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ReadoutIsLog': whether readout values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PLimits': limits of parameter axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'ReadoutLimits': limits of readout axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == []
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'PTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
%                   - 'PLabel': label for the parameter, 
%                               suppress by setting value to {'suppress'}
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'ColumnLabels': labels for the readout columns, 
%                               suppress by setting value to {'suppress'}
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Column #1', 'Column #2', ...}
%                   - 'PhaseLabels': phase labels if phase vectors are provided
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Phase #1', 'Phase #2', ...}
%                   - 'ColorByPhase': whether to color by phase
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ColorMap' - color map used when nColumnsToPlot > 1
%                   must be a 2-D numeric array with 3 columns
%                   default == jet(nColumnsToPlot) or 
%                               rgb('SkyBlue') == [0.5273, 0.8047, 0.9180]
%                                   if nColumnsToPlot == 1 or
%                               hsv(maxNPhases) if phaseVectors is provided
%                   - 'ConfIntColorMap': color map for confidence intervals
%                   must be a 3-element vector
%                   default == WHITE - (WHITE - colorMap) * confIntFadePercentage;
%                   - 'SelectedColorMap': color map for selected values
%                   must be a 3-element vector
%                   default == Same as ColorMap
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nTraces == 1 
%                               'northeast' if nTraces is 2~9
%                               'eastoutside' if nTraces is 10+
%                   - 'PlotOnly': whether to plot the curves only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotForCorel': whether to plot for CorelDraw
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotPhaseBoundaries': whether to plot phase boundaries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotPhaseAverages': whether to plot phase averages
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotIndSelected': whether to plot selected indices
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotAverageWindows': whether to plot average windows
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PBoundaries': parameter boundary values
%                       Note: each row is a set of boundaries
%                   must be a numeric array
%                   default == []
%                   - 'PBoundaryType': type of parameter boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'verticalLines'     - vertical dotted lines
%                       'horizontalBars'    - horizontal bars
%                       'verticalShades'    - vertical shades
%                   default == 'verticalLines'
%                   - 'RBoundaries': readout boundary values
%                       Note: each row is a set of boundaries
%                   must be a numeric array
%                   default == []
%                   - 'RBoundaryType': type of readout boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'horizontalLines'   - horizontal dotted lines
%                       'verticalBars'      - vertical bars
%                       'horizontalShades'  - horizontal shades
%                   default == 'horizontalLines'
%                   - 'PhaseVectors': phase information for each readout vector
%                   must be a numeric matrix or a cell array of numeric vectors
%                   default == {}
%                   - 'PhaseBoundaries': vector of phase boundaries
%                   must be a numeric vector
%                   default == set in parse_phase_info.m
%                   - 'AverageWindows': windows to average values
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == set in parse_phase_info.m
%                   - 'PhaseAverages': average values for each phase
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric 2-D array
%                   default == set in parse_phase_info.m
%                   - 'IndSelected': selected indices to mark differently
%                       Note: If a matrix cell array, 
%                           each column is for a curve and each row is for a phase
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == set in parse_phase_info.m
%                   - 'NLastOfPhase': number of values at the last of a phase
%                   must be empty or a positive integer scalar
%                   default == set in parse_phase_info.m
%                   - 'NToAverage': number of values to average
%                   must be empty or a positive integer scalar
%                   default == set in select_similar_values.m
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'          - set in select_similar_values.m
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'auto'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be empty or a nonnegative scalar
%                   default == set in select_similar_values.m
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == set in set_figure_properties.m
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == set in set_figure_properties.m
%                   - 'FigExpansion': expansion factor for figure position
%                   must be a must be a positive scalar or 2-element vector
%                   default == set in set_figure_properties.m
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'AxesHandle': axes handle for created axes
%                   must be a empty or a axes object handle
%                   default == set in set_axes_properties.m
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == {'png', 'epsc'}
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       ~/Downloaded_Functions/rgb.m
%       ~/Adams_Functions/argfun.m
%       ~/Adams_Functions/cell2num.m
%       ~/Adams_Functions/count_samples.m
%       ~/Adams_Functions/create_error_for_nargin.m
%       ~/Adams_Functions/create_labels_from_numbers.m
%       ~/Adams_Functions/decide_on_colormap.m
%       ~/Adams_Functions/extract_fileparts.m
%       ~/Adams_Functions/fill_markers.m
%       ~/Adams_Functions/force_matrix.m
%       ~/Adams_Functions/force_row_vector.m
%       ~/Adams_Functions/hold_off.m
%       ~/Adams_Functions/hold_on.m
%       ~/Adams_Functions/isfigtype.m
%       ~/Adams_Functions/islegendlocation.m
%       ~/Adams_Functions/islog2scale.m
%       ~/Adams_Functions/match_column_count.m
%       ~/Adams_Functions/parse_phase_info.m
%       ~/Adams_Functions/plot_horizontal_line.m
%       ~/Adams_Functions/plot_selected.m
%       ~/Adams_Functions/plot_test_result.m
%       ~/Adams_Functions/plot_vertical_line.m
%       ~/Adams_Functions/plot_window_boundaries.m
%       ~/Adams_Functions/remove_outliers.m
%       ~/Adams_Functions/save_all_figtypes.m
%       ~/Adams_Functions/set_axes_properties.m
%       ~/Adams_Functions/set_default_flag.m
%       ~/Adams_Functions/set_figure_properties.m
%       ~/Adams_Functions/test_normality.m
%       ~/Adams_Functions/unique_custom.m
%       ~/Adams_Functions/union_over_cells.m
%       ~/Adams_Functions/update_figure_for_corel.m
%
% Used by:
%       ~/Adams_Functions/m3ha_network_tuning_curves.m
%       ~/Adams_Functions/m3ha_plot_figure08.m
%       ~/Adams_Functions/parse_current_family.m
%       ~/Adams_Functions/plot_calcium_imaging_traces.m
%       ~/Adams_Functions/plot_chevron.m
%       ~/Adams_Functions/plot_table_parallel.m
%       ~/Adams_Functions/plot_measures.m
%       ~/Adams_Functions/plot_struct.m
%       ~/Adams_Functions/plot_traces.m
%       ~/Adams_Functions/plot_table.m
%       ~/Matts_Functions/contour_plot.m
%       ~/Marks_Functions/Adam/CLC2/markCLC2figures.m
%       ~/Marks_Functions/Katie/twoDGtimeSeries.m
%       /media/adamX/RTCl/tuning_curves.m

% File History:
% 2017-04-17 Moved from tuning_curves.m
% 2017-04-17 Simplified code
% 2017-04-17 Set default arguments
% 2017-04-17 Color map is now based on number of columns to plot
% 2017-05-09 Added 'FigTypes' as a parameter-value pair argument
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-09-25 Made almost all arguments parameter-value pairs
% 2018-12-15 Added 'LineSpec' as a parameter-value pair argument
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Changed lineSpec default to o and singleColorDefault to SkyBlue
% 2019-03-14 Added 'RemoveOutliers' as an optional argument
% 2019-03-25 Added 'PhaseVectors' as an optional argument
% 2019-03-25 Now expands the y limits by a little by default
% 2019-05-10 Now uses set_figure_properties.m
% 2019-06-10 Added 'PBoundaries' and 'RBoundaries' as optional arguments
% 2019-08-07 Now changes the pTickAngle only if the labels are too long
% 2019-08-07 Added 'PhaseLabels' as an optional argument
% 2019-08-07 Added 'ColorMap' as an optional argument
% 2019-08-07 Added 'ClearFig' as an optional argument
% 2019-08-07 Now accepts infinite values for readout limits
% 2019-08-07 Added 'RunTTest' as an optional argument
% 2019-08-07 Added 'RunRankTest' as an optional argument
% 2019-08-09 Updated confidence interval plots
% 2019-08-09 Combined SingleColor with ColorMap
% 2019-08-09 Fixed confidence interval plots for matrices
% 2019-08-09 Added 'IndSelected' as an optional argument
% 2019-08-21 Now outputs a handles structure
% 2019-08-21 Added 'PlotPhaseBoundaries', 'PlotPhaseAverages', 'PlotIndSelected'
% 2019-08-22 Made averageWindows an optional argument
% 2019-08-27 Fixed usage of plot flags
% 2019-08-27 Added 'PlotAverageWindows'
% 2019-10-02 Added 'PlotOnly' as an optional argument
% 2019-10-02 Now plots a star if significant
% 2019-10-04 Now plots 'NS' if not significant
% 2019-10-07 Now plots 'NS' or star in black
% 2019-11-24 Moved code to parse_phase_info.m
% 2019-11-28 Now plots confidence intervals with transparency
% 2019-12-18 Now allows pValues to be multiple vectors
% 2019-12-23 Fixed colorMap argument
% 2019-12-23 Added 'ReadoutIsLog' as an optional argument
% 2020-02-17 Added normality tests
% 2020-02-19 Added 'PlotForCorel' as an optional argument
% 2025-08-29 Updated to use plot_test_result.m by Gemini
% 2025-10-07 Fixed bug in axHandle usage
% 2025-10-09 Updated default p limits to use average spacing
% TODO: Make 'FillMarkers' an optional argument with default == true
% TODO: Allow inputs to be cell arrays (use force_matrix.m)
% TODO: Use test_difference.m?
% TODO: phaseBoundaries needs to be provided into parse_phase_info.m

%% Hard-coded parameters
validSelectionMethods = {'auto', 'notNaN', 'maxRange2Mean'};
validPBoundaryTypes = {'verticalLines', 'horizontalBars', 'verticalShades'};
validRBoundaryTypes = {'horizontalLines', 'verticalBars', 'horizontalShades'};

% TODO: Make optional arguments
sigLevel = 0.05;                    % significance level for tests
confIntFadePercentage = 50;         % fade percentage for confidence interval colors
confIntLineStyle = 'none';
confIntFaceAlpha = 0.25;
confIntEdgeAlpha = 0.25;
selectedLineWidth = 3;              % line width for selected values markers
selectedMarker = 'o';
outlierMethod = 'fiveStds';
pBoundaryColor = '';                % set in plot_window_boundaries.m
pBoundaryLineStyle = '--';
pBoundaryLineWidth = 0.5;
rBoundaryColor = '';                % set in plot_window_boundaries.m
rBoundaryLineStyle = '--';
rBoundaryLineWidth = 0.5;
averagesLineStyle = ':';
averagesLineWidth = 2;
avgWindowRelYValue = 0.1;
avgWindowColorMap = [];
avgWindowLineStyle = '-';
avgWindowLineWidth = 3;
testXLocRel = 0.5;
starXLocRel = 0.5;
tTestPString = 'p_t';
tTestYLocText = 0.2;
tTestYLocStar = 0.9;
rankTestPString = 'p_r';
rankTestYLocText = 0.1;
rankTestYLocStar = 0.8;

%% Default values for optional arguments
lowerCIDefault = [];
upperCIDefault = [];
removeOutliersDefault = false;      % don't remove outliers by default
runTTestDefault = false;            % don't run paired t-test by default
runRankTestDefault = false;         % don't run paired signed-rank test by default
columnsToPlotDefault = [];          % set later
lineSpecDefault = '-';
lineWidthDefault = 2;
pIsLogDefault = false;
readoutIsLogDefault = false;
pLimitsDefault = [];
readoutLimitsDefault = [];
pTicksDefault = [];
pTickLabelsDefault = {};
pTickAngleDefault = [];             % set later
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
columnLabelsDefault = '';           % set later
phaseLabelsDefault = '';            % set later
colorByPhaseDefault = false;        % don't color by phase by default
colorMapDefault = [];               % set later
confIntColorMapDefault = [];        % set later
selectedColorMapDefault = [];       % set later
legendLocationDefault = 'auto';     % set later
plotOnlyDefault = false;            % setup default labels by default
plotForCorelDefault = false;
plotPhaseBoundariesDefault = [];    % set later
plotPhaseAveragesDefault = [];      % set later
plotIndSelectedDefault = [];        % set later
plotAverageWindowsDefault = [];     % set later
pBoundariesDefault = [];
pBoundaryTypeDefault = 'verticalLines';
rBoundariesDefault = [];
rBoundaryTypeDefault = 'horizontalLines';
phaseVectorsDefault = {};           % no phase vectors by default
phaseBoundariesDefault = [];        % set in parse_phase_info.m
averageWindowsDefault = {};         % set in parse_phase_info.m
phaseAveragesDefault = [];          % set in parse_phase_info.m
indSelectedDefault = [];            % set in parse_phase_info.m
nLastOfPhaseDefault = [];           % set in parse_phase_info.m
nToAverageDefault = [];             % set in select_similar_values.m
selectionMethodDefault = 'auto';    % set in select_similar_values.m
maxRange2MeanDefault = [];          % set in select_similar_values.m
figTitleDefault = '';               % set later
figHandleDefault = [];              % no existing figure by default
figNumberDefault = [];              % no figure number by default
figExpansionDefault = [];           % no figure expansion by default
clearFigureDefault = false;         % don't clear figure by default
axHandleDefault = [];               % gca by default
figNameDefault = '';                % don't save figure by default
figTypesDefault = {'png', 'epsc'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an Input Parser
addRequired(iP, 'pValues', ...              % vector of parameter values
    @(x) validateattributes(x, {'numeric', 'datetime', 'duration'}, ...
                                {'2d'}));
addRequired(iP, 'readout', ...              % a readout matrix
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LowerCI', lowerCIDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'UpperCI', upperCIDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'2d'}));
addParameter(iP, 'RemoveOutliers', removeOutliersDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunTTest', runTTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunRankTest', runRankTestDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ColumnsToPlot', columnsToPlotDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'PIsLog', pIsLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'ReadoutIsLog', readoutIsLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'PLimits', pLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ReadoutLimits', readoutLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PTickAngle', pTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColumnLabels', columnLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PhaseLabels', phaseLabelsDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorByPhase', colorByPhaseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'ConfIntColorMap', confIntColorMapDefault);
addParameter(iP, 'SelectedColorMap', selectedColorMapDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d', 'numel', 3}));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'PlotOnly', plotOnlyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotForCorel', plotForCorelDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotPhaseBoundaries', plotPhaseBoundariesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotPhaseAverages', plotPhaseAveragesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotIndSelected', plotIndSelectedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAverageWindows', plotAverageWindowsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PBoundaries', pBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PBoundaryType', pBoundaryTypeDefault, ...
    @(x) any(validatestring(x, validPBoundaryTypes)));
addParameter(iP, 'RBoundaries', rBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RBoundaryType', rBoundaryTypeDefault, ...
    @(x) any(validatestring(x, validRBoundaryTypes)));
addParameter(iP, 'PhaseVectors', phaseVectorsDefault, ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['PhaseVectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'PhaseBoundaries', phaseBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'AverageWindows', averageWindowsDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'PhaseAverages', phaseAveragesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'IndSelected', indSelectedDefault, ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));
addParameter(iP, 'NLastOfPhase', nLastOfPhaseDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NLastOfPhase must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'NToAverage', nToAverageDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NToAverage must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) assert(isempty(x) || isscalar(x) && x >= 0, ...
                ['MaxRange2Mean must be either empty ', ...
                    'or a nonnegative scalar!']));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigExpansion', figExpansionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive'}));
addParameter(iP, 'ClearFigure', clearFigureDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, pValues, readout, varargin{:});
lowerCI = iP.Results.LowerCI;
upperCI = iP.Results.UpperCI;
removeOutliers = iP.Results.RemoveOutliers;
runTTest = iP.Results.RunTTest;
runRankTest = iP.Results.RunRankTest;
columnsToPlot = iP.Results.ColumnsToPlot;
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
pIsLog = iP.Results.PIsLog;
readoutIsLog = iP.Results.ReadoutIsLog;
pLimits = iP.Results.PLimits;
readoutLimits = iP.Results.ReadoutLimits;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pTickAngle = iP.Results.PTickAngle;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
columnLabels = iP.Results.ColumnLabels;
phaseLabels = iP.Results.PhaseLabels;
colorByPhase = iP.Results.ColorByPhase;
colorMap = iP.Results.ColorMap;
confIntColorMap = iP.Results.ConfIntColorMap;
selectedColorMap = iP.Results.SelectedColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
plotOnly = iP.Results.PlotOnly;
plotForCorel = iP.Results.PlotForCorel;
plotPhaseBoundaries = iP.Results.PlotPhaseBoundaries;
plotPhaseAverages = iP.Results.PlotPhaseAverages;
plotIndSelected = iP.Results.PlotIndSelected;
plotAverageWindows = iP.Results.PlotAverageWindows;
pBoundaries = iP.Results.PBoundaries;
pBoundaryType = validatestring(iP.Results.PBoundaryType, validPBoundaryTypes);
rBoundaries = iP.Results.RBoundaries;
rBoundaryType = validatestring(iP.Results.RBoundaryType, validRBoundaryTypes);
phaseVectors = iP.Results.PhaseVectors;
phaseBoundaries = iP.Results.PhaseBoundaries;
averageWindows = iP.Results.AverageWindows;
phaseAverages = iP.Results.PhaseAverages;
indSelected = iP.Results.IndSelected;
nLastOfPhase = iP.Results.NLastOfPhase;
nToAverage = iP.Results.NToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;
clearFigure = iP.Results.ClearFigure;
axHandle = iP.Results.AxesHandle;
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Prepare for tuning curve
% Initialize a handles structure
handles = struct;

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
        numel(pTicks) ~= numel(pTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    return
end

% Reorder pTicks if not increasing and update pTickLabels accordingly
if ~isempty(pTicks) && ~issorted(pTicks)
    [pTicks, sortOrder] = sort(pTicks);
    if ~isempty(pTickLabels)
        pTickLabels = pTickLabels(sortOrder);
    end
end

% Display 
if isempty(phaseVectors) && colorByPhase
    fprintf('Warning: Cannot color by phase if phase vectors not provided!\n');    
    colorByPhase = false;
end

% If plotting curve only, change some defaults
if plotOnly
    pLabel = 'suppress';
    readoutLabel = 'suppress';
    figTitle = 'suppress';
    pLimits = 'suppress';
    readoutLimits = 'suppress';
    legendLocation = 'suppress';
end

% Count number of entries
nEntries = size(readout, 1);

% Count number of columns
nCols = size(readout, 2);

% Make sure pValues as the correct number of entries
if size(pValues, 1) ~= nEntries
    error('Number of rows in pValues does not match that of readout!');
end

% Get all unique parameter values in original order
uniquePValuesOrig = unique_custom(pValues, 'stable', 'IgnoreNaN', true);

% Get all unique parameter values in increasing order
uniquePValuesIncr = unique_custom(pValues, 'sorted', 'IgnoreNaN', true);

% Decide whether to plot phase-related stuff
[plotPhaseBoundaries, plotPhaseAverages, ...
    plotIndSelected, plotAverageWindows] = ...
    argfun(@(x) set_default_flag(x, ~isempty(phaseVectors) || ...
                                    ~isempty(phaseBoundaries)), ...
            plotPhaseBoundaries, plotPhaseAverages, ...
            plotIndSelected, plotAverageWindows);

% Decide on compute flags
[computePhaseInfo, computePhaseBoundaries, computeAverageWindows, ...
        computePhaseAverages, computeIndSelected] = ...
    argfun(@(x) set_default_flag([], x), ...
            ~isempty(phaseVectors) || ~isempty(phaseBoundaries), ...
            plotPhaseBoundaries && isempty(phaseBoundaries), ...
            (plotAverageWindows || plotPhaseAverages) && ...
                isempty(averageWindows), ...
            plotPhaseAverages && isempty(phaseAverages), ...
            plotIndSelected && isempty(indSelected));

% If phase-related stuff are to be computed, phase vectors 
%   or phase boundaries must be provided
if isempty(phaseVectors) && isempty(phaseBoundaries) && ...
        (computePhaseBoundaries || computeAverageWindows || ...
        computePhaseAverages || computeIndSelected)
    fprintf('Phase vectors must be provided!\n');
    return
end

% Compute phase-related stuff
if computePhaseInfo
    % Store arguments list
    parsePhaseInfoArgs = ...
        {pValues, readout, phaseVectors, ...
        'PhaseBoundaries', phaseBoundaries, ...
        'AverageWindows', averageWindows, ...
        'PhaseAverages', phaseAverages, ...
        'IndSelected', indSelected, ...
        'NLastOfPhase', nLastOfPhase, ...
        'NToAverage', nToAverage, ...
        'SelectionMethod', selectionMethod, ...
        'MaxRange2Mean', maxRange2Mean};

    % Parse phase-related stuff only if needed
    if computeIndSelected
        [uniquePhases, phaseBoundaries, averageWindows, ...
                phaseAverages, indSelected] = ...
            parse_phase_info(parsePhaseInfoArgs{:});
    elseif computePhaseAverages
        [uniquePhases, phaseBoundaries, averageWindows, phaseAverages] = ...
            parse_phase_info(parsePhaseInfoArgs{:});
    elseif computeAverageWindows
        [uniquePhases, phaseBoundaries, averageWindows] = ...
            parse_phase_info(parsePhaseInfoArgs{:});
    elseif computePhaseBoundaries
        [uniquePhases, phaseBoundaries] = ...
            parse_phase_info(parsePhaseInfoArgs{:});
    else
        uniquePhases = ...
            parse_phase_info(parsePhaseInfoArgs{:});
    end

    % Generate phase boundaries to be plotted if requested
    if plotPhaseBoundaries
        % Pool all phase boundaries together
        if iscell(phaseBoundaries)
            phaseBoundaries = union_over_cells(phaseBoundaries);
        end

        % Force as a row vector
        phaseBoundaries = force_row_vector(phaseBoundaries);

        % Add to parameter boundaries to plot
        pBoundaries = [pBoundaries; phaseBoundaries];
    end

    % Count the number of phases for each readout column
    nPhases = count_samples(uniquePhases);

    % Compute the maximum number of phases
    maxNPhases = max(nPhases);
end

% Count the number of boundaries
nPBoundaries = size(pBoundaries, 2);
nRBoundaries = size(rBoundaries, 2);

% Create default phase labels
if isempty(phaseLabels) && computePhaseInfo
    phaseLabels = create_labels_from_numbers(1:maxNPhases, 'Prefix', 'Phase #');
end

% Decide on the number of curves to plot per phase
if colorByPhase
    nLinesPerPhase = maxNPhases;
else
    nLinesPerPhase = 1;
end

% Remove outliers when plotting if requested
if removeOutliers
    [readoutToPlot, rowsToKeep] = ...
        remove_outliers(readout, 'OutlierMethod', outlierMethod, ...
                                    'ReplaceWithNans', true);
    pValuesToPlot = pValues(rowsToKeep, :);
else
    readoutToPlot = readout;
    pValuesToPlot = pValues;
end

% Decide on whether parametric or nonparametric tests are more appropriate
if (runTTest || runRankTest)
    if size(readout, 1) > 1
        % Transpose the readout matrix
        readoutTransposed = transpose(readout);

        % Extract the before columns
        befores = readoutTransposed(:, 1:end-1);

        % Extract the after columns
        afters = readoutTransposed(:, 2:end);

        % Compute consecutive differences
        diffData = afters - befores;

        % Test whether each pair of differences is normally-distributed
        isNormal = ...
            arrayfun(@(x) test_normality(diffData(:, x), ...
                                'SigLevel', sigLevel), ...
                    1:size(afters, 2));
    else
        disp('Difference test can''t be conducted for less than 2 groups!');
        runTTest = false;
        runRankTest = false;
    end
end

% Run paired t-tests if requested
if runTTest
    % Run t-tests on each pair of columns
    % TODO: Make a function vecfun.m
    [~, tTestPValues] = ...
        arrayfun(@(x) ttest(afters(:, x), befores(:, x), ...
                                'Alpha', sigLevel), ...
                    1:size(afters, 2));
else
    tTestPValues = [];
end

% Run paired Wilcoxon signed rank tests if requested
if runRankTest
    % Run t-tests on each pair of columns
    % TODO: Make a function vecfun.m
    rankTestPValues = ...
        arrayfun(@(x) signrank(afters(:, x), befores(:, x), ...
                                'Alpha', sigLevel), ...
                    1:size(afters, 2));
else
    rankTestPValues = [];
end

% Set default columns to plot
if isempty(columnsToPlot)
    columnsToPlot = transpose(1:nCols);
end

% Set column labels
if isempty(columnLabels)
    columnLabels = cell(1, nCols);
    for iPlot = 1:nCols
        columnLabels{iPlot} = ['Column #', num2str(iPlot)];
    end
end

% Set legend location based on number of traces
if strcmpi(legendLocation, 'auto')
    if nCols > 1 && nCols < 10
        legendLocation = 'northeast';
    elseif nCols >= 10
        legendLocation = 'eastoutside';
    else
        legendLocation = 'suppress';
    end
end

% Set the default figure title
if isempty(figTitle) && ~strcmpi(figTitle, 'suppress')
    if ~strcmpi(readoutLabel, 'suppress') && ~strcmpi(pLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. ', pLabel], '_', '\_');
    elseif ~strcmpi(readoutLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. parameter'], '_', '\_');
    else
        figTitle = 'Readout vs. parameter';
    end
end

% Compute the number of parameter values that 
%   don't give infinite values
nNonInf = sum(~isinf(readout), 1);

% Count the number of columns to plot
nColumnsToPlot = length(columnsToPlot);

% Decide on the color map to use
if colorByPhase
    % Use the hsv color map by default
    if isempty(colorMap)
        colorMap = @hsv;
    end
    
    % Update the color map based on maxNPhases
    colorMap = decide_on_colormap(colorMap, maxNPhases);
else
    % Use the jet color map by default
    if isempty(colorMap)
        if nColumnsToPlot == 1
            colorMap = rgb('SkyBlue');
        else
            colorMap = @jet;
        end
    end
    
    % Update the color map based on nColumnsToPlot
    colorMap = decide_on_colormap(colorMap, nColumnsToPlot);
end

% Decide on the confidence interval color map to use
if isempty(confIntColorMap)
    % Color of the confidence interval
    confIntColorMap = decide_on_colormap(colorMap, 'OriginalNColors', true, ...
                            'FadePercentage', confIntFadePercentage);
end

% Decide on the selected values color map to use
if isempty(selectedColorMap)
    selectedColorMap = colorMap;
end

% Set the default parameter tick angle
if isempty(pTickAngle)
    if ~isempty(pTickLabels)
        maxCharPTickLabels = max(count_samples(pTickLabels, ...
                                    'TreatCellStrAsArray', false));
        if maxCharPTickLabels > 10
            pTickAngle = 60;
        elseif maxCharPTickLabels > 3
            pTickAngle = 45;
        else
            pTickAngle = 0;
        end
    else
        pTickAngle = 0;
    end
end

%% Plot tuning curve
% Decide on the figure to plot on and set properties
fig = set_figure_properties('FigHandle', figHandle, 'AxesHandle', axHandle, ...
                    'FigNumber', figNumber, 'FigExpansion', figExpansion, ...
                    'ClearFigure', clearFigure);

% Decide on the axes to plot on and set properties
axHandle = set_axes_properties('AxesHandle', axHandle);

% Initialize graphics objects
curves = gobjects(nColumnsToPlot, nLinesPerPhase);
if ~isempty(lowerCI) || ~isempty(upperCI)
    confInts = gobjects(nColumnsToPlot, nLinesPerPhase);
end

% Hold on
wasHold = hold_on(axHandle);

% Plot readout values against parameter values
for iPlot = 1:nColumnsToPlot
    % Get the column to plot
    col = columnsToPlot(iPlot);
    if size(pValuesToPlot, 2) > 1
        pValuesThis = pValuesToPlot(:, col);
    else
        pValuesThis = pValuesToPlot;
    end
    readoutThis = readoutToPlot(:, col);

    % Plot the tuning curve for this column
    if colorByPhase       
        % Get the current phase vector
        phaseVectorThis = phaseVectors{iPlot};
        uniquePhasesThis = uniquePhases{iPlot};
        nPhasesThis = nPhases(iPlot);

        % Get the distinct phase indices for this readout vector
        phaseIndices = arrayfun(@(x) find(phaseVectorThis == x), ...
                                uniquePhasesThis, 'UniformOutput', false);

        % Get the last index for this readout vector
        lastIndex = numel(phaseVectorThis);

        % Add the next index
        phaseIndices = cellfun(@(x) add_next_index(x, lastIndex), ...
                                phaseIndices, 'UniformOutput', false);

        % Plot readout vector for all phases
        curves(iPlot, 1:nPhasesThis) = ...
            cellfun(@(x) plot_one_line(axHandle, pValuesThis(x), readoutThis(x), ...
                                    lineSpec, lineWidth, otherArguments), ...
                    phaseIndices);
    else
        curves(iPlot, 1) = plot_one_line(axHandle, pValuesThis, readoutThis, ...
                                    lineSpec, lineWidth, otherArguments);
    end

    % If provided, plot a confidence interval for this column
    %   as a light-gray-shaded area
    if ~isempty(lowerCI) || ~isempty(upperCI)
        if ~isempty(lowerCI)
            lowerCIThis = lowerCI(:, col);
        else
            lowerCIThis = readoutThis;
        end
        if ~isempty(upperCI)
            upperCIThis = upperCI(:, col);
        else
            upperCIThis = readoutThis;
        end

        % Plot the confidence interval
        if colorByPhase || pIsLog || readoutIsLog
            fprintf('Not Supported Yet!\n');
        else
            % Get the current Y limits
            % yLimits = get(axHandle, 'YLim');

            % Compute the minimum y limits
            %TODO
            % minY = apply_iteratively(@min, {yLimits, readoutLimits});

            % TODO: use plot_vertical_shade.m
            % Remove tuning curve
            delete(curves(iPlot, 1));

            % The x and y values for the confidence intervals
            confIntXValues = [pValuesThis; flipud(pValuesThis)];
            confIntYValues = [upperCIThis; flipud(lowerCIThis)];

            % Fill the area between lowerCIThis and upperCIThis 
            %   with confIntColorMap
            confInts(iPlot, 1) = fill(axHandle, confIntXValues, confIntYValues, ...
                                        confIntColorMap(iPlot, :), ...
                                        'LineStyle', confIntLineStyle, ...
                                        'FaceAlpha', confIntFaceAlpha, ...
                                        'EdgeAlpha', confIntEdgeAlpha);

            % Plot tuning curve again
            curves(iPlot, 1) = plot_one_line(axHandle, pValuesThis, readoutThis, ...
                                        lineSpec, lineWidth, otherArguments);

            % Display tick marks and grid lines over graphics objects.
            set(axHandle, 'Layer', 'top');
        end
    end

    % Set color
    if colorByPhase
        % Set color by phase
        arrayfun(@(x) set(curves(iPlot, x), 'Color', colorMap(x, :)), ...
                1:nPhasesThis);
    else
        % Set color by columns
        set(curves(iPlot, 1), 'Color', colorMap(iPlot, :));
    end

    % Set display name
    if colorByPhase
        arrayfun(@(x) set(curves(iPlot, x), 'DisplayName', phaseLabels{x}), ...
                1:nPhasesThis);
    else        
        if ~strcmpi(columnLabels, 'suppress')
            set(curves(iPlot, 1), 'DisplayName', columnLabels{col});
        end
    end

    % If there is only one value for this column, mark with a circle
    % TODO: for colorByPhase?
    if nNonInf(col) == 1
        set(curves(iPlot, 1), 'Marker', 'o');
    end
end

% Fill markers if any
fill_markers('AxesHandle', axHandle);

% Restrict x axis if pLimits provided; 
%   otherwise expand the x axis by a little bit
if ~isempty(pLimits)
    if ~strcmpi(pLimits, 'suppress')
        % Use x limits
        xlim(axHandle, pLimits);
    else
        % Do nothing if user requested 'suppress'
    end
elseif numel(uniquePValuesIncr) > 1
    % Compute the average spacing between parameter values
    avgPSpacing = mean(diff(uniquePValuesIncr));

    % Set the new limits by extending by the average spacing
    xlim(axHandle, [uniquePValuesIncr(1) - avgPSpacing/2, uniquePValuesIncr(end) + avgPSpacing/2]);
else
    % Do nothing if there are less than 2 entries
end

% Restrict y axis if readoutLimits provided
%   otherwise expand the y axis by a little bit
if ~isempty(readoutLimits) && ~strcmpi(readoutLimits, 'suppress')
    yLimitsOrig = get(axHandle, 'YLim');
    if isinf(readoutLimits(1))
        readoutLimits(1) = yLimitsOrig(1);
    end
    if isinf(readoutLimits(2))
        readoutLimits(2) = yLimitsOrig(2);
    end
    ylim(axHandle, readoutLimits);
elseif ~strcmpi(readoutLimits, 'suppress')
    ylim(axHandle, compute_axis_limits(get(axHandle, 'YLim'), 'y', 'Coverage', 80));
end

% Set title and axes labels
if ~isempty(pTicks)
    xticks(axHandle, pTicks);
end
if ~isempty(pTickLabels)
    xticklabels(axHandle, pTickLabels);
end
if ~strcmpi(pLabel, 'suppress')
    xlabel(axHandle, pLabel);
end
if ~strcmpi(readoutLabel, 'suppress')
    ylabel(axHandle, readoutLabel);
end
if ~strcmpi(figTitle, 'suppress')
    title(axHandle, figTitle);
end

% Set the angle for parameter ticks
xtickangle(axHandle, pTickAngle);

% Plot parameter boundaries
if nPBoundaries > 0
    pLines = plot_window_boundaries(pBoundaries, ...
                                'BoundaryType', pBoundaryType, ...
                                'LineWidth', pBoundaryLineWidth, ...
                                'LineStyle', pBoundaryLineStyle, ...
                                'ColorMap', pBoundaryColor, ...
                                'AxesHandle', axHandle);
else
    pLines = gobjects;
end

% Plot readout boundaries
if nRBoundaries > 0
    rLines = plot_window_boundaries(rBoundaries, ...
                                'BoundaryType', rBoundaryType, ...
                                'LineWidth', rBoundaryLineWidth, ...
                                'LineStyle', rBoundaryLineStyle, ...
                                'ColorMap', rBoundaryColor, ...
                                'AxesHandle', axHandle);
else
    rLines = gobjects;
end

% Plot phaseAverages if any
if plotPhaseAverages && ~isempty(phaseAverages) && ~isempty(averageWindows)
    % Decide on the color map for each phase
    if colorByPhase
        colorMapEachPhase = ...
            arrayfun(@(x) repmat(colorMap(x, :), nColumnsToPlot, 1), ...
                        1:size(colorMap, 1), 'UniformOutput', false);
    else
        colorMapEachPhase = repmat({colorMap}, maxNPhases, 1);
    end

    % Plot the phase averages as horizontal lines
    averages = ...
        arrayfun(@(x) ...
            arrayfun(@(y) plot_horizontal_line(phaseAverages(x, columnsToPlot(y)), ...
                                'XLimits', averageWindows{x, columnsToPlot(y)}, ...
                                'ColorMap', colorMapEachPhase{x}(y, :), ...
                                'LineStyle', averagesLineStyle, ...
                                'LineWidth', averagesLineWidth, ...
                                'AxesHandle', axHandle), ...
                    transpose(1:nColumnsToPlot)), ...
            transpose(1:size(phaseAverages, 1)), 'UniformOutput', false);
            
    % Force the graphics array as a matrix
    %   Note: Each column is curve and each row is a phase
    averages = transpose(force_matrix(averages));
end

% Plot selected values if any
if plotIndSelected && ~isempty(indSelected)
    if iscell(indSelected)
        % Match the column count
        nColumns = size(readoutToPlot, 2);
        pValuesToPlot = match_column_count(pValuesToPlot, nColumns);
        
        % Color arbitrarily first
        selectedCell = ...
            arrayfun(@(x) ...
                cellfun(@(y) plot_selected(...
                            pValuesToPlot(:, columnsToPlot(x)), ...
                            readoutToPlot(:, columnsToPlot(x)), y, ...
                            'Marker', selectedMarker, 'ColorMap', 'r', ...
                            'LineWidth', selectedLineWidth, ...
                            'AxesHandle', axHandle), ...
                        indSelected(:, columnsToPlot(x))), ...
                1:nColumnsToPlot, 'UniformOutput', false);            

        % Force the graphics array as a matrix
        %   Note: Each column is curve and each row is a phase
        selected = force_matrix(selectedCell);

        % Change the color 
        if colorByPhase
            for iPhase = 1:maxNPhases
                colorThis = selectedColorMap(iPhase, :);

                for iCol = 1:nColumnsToPlot
                    x = selected(iPhase, iCol);
                    if isgraphics(x)
                        x.Color = colorThis;
                    end
                end
            end
        else
            for iCol = 1:nColumnsToPlot
                colorThis = selectedColorMap(iCol, :);
                for iPhase = 1:maxNPhases
                    x = selected(iPhase, iCol);
                    if isgraphics(x)
                        x.Color = colorThis;
                    end
                end
            end
        end
    else
        selected = plot_selected(pValuesToPlot, readoutToPlot, indSelected, ...
                                'Marker', selectedMarker, ...
                                'ColorMap', selectedColorMap(1, :), ...
                                'LineWidth', selectedLineWidth, ...
                                'AxesHandle', axHandle);
    end
end

% Plot averageWindows if requested
if plotAverageWindows && ~isempty(averageWindows)
    % Decide on the color map
    avgWindowColorMap = decide_on_colormap(avgWindowColorMap, maxNPhases, ...
                                            'ColorMapFunc', @hsv);

    % Plot the average windows as horizontal lines
    avgWindows = ...
        arrayfun(@(x) plot_window_boundaries(averageWindows{x, 1}, ...
                                'BarRelValue', avgWindowRelYValue, ...
                                'BoundaryType', 'horizontalBar', ...
                                'ColorMap', avgWindowColorMap(x, :), ...
                                'LineStyle', avgWindowLineStyle, ...
                                'LineWidth', avgWindowLineWidth, ...
                                'AxesHandle', axHandle), ...
                transpose(1:maxNPhases), 'UniformOutput', false);
end

% Decide on the x locations for the test results
if ~isempty(tTestPValues) || ~isempty(rankTestPValues)
    xIntervals = diff(uniquePValuesOrig);
    xLocText = uniquePValuesOrig(1:end-1) + xIntervals .* testXLocRel;
    xLocStar = uniquePValuesOrig(1:end-1) + xIntervals .* starXLocRel;
end

% Plot t-test p values if any
if ~isempty(tTestPValues)
    plot_test_result(tTestPValues, 'PString', tTestPString, ...
                    'YLocTextRel', tTestYLocText, 'YLocStarRel', tTestYLocStar, ...
                    'XLocText', xLocText, 'XLocStar', xLocStar, ...
                    'SigLevel', sigLevel, 'IsAppropriate', isNormal, ...
                    'AxesHandle', axHandle);
end

% Plot rank test p values if any
if ~isempty(rankTestPValues)
    plot_test_result(rankTestPValues, 'PString', rankTestPString, ...
                    'YLocTextRel', rankTestYLocText, 'YLocStarRel', rankTestYLocStar, ...
                    'XLocText', xLocText, 'XLocStar', xLocStar, ...
                    'SigLevel', sigLevel, 'IsAppropriate', ~isNormal, ...
                    'AxesHandle', axHandle);
end

% Hold off
hold_off(wasHold, axHandle);

% Modify axes scale if requested
%   Note: this must occur after holding off
if pIsLog || readoutIsLog
    [xScale, yScale] = argfun(@islog2scale, pIsLog, readoutIsLog);
    set(axHandle, 'XScale', xScale, 'YScale', yScale);
end

%% Post-plotting
% Generate a legend for the curves only if there is more than one trace
if ~strcmpi(legendLocation, 'suppress') && nColumnsToPlot > 1
    if colorByPhase
        lgd = legend(axHandle, curves(1, :), 'location', legendLocation);
    else
        lgd = legend(axHandle, curves, 'location', legendLocation);
    end

    set(lgd, 'AutoUpdate', 'off', 'Interpreter', 'none');
end

% Save figure if figName provided
if ~isempty(figName)
    % Plot for CorelDraw
    if plotForCorel
        % Create path for original figure
        figNameOrig = [extract_fileparts(figName, 'pathbase'), '_orig'];

        % Save original figure as png
        save_all_figtypes(fig, figNameOrig, 'png');

        % Update figure for CorelDraw
        update_figure_for_corel(fig);
    end

    % Save figure
    save_all_figtypes(fig, figName, figTypes);
end

%% Output handles
handles.fig = fig;
handles.ax = axHandle;
handles.curves = curves;
if ~isempty(lowerCI) || ~isempty(upperCI)
    handles.confInts = confInts;
end
if ~isempty(pLines) || ~isempty(rLines)
    handles.boundaries = transpose(vertcat(pLines, rLines));
end
if plotIndSelected && ~isempty(indSelected)
    handles.selected = selected;
end
if plotPhaseAverages && ~isempty(phaseAverages) && ~isempty(averageWindows)
    handles.averages = averages;
end
if plotAverageWindows && ~isempty(averageWindows)
    handles.avgWindows = avgWindows;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = plot_one_line(axHandle, pValues, readout, lineSpec, ...
                            lineWidth, otherArguments)

p = plot(axHandle, pValues, readout, lineSpec, ...
                'LineWidth', lineWidth, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = add_next_index(indices, lastIndex)
%% Add the next index if not the last
% TODO: What if indices not contiguous?

if indices(end) < lastIndex
    indices = [indices; indices(end) + 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Usage: plot_tuning_curve(pValues, readout, columnsToPlot, pIsLog, pLabel, ...
            readoutLabel, columnLabels, pLimits, readoutLimits, figName, varargin)

if ~isequal(columnLabels, {'suppress'})

if isequal(pLimits, -1)

if ~isequal(pLabel, 'suppress')
if ~isequal(readoutLabel, 'suppress')
if ~isequal(pLabel, 'suppress') && ~isequal(readoutLabel, 'suppress')

singleColorDefault = [0, 0, 1];
lineSpecDefault = '-';

set(fig, 'Visible', 'Off');
fig = figure(floor(rand()*10^4)+1);

if pIsLog
    % Note: can't have hold on before semilogx
    p = semilogx(pValues, readout(:, col), lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
else
    p = plot(pValues, readout(:, col), lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
end

if ~isempty(figName)
    % Create a figure
    if ~isempty(figNumber)
        % Create an invisible figure
        fig = figure(figNumber);
        set(fig, 'Visible', 'Off');
    else
        % Create a new figure
        fig = figure;
    end

    % Clear the figure
    clf(fig);
else
    % Get the current figure
    fig = gcf;
end

set(axHandle, 'XTick', pTicks);
set(axHandle, 'XTickLabel', pTickLabels);
pTickAngle = 60;                % x tick angle in degrees

phaseVectorsNoNaN = cellfun(@(x) x(~isnan(x)), phaseVectors, ...
                            'UniformOutput', false);
uniquePhases = cellfun(@(x) unique(x, 'stable'), phaseVectorsNoNaN, ...
                        'UniformOutput', false);

confInts = gobjects(nColumnsToPlot, 2);

% Make the area under upperCIThis light gray
confInts(iPlot, 1) = area(pValues, upperCIThis, minY, ...
                    'LineStyle', 'none', 'FaceColor', [0.9, 0.9, 0.9]);

% Make the area under lowerCIThis white
confInts(iPlot, 2) = area(pValues, lowerCIThis, minY, ...
                    'LineStyle', 'none', 'FaceColor', [1, 1, 1]);

%                   - 'SingleColor': color when nColumnsToPlot == 1
%                   must be a 3-element vector
%                   default == rgb('SkyBlue') == [0.5273, 0.8047, 0.9180]
singleColorDefault = rgb('SkyBlue');
addParameter(iP, 'SingleColor', singleColorDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3}));
singlecolor = iP.Results.SingleColor;
if nColumnsToPlot > 1
    set(curves(iPlot, 1), 'Color', colorMap(iPlot, :));
elseif nColumnsToPlot == 1
    set(curves(iPlot, 1), 'Color', singlecolor);
end

confIntColorMapDefault = rgb('LightGray');  
                                    % light gray confidence intervals by default

if colorByPhase
    selectedCell = ...
        arrayfun(@(x) ...
            cellfun(@(y) plot_selected(pValues, ...
                        readoutToPlot(:, columnsToPlot(x)), y, ...
                        selectedMarker, selectedColorMap(y, :), ...
                        selectedLineWidth), ...
                    indSelected(:, columnsToPlot(x))), ...
            1:nColumnsToPlot, 'UniformOutput', false);
else
    selectedCell = ...
        arrayfun(@(x) ...
            cellfun(@(y) plot_selected(pValues, ...
                        readoutToPlot(:, columnsToPlot(x)), y, ...
                        selectedMarker, selectedColorMap(x, :), ...
                        selectedLineWidth), ...
                    indSelected(:, columnsToPlot(x))), ...
            1:nColumnsToPlot, 'UniformOutput', false);            
end


% Decide on the average window y value
if isempty(avgWindowYValue)
    % Get the current y axis limits
    yLimitsNow = get(axHandle, 'YLim');

    % Compute a default window bar y value
    avgWindowYValue = yLimitsNow(1) + 0.1 * (yLimitsNow(2) - yLimitsNow(1));
end
avgWindows = ...
    arrayfun(@(x) plot_horizontal_line(avgWindowYValue, ...
                            'XLimits', averageWindows{x, 1}, ...
                            'ColorMap', avgWindowColorMap(x, :), ...
                            'LineStyle', avgWindowLineStyle, ...
                            'LineWidth', avgWindowLineWidth), ...
            transpose(1:maxNPhases), 'UniformOutput', false);


pLines = plot_vertical_line(pBoundaries, 'LineWidth', 0.5, ...
                            'LineStyle', pBoundaryLineStyle, 'Color', 'g');
rLines = plot_horizontal_line(rBoundaries, 'LineWidth', 0.5, ...
                            'LineStyle', rBoundaryLineStyle, 'Color', 'r');

% Hold on if more than one column
if nColumnsToPlot > 1
    hold on
end
        hold on;
            hold on;
    hold on
    hold on
% Hold off if more than one column
if nColumnsToPlot > 1
    hold off
end

for iPhase = 1:nPhasesThis
    set(curves(iPlot, iPhase), 'Color', colorMap(iPhase, :));
end

set(curves(iPlot, iPhase), 'DisplayName', ...
    replace(phaseLabels{iPhase}, '_', '\_'));

set(curves(iPlot, 1), 'DisplayName', ...
    replace(columnLabels{col}, '_', '\_'));

nBoundaries = nPBoundaries + nRBoundaries;

% Note: can't have hold on before loglog, semilogx or semilogy
if pIsLog && readoutIsLog
    p = loglog(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
elseif pIsLog && ~readoutIsLog
    p = semilogx(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
elseif ~pIsLog && readoutIsLog
    p = semilogy(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
else
    p = plot(pValues, readout, lineSpec, ...
                    'LineWidth', lineWidth, otherArguments);
end

if pIsLog && readoutIsLog
    set(axHandle, 'XScale', 'log', 'YScale', 'log');
elseif pIsLog && ~readoutIsLog
    set(axHandle, 'XScale', 'log');
elseif ~pIsLog && readoutIsLog
    set(axHandle, 'YScale', 'log');
end

%% Hard-coded constants
WHITE = [1, 1, 1];
confIntColorMap = WHITE - (WHITE - colorMap) * confIntFadePercentage;

function plot_test_result (testPValues, pString, yLocTextRel, yLocStarRel, ...
                            xLocTextRel, xLocStarRel, uniquePValuesOrig, ...
                            sigLevel, isAppropriate)
%% Plots p values and star if significant

% Decide on the x locations
xLocText = uniquePValuesOrig(1:end-1) + (uniquePValuesOrig(2) - uniquePValuesOrig(1)) * xLocTextRel;
xLocStar = uniquePValuesOrig(1:end-1) + (uniquePValuesOrig(2) - uniquePValuesOrig(1)) * xLocStarRel;

% Get current y axis limits
yLimitsNow = get(axHandle, 'YLim');

% Decide on the y location for texts
yLocText = yLimitsNow(1) + (yLimitsNow(2) - yLimitsNow(1)) * yLocTextRel;

% Decide on the y location for stars
yLocStar = yLimitsNow(1) + (yLimitsNow(2) - yLimitsNow(1)) * yLocStarRel;

% Plot texts
for iValue =  1:numel(testPValues)
    % Get the current values
    testPValueThis = testPValues(iValue);
    xLocTextThis = xLocText(iValue);
    xLocStarThis = xLocStar(iValue);
    isAppropriateThis = isAppropriate(iValue);

    % Create a p value string to 2 significant digits
    pValueString = [pString, ' = ', num2str(testPValueThis, 2)];

    % Plot gray if inappropriate, red if significant
    if ~isAppropriateThis
        pColor = [0.5, 0.5, 0.5];       % rgb('Gray')
    elseif testPValueThis < sigLevel
        pColor = 'r';
    else
        pColor = 'k';
    end

    % Plot text
    % TODO: Make function plot_text.m
    text(xLocTextThis, yLocText, pValueString, 'Color', pColor, ...
            'HorizontalAlignment', 'center');

    % Plot star if significant, 'NS' if not
    if testPValueThis < sigLevel
        plot(xLocStarThis, yLocStar, '*', 'Color', [0, 0, 0], ...
            'MarkerSize', 4);
    else
        text(xLocStarThis, yLocStar, 'NS', 'Color', [0, 0, 0], ...
            'HorizontalAlignment', 'center');
    end
end

% Restrict x axis if pLimits provided; 
%   otherwise expand the x axis by a little bit
if ~isempty(pLimits)
    if ~strcmpi(pLimits, 'suppress')
        % Use x limits
        xlim(pLimits);
    end
else
    if nEntries > 1 && nEntries < 4
        xlim(compute_axis_limits(uniquePValuesOrig, 'x', 'Coverage', 90));
    elseif nEntries >= 4
        % Compute the average spacing between parameter values
        avgPSpacing = mean(diff(uniquePValuesOrig));

        % Set the new limits by extending by the average spacing
        xlim([uniquePValuesOrig(1) - avgPSpacing/2, uniquePValuesOrig(end) + avgPSpacing/2]);
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
