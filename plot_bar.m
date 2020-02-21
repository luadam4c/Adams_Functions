function handles = plot_bar (val, varargin)
%% Plots a bar graph (grouped or not) with or without confidence intervals
% Usage: handles = plot_bar (val, varargin)
% Explanation:
%       TODO
%
% Example:
%       load_examples;
%       handles = plot_bar(val1, low1, high1);
%       handles = plot_bar(val2, low2, high2);
%       handles = plot_bar(val2, low2, high2, 'ForceVectorAsRow', true);
%       handles = plot_bar(val2, low2, high2, 'BarDirection', 'horizontal');
%       handles = plot_bar(val1, low1, high1, 'BarDirection', 'horizontal', 'ForceVectorAsRow', false);
%       handles = plot_bar(val3, low3, high3, 'PTickLabels', {'Mark', 'Ashley', 'Katie', 'Adam'});
%       handles = plot_bar(val3, low3, high3, 'BarDirection', 'horizontal', 'PTickLabels', {'Mark', 'Ashley', 'Katie', 'Adam'});
%       handles = plot_bar(val3, low3, high3, 'ReverseOrder', true, 'PTickLabels', {'Mark', 'Ashley', 'Katie', 'Adam'});
%       handles = plot_bar(val1, low1, high1, 'ReverseOrder', true, 'BarDirection', 'horizontal');
%       plot_bar([1, 2], 'ColorMap', {'Red', 'Green'})
%       testValues = (20:-1:1);
%       plot_bar(testValues, 'ReverseOrder', true, 'BarDirection', 'horizontal', 'PTickLabels', create_labels_from_numbers(testValues))
%
% Outputs:
%       handles - handles structure with fields:
%                   bars    - bar objects (bars for each group or 
%                                           column is one bar object)
%                           specified as a vector of Bar object handles
%                   lines   - error bar lines
%                               1st dim: connecting (1), upper limit (2), lower limit (3)
%                               2nd dim: sample number
%                               3rd dim: group number
%                           specified as an array of Primitive Line object handles
%                   fig     - figure handle
%                           specified as a Figure object handle
%               specified as a structure
%
% Arguments:
%       val     - mean values for the bar() function
%                   each column is a different group
%                   each row is a different sample number
%               must be a numeric array accepted by the bar() function
%       low     - (opt) lower limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       high    - (opt) upper limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       varargin    - 'ForceVectorAsRow': whether to force a vector 
%                                       as a row vector (multiple groups)
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ReverseOrder': whether to reverse the order of the traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'BarDirection': bar direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'vertical'   - vertical bars
%                       'horizontal' - horizontal bars
%                   default == 'vertical'
%                   - 'BarWidth': relative bar width
%                   must be a numeric scalar
%                   default == 0.8
%                   - 'GroupStyle': bar group style
%                   must be recognizable as the 'style' parameter of the 
%                       built-in bar() function
%                   default == 'grouped'
%                   - 'CIRelativeBarWidth': TODO
%                   - 'CIBarWidth': TODO
%                   - 'CILineWidth': TODO
%                   - 'CIColor': TODO
%                   - 'PValues': parameter axis values for bar placement
%                   must be a numeric array
%                   default == 1:nCols
%                   - 'PLimits': limits of parameter axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - 'ReadoutLimits': limits of readout axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == []
%                   - 'PTicks': tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'PTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'ColumnLabels': labels for the readout columns (groups), 
%                               suppress by setting value to {'suppress'}
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Column #1', 'Column #2', ...}
%                   - 'PhaseLabels': phase labels if phase vectors are provided TODO
%                   must be a scalartext 
%                       or a cell array of strings or character vectors
%                   default == {'Phase #1', 'Phase #2', ...}
%                   - 'ColorMap' - color map used when nColumnsToPlot > 1
%                   must be a 2-D numeric array with 3 columns
%                   default == set by the bar() function
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nGroups == 1 or oneSamplePerGroup
%                               'northeast' if nGroups is 2~9
%                               'eastoutside' if nGroups is 10+
%                   - 'PlotOnly': whether to plot the bars only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotPhaseBoundaries': whether to plot phase boundaries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotPhaseAverages': whether to plot phase averages TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotIndSelected': whether to plot selected indices TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PlotAverageWindows': whether to plot average windows TODO
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true if PhaseVectors provided, false otherwise
%                   - 'PBoundaries': parameter boundary values TODO
%                       Note: each row is a set of boundaries
%                   must be a numeric array
%                   default == []
%                   - 'PBoundaryType': type of parameter boundaries TODO
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'horizontalLines'   - horizontal dotted lines
%                       'verticalBars'      - vertical bars
%                       'horizontalShades'  - horizontal shades
%                   default == 'horizontalLines'
%                   - 'RBoundaries': readout boundary values
%                       Note: each row is a set of boundaries
%                   must be a numeric array
%                   default == []
%                   - 'RBoundaryType': type of readout boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'verticalLines'     - vertical dotted lines
%                       'horizontalBars'    - horizontal bars
%                       'verticalShades'    - vertical shades
%                   default == 'verticalLines'
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
%                   default == TODO: ['Bar graph for ', figName]
%                               or [readoutLabel, ' across ', pLabel]
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the bar() function
%   
%
% Requires:
%       cd/argfun.m
%       cd/decide_on_colormap.m
%       cd/set_figure_properties.m
%       cd/force_column_vector.m
%       cd/isfigtype.m
%       cd/islegendlocation.m
%       cd/plot_horizontal_line.m
%       cd/plot_selected.m
%       cd/plot_vertical_line.m
%       cd/save_all_figtypes.m
%       cd/struct2arglist.m
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:
%       cd/m3ha_compute_and_compare_lts_statistics.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_rank_neurons.m
%       cd/parse_multiunit.m
%       cd/plot_chevron_bar_inset.m
%       cd/plot_struct.m
%       cd/ZG_fit_IEI_distributions.m
%       /media/adamX/Paula_IEIs/paula_iei3.m

% File History: 
% 2017-10-19 Moved from paula_iei3.m
% 2017-10-19 Added input parser and various optional arguments
% 2017-12-01 Added figure(h)
% 2018-03-08 Changed figure(h) to set(0, 'CurrentFigure', h)
%               to prevent the display of invisible figures
% 2019-01-15 Renamed bar_w_CI.m -> plot_bar.m
% 2019-01-15 Added otherArguments
% 2019-01-15 Made h -> 'FigHandle' an optional argument
% 2019-05-08 Added 'barDirection' as an optional argument
% 2019-05-08 Made low and high optional arguments
% 2019-05-10 Added 'ForceVectorAsRow' as an optional argument
% 2019-05-10 Fixed bugs in logic
% 2019-05-10 Fix error bars when 'barDirection' is 'horizontal'
% 2019-05-10 Now uses plot_error_bar.m
% 2019-05-10 Now grabs XOffset or YOffset and uses it to plot error bars
% 2019-05-10 Now uses set_figure_properties.m
% 2019-05-11 Added 'ReverseOrder' as an optional argument
% 2019-05-11 Added 'FigTitle' as an optional argument
% 2019-06-10 Added 'PBoundaries' and 'RBoundaries' as optional arguments
% 2019-11-23 Added 'PhaseVectors' and other dependent optional arguments
% 2019-11-24 Now returns handles structure as output
% 2019-12-03 Added 'ColorMap' as an optional argument
% 2019-12-18 Added 'GroupStyle' as an optional argument
% 2019-12-18 Added 'LegendLocation' as an optional argument
% 2019-12-18 Added 'ColumnLabels' as an optional argument
% 2020-02-07 Added 'BarWidth' as an optional argument
% TODO: Make sure there are enough ticks for the labels!
% TODO: phaseBoundaries needs to be provided into parse_phase_info.m
% TODO: Finish implementation of 'PhaseVectors' as in plot_tuning_curve
% TODO: Change usage in all functions using this
% TODO: Expand to use bar3()

%% Hard-coded parameters
validBarDirections = {'vertical', 'horizontal'};
validSelectionMethods = {'auto', 'notNaN', 'maxRange2Mean'};
validPBoundaryTypes = {'horizontalLines', 'verticalBars', 'horizontalShades'};
validRBoundaryTypes = {'verticalLines', 'horizontalBars', 'verticalShades'};
averagesLineStyle = ':';
averagesLineWidth = 2;

%% Default values for optional arguments
lowDefault = [];
highDefault = [];
forceVectorAsRowDefault = false;
reverseOrderDefault = false;        % don't reverse order by default
barDirectionDefault = 'vertical';
barWidthDefault = 0.8;
groupStyleDefault = 'grouped';
cIRelativeBarWidthDefault = 0.33;   % default bar width relative to separation
cIBarWidthDefault = [];
cILineWidthDefault = 2;             % default line width for CIs
cIColorDefault = '';
pValuesDefault = [];
pLimitsDefault = [];
readoutLimitsDefault = [];
pTicksDefault = [];
pTickLabelsDefault = {};
pTickAngleDefault = [];
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
columnLabelsDefault = '';           % set later
phaseLabelsDefault = '';            % set later
colorMapDefault = [];               % set later
legendLocationDefault = 'auto';     % set later
plotOnlyDefault = false;            % setup default labels by default
plotPhaseBoundariesDefault = [];    % set later
plotPhaseAveragesDefault = [];      % set later
plotIndSelectedDefault = [];        % set later
plotAverageWindowsDefault = [];     % set later
pBoundariesDefault = [];
pBoundaryTypeDefault = 'horizontalLines';
rBoundariesDefault = [];
rBoundaryTypeDefault = 'verticalLines';
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
figNameDefault = '';                % don't save figure by default
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'val', ...                      % values
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'low', lowDefault, ...          % low limit of CI
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addOptional(iP, 'high', highDefault, ...        % high limit of CI
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceVectorAsRow', forceVectorAsRowDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ReverseOrder', reverseOrderDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BarDirection', barDirectionDefault, ...
    @(x) any(validatestring(x, validBarDirections)));
addParameter(iP, 'BarWidth', barWidthDefault);
addParameter(iP, 'GroupStyle', groupStyleDefault);
addParameter(iP, 'CIRelativeBarWidth', cIRelativeBarWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CIBarWidth', cIBarWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CILineWidth', cILineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CIColor', cIColorDefault);
addParameter(iP, 'PValues', pValuesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PLimits', pLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ReadoutLimits', readoutLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
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
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
addParameter(iP, 'PlotOnly', plotOnlyDefault, ...
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
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, val, varargin{:});
low = iP.Results.low;
high = iP.Results.high;
forceVectorAsRow = iP.Results.ForceVectorAsRow;
reverseOrder = iP.Results.ReverseOrder;
barDirection = validatestring(iP.Results.BarDirection, validBarDirections);
barWidth = iP.Results.BarWidth;
groupStyle = iP.Results.GroupStyle;
cIRelativeBarWidth = iP.Results.CIRelativeBarWidth;
cIBarWidth = iP.Results.CIBarWidth;
cILineWidth = iP.Results.CILineWidth;
cIColor = iP.Results.CIColor;
pValues = iP.Results.PValues;
pLimits = iP.Results.PLimits;
readoutLimits = iP.Results.ReadoutLimits;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pTickAngle = iP.Results.PTickAngle;
pLabel = iP.Results.PLabel;
readoutLabel = iP.Results.ReadoutLabel;
columnLabels = iP.Results.ColumnLabels;
phaseLabels = iP.Results.PhaseLabels;
colorMap = iP.Results.ColorMap;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
plotOnly = iP.Results.PlotOnly;
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
figName = iP.Results.FigName;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the bar() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Initialize output
handles = struct;

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
        numel(pTicks) ~= numel(pTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    return
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

% Either force all vectors as row vectors
%   or make the vectors consistent
%   Note: This will cause each value of a vector to be plotted as separately
%           colored bars
if forceVectorAsRow
    [val, low, high] = ...
        argfun(@(x) force_row_vector(x, 'IgnoreNonVectors', true), ...
                val, low, high);
else
    if iscolumn(val)
        [low, high] = ...
            argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
                    low, high);
    elseif isrow(val)
        [low, high] = ...
            argfun(@(x) force_row_vector(x, 'IgnoreNonVectors', true), ...
                    low, high);
    end
end

% Count the number of rows (samples)
nRows = size(val, 1);

% Count the number of columns (groups)
nCols = size(val, 2);

% Decide whether there is only one row
singleGroup = nCols == 1;

% Decide whether there is one sample per group
oneSamplePerGroup = nRows == 1;

% Decide on the color map
%   Note: Only do this when colorMap is nonempty so that
%           the default of bar() will be applied otherwise
if ~isempty(colorMap)
    colorMap = decide_on_colormap(colorMap, nCols);
end

% Set the default confidence interval line color
if isempty(cIColor)
    if singleGroup
        % All bars are the same color, so use red
        cIColor = 'r';
    else
        % All bars are different colors, so use black
        cIColor = 'k';
    end
end

% Set the default parameter values
if isempty(pValues)
    if oneSamplePerGroup
        pValues = 1:nCols;
    else
        pValues = transpose(1:nRows);
    end
end

% Set the default parameter tick angle
if isempty(pTickAngle)
    if singleGroup && ~isempty(pTickLabels)
        % One group only, so make the tick labels slanted
        pTickAngle = 75;
    else
        % Multiple groups, so no need to slant tick labels
        pTickAngle = 0;
    end
end

% Set bar direction-dependent parameters
switch barDirection
    case 'vertical'
        xLimits = pLimits;
        yLimits = readoutLimits;
        xLabel = pLabel;
        yLabel = readoutLabel;
    case 'horizontal'
        xLimits = readoutLimits;
        yLimits = pLimits;
        xLabel = readoutLabel;
        yLabel = pLabel;
end

% Set column (group) labels
if isempty(columnLabels)
    columnLabels = cell(1, nCols);
    for iGroup = 1:nCols
        columnLabels{iGroup} = ['Group #', num2str(iGroup)];
    end
end

% Set legend location based on number of groups
if strcmpi(legendLocation, 'auto')
    if oneSamplePerGroup || nCols == 1
        legendLocation = 'suppress';
    elseif nCols > 1 && nCols < 10
        legendLocation = 'northeast';
    else
        legendLocation = 'eastoutside';
    end
end

% Set the default figure title
if isempty(figTitle)
    if ~strcmpi(readoutLabel, 'suppress') && ~strcmpi(pLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. ', pLabel], '_', '\_');
    elseif ~strcmpi(readoutLabel, 'suppress')
        figTitle = strrep([readoutLabel, ' vs. parameter'], '_', '\_');
    else
        figTitle = 'Readout vs. parameter';
    end
end

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
        {pValues, val, phaseVectors, ...
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

%% Plot bars
% Decide on the figure to plot on
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber);

% Reverse the order of pValues on the plot if requested
if reverseOrder
    pValuesToPlot = flip(pValues);
else
    pValuesToPlot = pValues;
end

% Draw bar graph
switch barDirection
    case 'vertical'
        if oneSamplePerGroup && strcmp(groupStyle, 'grouped')
            % One sample per group, but the bar() function
            %   will automatically construe it as one group. 
            %   Therefore, plot a stacked grouped bar graph to do the trick
            bars = bar(pValuesToPlot, diag(val), barWidth, 'stacked', ...
                        otherArguments{:});
        else
            % Just use the bar() function
            bars = bar(pValuesToPlot, val, barWidth, groupStyle, ...
                        otherArguments{:});
        end
    case 'horizontal'
        if oneSamplePerGroup && strcmp(groupStyle, 'grouped')
            bars = barh(pValuesToPlot, diag(val), barWidth, 'stacked', ...
                        otherArguments{:});
        else
            % Just use the barh() function
            bars = barh(pValuesToPlot, val, barWidth, groupStyle, ...
                        otherArguments{:});
        end
    otherwise
        error('barDirection unrecognized!');
end

% Change pTicks if provided
if ~isempty(pTicks)
    switch barDirection
        case 'vertical'
            xticks(pTicks);
        case 'horizontal'
            yticks(pTicks);
    end
end

% Change pTickLabels if provided
if ~isempty(pTickLabels)
    % Get current tick positions
    switch barDirection
        case 'vertical'
            pTicksNow = get(gca, 'XTick');
        case 'horizontal'
            pTicksNow = get(gca, 'YTick');
    end

    % Extract the corresponding tick labels
    if numel(pTickLabels) == numel(pValues)
        pTickLabelsNeeded = match_positions(pTickLabels, pValues, pTicksNow);
    else
        pTickLabelsNeeded = pTickLabels;
    end

    % Plot the tick locations
    switch barDirection
        case 'vertical'
            xticklabels(pTickLabelsNeeded);
        case 'horizontal'
            yticklabels(pTickLabelsNeeded);
    end
end

% Set display name
if ~strcmpi(columnLabels, 'suppress')
    for iGroup = 1:nCols
        set(bars(iGroup), 'DisplayName', columnLabels{iGroup});
    end
end

% Reverse pTicks and pTickLabels if requested
if reverseOrder
    switch barDirection
        case 'vertical'
            % Get current ticks and labels
            xTicksOld = get(gca, 'XTick');
            xTickLabelsOld = get(gca, 'XTickLabel');

            % Flip the ticks and labels
            xTicksNew = flip(pValues(1) + pValues(end) - xTicksOld);
            xTickLabelsNew = flip(xTickLabelsOld);

            % Set new ticks and labels
            xticks(xTicksNew);
            xticklabels(xTickLabelsNew);
        case 'horizontal'
            % Get current ticks and labels
            yTicksOld = get(gca, 'YTick');
            yTickLabelsOld = get(gca, 'YTickLabel');

            % Flip the ticks and labels
            yTicksNew = flip(pValues(1) + pValues(end) - yTicksOld);
            yTickLabelsNew = flip(yTickLabelsOld);

            % Set new ticks and labels
            yticks(yTicksNew);
            yticklabels(yTickLabelsNew);
    end
end

% Change the parameter tick angle
switch barDirection
    case 'vertical'
        xtickangle(pTickAngle);
    case 'horizontal'
        ytickangle(pTickAngle);
end

%% Plot error bars
% Set the default confidence interval bar width
if isempty(cIBarWidth)
    if singleGroup || oneSamplePerGroup
        barSeparation = 1;
    else
        % Extract the parameter offsets
        %   Note: this is XOffset for BOTH bar() and barh()
        offsets = arrayfun(@(x) bars(x).XOffset, 1:nCols);
        
        % Compute the bar separation
        barSeparation = mean(diff(offsets));
    end

    % Compute the confidence interval bar width
    cIBarWidth = cIRelativeBarWidth * barSeparation;
end

% Set the color for each Bar object
% TODO: Make edgeColor an optional argument
% TODO: Edge color can't seem to be set different for each bar?
if ~isempty(colorMap)
    for iBar = 1:numel(bars)
        set(bars(iBar), 'FaceColor', colorMap(iBar, :), ...
            'EdgeColor', 'none');
    end
end

% Plot error bars
if ~isempty(low) || ~isempty(high)
    hold on;
    if singleGroup || oneSamplePerGroup
        % Draw error bars
        lines = plot_error_bar(pValuesToPlot, low, high, ...
                                'BarDirection', barDirection, ...
                                'BarWidth', cIBarWidth, ...
                                'Color', cIColor, 'LineWidth', cILineWidth);
    else                % Data is grouped
        lines = gobjects(3, nRows, nCols);
        for iCol = 1:nCols              % for each group
            % Get the parameter offset
            %   Note: this is XOffset for BOTH bar() and barh()
            pOffset = bars(iCol).XOffset;

            for iRow = 1:nRows                  % for each sample
                % Compute parameter position
                pPos = pValuesToPlot(iRow) + pOffset;

                % Draw error bar
                lines(:, iRow, iCol) = ...
                    plot_error_bar(pPos, low(iRow, iCol), high(iRow, iCol), ...
                                    'BarDirection', barDirection, ...
                                    'BarWidth', cIBarWidth, ...
                                    'Color', cIColor, 'LineWidth', cILineWidth);
            end
        end
    end
else
    lines = gobjects(1);
end

% Set axes limits
if ~isempty(xLimits) && ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end
if ~isempty(yLimits) && ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Set axes labels and title
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

% Plot parameter boundaries
if nPBoundaries > 0
    % Compute flipped boundaries if requested
    if reverseOrder
        pBoundariesToPlot = pValues(1) + pValues(end) - pBoundaries;
    else
        pBoundariesToPlot = pBoundaries;
    end

    hold on
    pLines = plot_horizontal_line(pBoundariesToPlot, 'LineWidth', 0.5, ...
                                        'LineStyle', '--', 'Color', 'g');
else
    pLines = gobjects;
end

% Plot readout boundaries
if nRBoundaries > 0
    hold on
    rLines = plot_vertical_line(rBoundariesToPlot, 'LineWidth', 0.5, ...
                                        'LineStyle', '--', 'Color', 'r');
else
    rLines = gobjects;
end

% Plot phase averages if any
if plotPhaseAverages && ~isempty(phaseAverages) && ~isempty(averageWindows)
    % Compute flipped indices if requested
    if reverseOrder
        averageWindowsToPlot = cellfun(@(x) pValues(1) + pValues(end) - x, ...
                                        averageWindows, 'UniformOutput', false);    
    else
        averageWindowsToPlot = averageWindows;
    end

    % Plot the phase averages as horizontal lines
    averages = ...
        arrayfun(@(x) plot_vertical_line(phaseAverages(x), ...
                                'YLimits', averageWindowsToPlot{x}, ...
                                'LineStyle', averagesLineStyle, ...
                                'LineWidth', averagesLineWidth), ...
            transpose(1:size(phaseAverages, 1)), 'UniformOutput', false);
end

% Plot selected values if any
if plotIndSelected && ~isempty(indSelected)
    hold on
    % Plot selected values
    if iscell(indSelected)
        % TODO: Make cellfun_custom.m
        selectedCell = cellfun(@(y) plot_selected(val, pValuesToPlot, y), ...
                            indSelected, 'UniformOutput', false);
        selected = vertcat(selectedCell{:});
    else
        selected = plot_selected(val, pValuesToPlot, indSelected);
    end
else
    selected = gobjects;
end

%% Post-plotting
% Generate a legend if requested
if ~strcmpi(legendLocation, 'suppress')
    lgd = legend(bars, 'location', legendLocation);

    set(lgd, 'AutoUpdate', 'off', 'Interpreter', 'none');
end

% Save in output
handles.fig = fig;
handles.bars = bars;
if exist('lines', 'var')
    handles.lines = lines;
end
if exist('pLines', 'var') && exist('rLines', 'var')
    handles.boundaries = transpose(vertcat(pLines, rLines));
end
if exist('averages', 'var')
    handles.averages = averages;
end
if exist('selected', 'var')
    handles.selected = selected;
end

% Save figure if figName provided
if ~isempty(figName)
    save_all_figtypes(fig, figName, figTypes);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

figure(h);

    % TODO: Make each bar a different color 
barColor = 'blue';
    %% bars.FaceColor = rgb(barColor);

%% bars.CData = colormap(lines(nCols));    % TODO: Not working!

% Set the default confidence interval bar color
if isempty(cIColor)
    if nRows == 1
        cIColor = 'r';
    else
        cIColor = 'k';
    end
end

if nCols == 1
    % One sample per group, so use group numbers
else
    % Many samples per group, so use sample numbers
    pValues = 1:nCols;
end

for iCol = 1:nCols                  % for each sample
    pPos = iCol; 
    lines(1, iCol) = line(pPos * ones(1, 2), ...
                    [low(iCol), high(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
    lines(2, iCol) = line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [low(iCol), low(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
    lines(3, iCol) = line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [high(iCol), high(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
end

pPos = iRow + (iCol - (nCols + 1) / 2) * barSeparation; 
lines(1, iCol, iRow) = ...
    line(pPos * ones(1, 2), ...
        [low(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(2, iCol, iRow) = ...
    line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [low(iRow, iCol), low(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(3, iCol, iRow) = ...
    line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [high(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);

[bars, lines, fig] = plot_bar(val1, low1, high1);

% Decide whether there is only one group before any change
singleGroup = nRows == 1;
% Count the number of rows (groups)
nRows = size(val, 1);

lines(1, iCol, iRow) = ...
    line(pPos * ones(1, 2), ...
        [low(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(2, iCol, iRow) = ...
    line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [low(iRow, iCol), low(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(3, iCol, iRow) = ...
    line(pPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [high(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);

pPos = iRow + (iCol - (nCols + 1) / 2) * barSeparation; 

%                   each row is a different group
%                   each column is a different sample number

singleGroup = nRows == 1;
oneSamplePerGroup = nCols == 1;

%                   - 'BarSeparation': TODO
barSeparationDefault = [];
addParameter(iP, 'BarSeparation', barSeparationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
barSeparation = iP.Results.BarSeparation;
% Set the default bar separation
if isempty(barSeparation)
    if singleGroup || oneSamplePerGroup
        barSeparation = 1;
    else
        barSeparation = 1/(nCols + 2);
    end
end

% Set figure as current figure
if isempty(figHandle)
    set(0, 'CurrentFigure', figHandle);
else
    figure(figHandle);
end
fig = gcf;

if oneSamplePerGroup
    pValues = fliplr(pValues);
else
    pValues = flipud(pValues);
end

% Reverse the order of pValues, pTicks and pTick labels
[pValues, pTickLabels] = ...
    argfun(@flip, pValues, pTickLabels);

set(gca, 'XTick', pTicks);
set(gca, 'YTick', pTicks);
set(gca, 'XTickLabel', pTickLabels);
set(gca, 'YTickLabel', pTickLabels);

% Compute the minimum and maximum pValue
minPValue = min(pValues);
maxPValue = max(pValues);

nBoundaries = nPBoundaries + nRBoundaries;

% Get the selected values
valSelected = val(indSelected);

% Plot red crosses
hold on
plot(valSelected, indSelected, 'rx', 'LineWidth', 2, 'MarkerSize', 6);

~all(isnan(indSelected))

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
