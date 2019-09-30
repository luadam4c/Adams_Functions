function h = plot_window_boundaries (win, varargin)
%% Plots window boundaries as separating lines, duration bars or background shades
% Usage: h = plot_window_boundaries (win, varargin)
% Explanation:
%       Plots a vector of window boundaries in potentially many different ways
%
% Example(s):
%       figure(2); clf; load_examples; plot(myTimeVec, myRandomSignals1)
%       plot_window_boundaries([5 10 20 25])
%       plot_window_boundaries([1 2 3], 'BoundaryType', 'horizontalLines')
%       plot_window_boundaries([1.5 2 3 3.5], 'BoundaryType', 'verticalBars')
%       plot_window_boundaries([5 10 20 25], 'BoundaryType', 'horizontalBars')
%       plot_window_boundaries([5 10 20 25], 'BoundaryType', 'verticalShades')
%       plot_window_boundaries([2 3], 'BoundaryType', 'horizontalShades')
%
% Outputs:
%       h           - handles to each line object (left, right)
%                   specified as a 2-element column array 
%                       of primitive line object handles
% Arguments:
%       win         - window(s) to plot boundaries for
%                   must be a numeric vector
%       varargin    - 'BoundaryType': type of boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'verticalLines'     - vertical dotted lines
%                       'horizontalLines'   - horizontal dotted lines
%                       'verticalBars'      - vertical bars
%                       'horizontalBars'    - horizontal bars
%                       'verticalShades'    - vertical shades
%                       'horizontalShades'  - horizontal shades
%                   default == 'verticalLines'
%                   - 'LineStyle': line style of boundaries
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - 'LineWidth': color of boundaries
%                   must be empty or a positive scalar
%                   default == 3 for bars and 2 for lines
%                   - 'BarRelValue': value for bars relative to current axis
%                   must be empty or a positive scalar
%                   default == barRelValue relative to current axis
%                   - 'BarValue': value for bars
%                   must be empty or a numeric scalar
%                   default == barRelValue relative to current axis
%                   - 'ColorMap': color map passed in
%                   must be empty or a string/character vector
%                       or an n-by-3 numeric array
%                   default == depends on BoundaryType
%                   - Any other parameter-value pair for the line() function
% 
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/islinestyle.m
%       cd/plot_horizontal_line.m
%       cd/plot_vertical_line.m
%       cd/plot_horizontal_shade.m
%       cd/plot_vertical_shade.m
%
% Used by:    
%       cd/m3ha_plot_individual_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/plot_tuning_curve.m

% File History:
% 2018-10-29 Created by Adam Lu
% 2018-12-19 Now passes extra arguments
% 2018-12-19 Now returns object handles
% 2019-08-27 Added 'BoundaryType' as an optional argument with values:
%               'verticalLines' (default), 'horizontalBars'
% 2018-08-30 Added 'verticalShades', 'horizontalShades'
% 

%% Hard-coded parameters
validBoundaryTypes = {'verticalLines', 'horizontalLines', ...
                        'horizontalBars', 'verticalBars', ...
                        'verticalShades', 'horizontalShades'};
lineLineStyle = '--';
lineLineWidth = 2;
barLineStyle = '-';
barLineWidth = 3;
shadeLineStyle = 'none';
shadeLineWidth = 0.5;

%% Default values for optional arguments
boundaryTypeDefault = 'verticalLines';
lineStyleDefault = '';      % set later
lineWidthDefault = '';      % set later
barRelValueDefault = 0.1;   % set later
barValueDefault = [];       % set later
colorMapDefault = '';       % set later

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
addRequired(iP, 'win', ...
    @(x) validateattributes(x, {'numeric', 'duration', ...
                                'logical', 'datetime'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BoundaryType', boundaryTypeDefault, ...
    @(x) any(validatestring(x, validBoundaryTypes)));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BarRelValue', barRelValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'BarValue', barValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'ColorMap', colorMapDefault);

% Read from the Input Parser
parse(iP, win, varargin{:});
boundaryType = validatestring(iP.Results.BoundaryType, validBoundaryTypes);
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
barRelValue = iP.Results.BarRelValue;
barValue = iP.Results.BarValue;
colorMap = iP.Results.ColorMap;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output
h = gobjects;

% Set default line style
if isempty(lineStyle)
    switch boundaryType
        case {'verticalLines', 'horizontalLines'}
            lineStyle = lineLineStyle;
        case {'verticalBars', 'horizontalBars'}
            lineStyle = barLineStyle;
        case {'verticalShades', 'horizontalShades'}
            lineStyle = shadeLineStyle;
        otherwise
            error('boundaryType unrecognized!');
    end
end

% Set default line width
if isempty(lineWidth)
    switch boundaryType
        case {'verticalLines', 'horizontalLines'}
            lineWidth = lineLineWidth;
        case {'verticalBars', 'horizontalBars'}
            lineWidth = barLineWidth;
        case {'verticalShades', 'horizontalShades'}
            lineWidth = shadeLineWidth;
        otherwise
            error('boundaryType unrecognized!');
    end
end

% Set default bar value
if isempty(barValue)
    switch boundaryType
        case 'verticalBars'
            % Get the current x axis limits
            xLimitsNow = get(gca, 'XLim');

            % Compute a default y value for the bar
            barValue = xLimitsNow(1) + barRelValue * range(xLimitsNow);
        case 'horizontalBars'
            % Get the current y axis limits
            yLimitsNow = get(gca, 'YLim');

            % Compute a default y value for the bar
            barValue = yLimitsNow(1) + barRelValue * range(yLimitsNow);
        case {'verticalLines', 'horizontalLines', ...
                'verticalShades', 'horizontalShades'}
            % Keep empty
        otherwise
            error('boundaryType unrecognized!');
    end    
end

% Set default color map
if isempty(colorMap)
    switch boundaryType
        case 'verticalLines' 
            colorMap = 'g';
        case 'horizontalLines' 
            colorMap = 'r';
        case 'verticalBars'
            colorMap = 'DarkRed';
        case 'horizontalBars'
            colorMap = 'DarkGreen';
        case 'verticalShades'
            colorMap = 'LightGreen';
        case 'horizontalShades'
            colorMap = 'LightSalmon';
        otherwise
            error('boundaryType unrecognized!');
    end
end

% Force as a column
win = force_column_vector(win);

% Make sure the number of elements is allowed
switch boundaryType
    case {'verticalBars', 'horizontalBars', ...
            'verticalShades', 'horizontalShades'}
        if mod(numel(win), 2) ~= 0
            fprintf(['The number of elements in the ', ...
                    'first argument must be even!\n']);
            return
        end
    otherwise
        % Do nothing
end

% Transform into a cell array
switch boundaryType
    case {'verticalShades', 'horizontalShades'}
        if numel(win) > 2
            % Reshape as 2 rows
            win = reshape(win, 2, []);

            % Place into a cell array
            win = force_column_cell(win);
        end
    otherwise
        % Do nothing
end

%% Do the job
% Plot lines
switch boundaryType
    case 'verticalLines' 
        h = plot_vertical_line(win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    case 'horizontalLines' 
        h = plot_horizontal_line(win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    case 'verticalBars'
        h = plot_vertical_line(barValue, 'YLimits', win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    case 'horizontalBars'
        h = plot_horizontal_line(barValue, 'XLimits', win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    case 'verticalShades'
        h = plot_vertical_shade(win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    case 'horizontalShades'
        h = plot_horizontal_shade(win, ...
                            'LineStyle', lineStyle, 'LineWidth', lineWidth, ...
                            'ColorMap', colorMap, otherArguments);
    otherwise
        error('boundaryType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
