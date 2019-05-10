function [bars, lines, fig] = plot_bar(val, varargin)
%% Plots a bar graph (grouped or not) with or without confidence intervals
% Usage: [bars, lines, fig] = plot_bar(val, varargin)
% Explanation:
%       TODO
% Example:
%       val1 = [2, 5, 3, 13]; low = val1 - 1; high = val1 + 1;
%       val2 = [2 2 3; 2 5 6; 2 8 9; 2 11 12]; low = val2 - 1; high = val2 + 1;
%       [bars, lines, fig] = plot_bar(val1, low, high);
%       [bars, lines, fig] = plot_bar(val2, low, high, 'PTickLabels', {'Mark', 'Ashley', 'Katie', 'Adam'});
%       [bars, lines, fig] = plot_bar(val2, low, high, 'BarDirection', 'horizontal');
% Arguments:
%       val     - mean values for the bar() function
%                   each row is a different group
%                   each column is a different sample number
%               must be a numeric array accepted by the bar() function
%       low     - (opt) lower limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       high    - (opt) upper limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       varargin    - 'BarDirection': bar direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'vertical'   - vertical bars
%                       'horizontal' - horizontal bars
%                   default == 'vertical'
%                   - 'BarSeparation': TODO
%                   - 'CIBarWidth': TODO
%                   - 'CILineWidth': TODO
%                   - 'CIColor': TODO
%                   - 'PValues': TODO
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
%                   - 'PTickAngle': TODO
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == 'Readout'
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - Any other parameter-value pair for the bar() function
%   
%
% Requires:
%       cd/argfun.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_struct.m
%       cd/ZG_fit_IEI_distributions.m
%       /media/adamX/Paula_IEIs/paula_iei3.m
%
% File History: 
% 2017-10-19 - Moved from paula_iei3.m
% 2017-10-19 - Added input parser and various optional arguments
% 2017-12-01 - Added figure(h)
% 2018-03-08 - Changed figure(h) to set(0, 'CurrentFigure', h)
%               to prevent the display of invisible figures
% 2019-01-15 - Renamed bar_w_CI.m -> plot_bar.m
% 2019-01-15 - Added otherArguments
% 2019-01-15 - Made h -> 'FigHandle' an optional argument
% 2019-05-08 - Added 'barDirection' as an optional argument
% 2019-05-08 - Made low and high optional arguments
% 2019-05-10 - Fixed bugs
% TODO: Fix error bars when 'barDirection' is 'horizontal'
% TODO: Update the code to use plot_horizontal_line.m and plot_vertical_line.m
% TODO: If 'TreatVectorAsArray' is true, don't force all vectors as column vectors
% TODO: Add 'BarColors' as an optional argument
% TODO: Change usage in all functions using this
% 

%% Hard-coded parameters
validBarDirections = {'vertical', 'horizontal'};

% TODO: Make this an optional parameter
treatVectorAsArray = false;

%% Default values for optional arguments
lowDefault = [];
highDefault = [];
barDirectionDefault = 'vertical';
barSeparationDefault = [];
cIBarWidthDefault = [];
cILineWidthDefault = 2;                 % default line width for CIs
cIColorDefault = '';
pValuesDefault = [];
pLimitsDefault = [];
readoutLimitsDefault = [];
pTicksDefault = [];
pTickLabelsDefault = {};
pTickAngleDefault = [];
pLabelDefault = 'Parameter';
readoutLabelDefault = 'Readout';
figHandleDefault = [];                  % no existing figure by default

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
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'low', lowDefault, ...          % low limit of CI
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addOptional(iP, 'high', highDefault, ...        % high limit of CI
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarDirection', barDirectionDefault, ...
    @(x) any(validatestring(x, validBarDirections)));
addParameter(iP, 'BarSeparation', barSeparationDefault, ...
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
addParameter(iP, 'FigHandle', figHandleDefault);

% Read from the Input Parser
parse(iP, val, varargin{:});
low = iP.Results.low;
high = iP.Results.high;
barDirection = validatestring(iP.Results.BarDirection, validBarDirections);
barSeparation = iP.Results.BarSeparation;
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
figHandle = iP.Results.FigHandle;

% Keep unmatched arguments for the bar() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
        numel(pTicks) ~= numel(pTickLabels)
    fprintf('PTicks and PTickLabels must have the same number of elements!\n');
    bars = gobjects;
    lines = gobjects;
    fig = gobjects;
    return
end

%% Preparation
% Count the number of rows (groups)
nRows = size(val, 1);

% Decide whether there is only one group before any change
singleGroup = nRows == 1;

% Force row vectors as column vectors
%   Note: This will cause each value of a vector to be plotted as separately
%           colored bars
if ~treatVectorAsArray
    [val, low, high] = ...
        argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
                val, low, high);
end

% Count the number of columns (samples)
nCols = size(val, 2);

% Count the number of rows (groups)
nRows = size(val, 1);

% Decide whether there is only one row
singleRow = nRows == 1;

% Decide whether there is one sample per group
oneSamplePerGroup = nCols == 1;

% Set the default bar separation
if isempty(barSeparation)
    barSeparation = 1/(nCols + 2);
end

% Set the default confidence interval bar width
if isempty(cIBarWidth)
    if singleGroup
        cIBarWidth = 0.5;
    else
        cIBarWidth = 0.05;
    end
end

% Set the default confidence interval line color
if isempty(cIColor)
    if singleRow
        % All bars are the same color, so use red
        cIColor = 'r';
    else
        % All bars are different colors, so use black
        cIColor = 'k';
    end
end

% Set the default x values
if isempty(pValues)
    pValues = 1:nRows;
end

% Set the default x tick angle
if isempty(pTickAngle)
    if singleGroup
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

%% Plot things
% Set figure as current figure
if isempty(figHandle)
    set(0, 'CurrentFigure', figHandle);
else
    figure(figHandle);
end
fig = gcf;

% Draw bar graph
switch barDirection
    case 'vertical'
        if oneSamplePerGroup
            % One sample per group, but the bar() function
            %   will automatically construe it as one group. 
            %   Therefore, plot a stacked grouped bar graph to do the trick
            bars = bar(pValues, diag(val), 'stacked', otherArguments{:});
        else
            % Just use the bar() function
            bars = bar(pValues, val, otherArguments{:});
        end
    case 'horizontal'
        if oneSamplePerGroup
            bars = barh(pValues, diag(val), 'stacked', otherArguments{:});
        else
            % Just use the barh() function
            bars = barh(pValues, val, otherArguments{:});
        end
    otherwise
        error('barDirection unrecognized!');
end

% Set the color for each Bar object
% TODO
% for iBar = 1:numel(bars)
%     set(bars(iBar), 'CData', barColors{iBar});
% end

% Change pTickLabels if provided
if ~isempty(pTicks)
    switch barDirection
        case 'vertical'
            set(gca, 'XTick', pTicks);
            % xticks(pTicks);
        case 'horizontal'
            set(gca, 'YTick', pTicks);
            % yticks(pTicks);
    end
end

% Change pTickLabels if provided
if ~isempty(pTickLabels)
    switch barDirection
        case 'vertical'
            set(gca, 'XTickLabel', pTickLabels);
            % xticklabels(pTickLabels);
        case 'horizontal'
            set(gca, 'YTickLabel', pTickLabels);
            % yticklabels(pTickLabels);
    end
end

% Change the X Tick Angle
xtickangle(pTickAngle);

% Plot error bars
if ~isempty(low) || ~isempty(high)
    hold on;
    if singleGroup
            % Draw error bar
            lines = plot_error_bar(pValues, low, high, 'BarWidth', cIBarWidth);
    else                % Data is grouped
    %{
        for iRow = 1:nRows                  % for each group
            for iCol = 1:nCols              % for each sample
                % Draw error bar
                xPos = iRow + (iCol - (nCols + 1) / 2) * barSeparation; 
                lines(1, iCol, iRow) = ...
                    line(xPos * ones(1, 2), ...
                        [low(iRow, iCol), high(iRow, iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
                lines(2, iCol, iRow) = ...
                    line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                        [low(iRow, iCol), low(iRow, iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
                lines(3, iCol, iRow) = ...
                    line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                        [high(iRow, iCol), high(iRow, iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
            end
        end
    %}
        lines = gobjects(1);
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

% Set axes labels
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
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
    xPos = iCol; 
    lines(1, iCol) = line(xPos * ones(1, 2), ...
                    [low(iCol), high(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
    lines(2, iCol) = line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [low(iCol), low(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
    lines(3, iCol) = line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [high(iCol), high(iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
end

xPos = iRow + (iCol - (nCols + 1) / 2) * barSeparation; 
lines(1, iCol, iRow) = ...
    line(xPos * ones(1, 2), ...
        [low(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(2, iCol, iRow) = ...
    line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [low(iRow, iCol), low(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);
lines(3, iCol, iRow) = ...
    line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
        [high(iRow, iCol), high(iRow, iCol)], ...
        'Color', cIColor, 'LineWidth', cILineWidth);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
