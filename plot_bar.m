function [bars, lines, fig] = plot_bar(val, low, high, varargin)
%% Plots a bar graph (grouped or not) with confidence intervals
% Usage: [bars, lines, fig] = plot_bar(val, low, high, varargin)
% Explanation:
%       TODO
% Example:
%       val = [2 2 3; 2 5 6; 2 8 9; 2 11 12]; low = val - 1; high = val + 1;
%       [bars, lines, fig] = plot_bar(val, low, high);
%       [bars, lines, fig] = plot_bar(val, low, high, 'XTickLabel', {'Mark', 'Ashley', 'Katie', 'Adam'});
% Arguments:
%       val   - mean values for the bar() function
%                   each row is a different group
%                   each column is a different sample number
%               must be a numeric array accepted by the bar() function
%       low     - lower limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       high    - upper limits of the confidence intervals
%               must be a vector with numel same as the number of columns in val
%       varargin    - 'BarSeparation': TODO
%                   - 'CIBarWidth': TODO
%                   - 'CILineWidth': TODO
%                   - 'CIColor': TODO
%                   - 'XValues': TODO
%                   - 'XTickLabel': TODO
%                   - 'XTickAngle': TODO
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
%       TODO: cd/m3ha_neuron_run_and_analyze.m
%       TODO: cd/ZG_fit_IEI_distributions.m
%       TODO: /media/adamX/Paula_IEIs/paula_iei3.m
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
% TODO: Update the code to use plot_horizontal_line.m and plot_vertical_line.m
% TODO: If 'TreatVectorAsArray' is true, don't force all vectors as column vectors
% TODO: Add 'BarColors' as an optional argument
% TODO: Change usage in all functions using this
% 

%% Hard-coded parameters
% TODO: Make this an optional parameter
treatVectorAsArray = false;

%% Default values for optional arguments
cILineWidthDefault = 2;                 % default line width for CIs
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
addRequired(iP, 'val', ...                  % values
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addRequired(iP, 'low', ...                  % low limit of CI
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addRequired(iP, 'high', ...                 % high limit of CI
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarSeparation', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CIBarWidth', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CILineWidth', cILineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CIColor', '');
addParameter(iP, 'XValues', [], ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addParameter(iP, 'XTickLabel', '', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'XTickAngle', '', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'FigHandle', figHandleDefault);

% Read from the Input Parser
parse(iP, val, low, high, varargin{:});
barSeparation = iP.Results.BarSeparation;
cIBarWidth = iP.Results.CIBarWidth;
cILineWidth = iP.Results.CILineWidth;
cIColor = iP.Results.CIColor;
xValues = iP.Results.XValues;
xTickLabel = iP.Results.XTickLabel;
xTickAngle = iP.Results.XTickAngle;
figHandle = iP.Results.FigHandle;

% Keep unmatched arguments for the bar() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force column vectors as column vectors
%   Note: This will cause each value of a vector to be plotted as separately
%           colored bars
if ~treatVectorAsArray
    [val, low, high] = ...
        argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
                val, low, high);
end

% Count the number of rows (groups)
nRows = size(val, 1);

% Count the number of columns (samples)
nCols = size(val, 2);

% Decide whether there is only one group
singleGroup = nRows == 1;

% Decide whether there is one sample per group
oneSamplePerGroup = nCols == 1;

% Set the default bar separation
if isempty(barSeparation)
    barSeparation = 1/(nCols + 2);
end

% Set the default confidence interval bar width
if isempty(cIBarWidth)
    if nRows == 1
        cIBarWidth = 0.5;
    else
        cIBarWidth = 0.05;
    end
end

% Set the default confidence interval line color
if isempty(cIColor)
    if nRows == 1
        % All bars are the same color, so use red
        cIColor = 'r';
    else
        % All bars are different colors, so use black
        cIColor = 'k';
    end
end

% Set the default x values
if isempty(xValues)
    if nCols == 1
        % One sample per group, so use group numbers
        xValues = 1:nRows;
    else
        % Many samples per group, so use sample numbers
        xValues = 1:nCols;
    end
end

% Set the default x tick angle
if isempty(xTickAngle)
    if singleGroup
        % One group only, so make the tick labels slanted
        xTickAngle = 75;
    else
        % Multiple groups, so no need to slant tick labels
        xTickAngle = 0;
    end
end

%% Plot things
% Set figure as current figure
if isempty(figHandle)
    set(0, 'CurrentFigure', figHandle);
end
fig = gcf;

% Draw bar graph
if oneSamplePerGroup
    % One sample per group, but the bar() function
    %   will automatically construe it as one group. 
    %   Therefore, plot a stacked grouped bar graph to do the trick
    bars = bar(xValues, diag(val), 'stacked', otherArguments{:});
else
    % Just use the bar() function
    bars = bar(xValues, val, otherArguments{:});
end

% Set the color for each Bar object
% TODO
% for iBar = 1:numel(bars)
%     set(bars(iBar), 'CData', barColors{iBar});
% end

% Change xTickLabel if provided
if ~isempty(xTickLabel)
    set(gca, 'XTickLabel', xTickLabel);
end

% Change the X Tick Angle
xtickangle(xTickAngle);

% Plot error bars
hold on;
if singleGroup
    for iCol = 1:nCols                  % for each sample
        % Draw error bar
        % TODO: function plot_error_bar(x, yLow, yHigh, 'BarWidth', barWidth)
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
else                % Data is grouped
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
