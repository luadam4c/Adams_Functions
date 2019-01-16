function [bars, lines, fig] = plot_bar(means, low, high, varargin)
%% Plot bar graph (esp. grouped) with confidence intervals
% Usage: [bars, lines, fig] = plot_bar(means, low, high, varargin)
% Explanation:
%       TODO
% Example:
%       means = [2 2 3; 2 5 6; 2 8 9; 2 11 12];
%       low   = means - 1;
%       high  = means + 1;
%       For a means mu with 2 rows and 7 columns, 
%       [bars, lines] = plot_bar(means, low, high, ...
%                               'BarSeparation', 0.115, ...
%                               'CIBarWidth', 0.05, ...
%                               'CILineWidth', 2, ...
%                               'CIColor', 'k', ...
%                               'XTickLabel', {'Mark', 'Peter', 'Katie'});
% Arguments: TODO
%       means   - mean values for the bar() function
%                   each row is a different group
%                   each column is a different sample number
%               must be a numeric array accepted by the bar() function
%       low     - lower limits of the confidence intervals
%               must be a vector with numel same as the number of columns in means
%       high    - upper limits of the confidence intervals
%               must be a vector with numel same as the number of columns in means
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
% TODO: Change usage in all functions using this
% 

%% Parameters
barColor = 'blue';

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
addRequired(iP, 'h')                        % figure handle
addRequired(iP, 'means', ...                 % means
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
addParameter(iP, 'CIColor', '');    % TODO: validation
addParameter(iP, 'XValues', [], ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addParameter(iP, 'XTickLabel', '', ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'XTickAngle', '', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'FigHandle', figHandleDefault);

% Read from the Input Parser
parse(iP, h, means, low, high, varargin{:});
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
% Count the number of rows (groups)
nRows = size(means, 1);

% Count the number of columns (samples)
nCols = size(means, 2);

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

% Set the default confidence interval bar color
if isempty(cIColor)
    if nRows == 1
        cIColor = 'r';
    else
        cIColor = 'k';
    end
end

% Set the default x tick angle
if isempty(xTickAngle)
    if nRows == 1
        xTickAngle = 75;
    else
        xTickAngle = 0;
    end
end

%% Plot things
% Set figure as current figure
if isempty(figHandle)
    set(0, 'CurrentFigure', figHandle);
end
fig = gcf

% Draw bar graph
if isempty(xValues)
    bars = bar(means, otherArguments{:});
else
    bars = bar(xValues, means, otherArguments{:});
end

% Change xTickLabel if provided
if ~isempty(xTickLabel)
    set(gca, 'XTickLabel', xTickLabel);
end

% Change the X Tick Angle
xtickangle(xTickAngle);

% Plot error bars
hold on;
if nRows == 1       % Data is not grouped
    % Make each bar a different color 
    %% b.CData = colormap(lines(nCols));    % TODO: Not working!
    b.FaceColor = rgb(barColor);
       
    for iCol = 1:nCols              % for each column
        % Draw error bar
        xPos = iCol; 
        lines(iCol, 1) = line(xPos * ones(1, 2), ...
                        [low(iCol), high(iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
        lines(iCol, 2) = line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                        [low(iCol), low(iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
        lines(iCol, 3) = line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                        [high(iCol), high(iCol)], ...
                        'Color', cIColor, 'LineWidth', cILineWidth);
    end
else                % Data is grouped
    for iRow = 1:nRows                  % for each row
        for iCol = 1:nCols              % for each column
            % Draw error bar
            xPos = iRow + (iCol - (nCols+1)/2) * barSeparation; 
            lines(iRow, iCol, 1) = ...
                line(xPos * ones(1, 2), ...
                    [low(iRow, iCol), high(iRow, iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
            lines(iRow, iCol, 2) = ...
                line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [low(iRow, iCol), low(iRow, iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
            lines(iRow, iCol, 3) = ...
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
