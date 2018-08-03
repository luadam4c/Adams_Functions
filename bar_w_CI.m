function h = bar_w_CI(h, data, low, high, varargin)
%% Plot bar graph (esp. grouped) with confidence intervals
% Usage: h = bar_w_CI(h, data, low, high, varargin)
% Example:
%       For a data mu with 2 rows and 7 columns, 
%       h = bar_w_CI(h, mu, muLow, muHigh, ...
%                       'BarSeparation', 0.115, ...
%                       'CIBarWidth', 0.05, ...
%                       'CILineWidth', 2, ...
%                       'CIColor', 'k', ...
%                       'XValues', 1:2);
%                       'XTickLabel', {'Intra-burst', 'Inter-burst'});
% Arguments: TODO
%
% Requires:
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:
%       /home/Matlab/Adams_Functions/ZG_fit_IEI_distributions.m
%       /media/adamX/Paula_IEIs/paula_iei3.m
%       /media/adamX/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m
%
% File History: 
% 2017-10-19 - Moved from paula_iei3.m
% 2017-10-19 - Added input parser and various optional arguments
% 2017-12-01 - Added figure(h)
% 2019-03-08 - Changed figure(h) to set(0, 'CurrentFigure', h)
%               to prevent the display of invisible figures
% 

%% Parameters
barColor = 'blue';

%% Default values for optional arguments
cILineWidthDefault = 2;             % default line width for CIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(['Not enough input arguments, ', ...
            'type ''help bar_w_CI'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'bar_w_CI';

% Add required inputs to the Input Parser
addRequired(iP, 'h')                        % figure handle
addRequired(iP, 'data', ...                 % data
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

% Read from the Input Parser
parse(iP, h, data, low, high, varargin{:});
barSeparation = iP.Results.BarSeparation;
cIBarWidth = iP.Results.CIBarWidth;
cILineWidth = iP.Results.CILineWidth;
cIColor = iP.Results.CIColor;
xValues = iP.Results.XValues;
xTickLabel = iP.Results.XTickLabel;
xTickAngle = iP.Results.XTickAngle;

% Set dependent argument defaults
nRows = size(data, 1);
nCols = size(data, 2);
if isempty(barSeparation)
    barSeparation = 1/(nCols+2);
end
if isempty(cIBarWidth)
    if nRows == 1
        cIBarWidth = 0.5;
    else
        cIBarWidth = 0.05;
    end
end
if isempty(cIColor)
    if nRows == 1
        cIColor = 'r';
    else
        cIColor = 'k';
    end
end
if isempty(xTickAngle)
    if nRows == 1
        xTickAngle = 75;
    else
        xTickAngle = 0;
    end
end

%% Plot things
% Set figure as current figure
set(0, 'CurrentFigure', h);

% Draw bar graph
if isempty(xValues)
    b = bar(data);
else
    b = bar(xValues, data);
end

% Change xTickLabel if provided
if ~isempty(xTickLabel)
    set(gca, 'XTickLabel', xTickLabel);
end

% Change the X Tick Angle
% xtickangle(xTickAngle);

% Plot error bars
hold on;
if nRows == 1       % Data is not grouped
    % Make each bar a different color 
    %% b.CData = colormap(lines(nCols));    % TODO: Not working!
    b.FaceColor = rgb(barColor);
       
    for iCol = 1:nCols              % for each column
        % Draw error bar
        xPos = iCol; 
        line(xPos * ones(1, 2), ...
                [low(iCol), high(iCol)], ...
                'Color', cIColor, 'LineWidth', cILineWidth);
        line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                [low(iCol), low(iCol)], ...
                'Color', cIColor, 'LineWidth', cILineWidth);
        line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                [high(iCol), high(iCol)], ...
                'Color', cIColor, 'LineWidth', cILineWidth);
    end
else                % Data is grouped
    for iRow = 1:nRows                  % for each row
        for iCol = 1:nCols              % for each column
            % Draw error bar
            xPos = iRow + (iCol - (nCols+1)/2) * barSeparation; 
            line(xPos * ones(1, 2), ...
                    [low(iRow, iCol), high(iRow, iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
            line(xPos * ones(1, 2) + [-cIBarWidth/2, cIBarWidth/2], ...
                    [low(iRow, iCol), low(iRow, iCol)], ...
                    'Color', cIColor, 'LineWidth', cILineWidth);
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
