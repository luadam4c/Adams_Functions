function [supAx, axHandles] = resize_subplots_for_labels (varargin)
%% Resizes subplots to make room for overarching labels and title
% Usage: [supAx, axHandles] = resize_subplots_for_labels (varargin)
% Explanation:
%       This function finds all subplot axes in the current figure, 
%       resizes them to create margins, and then adds an invisible axes 
%       on top to hold an overarching x-label, y-label, and/or title.
%
% Example(s):
%       figure; subplot(2,2,1); plot(1:10); subplot(2,2,2); plot(rand(10,1)); subplot(2,2,3); plot(sin(1:10)); subplot(2,2,4); plot(cos(1:10));
%       [supAx, axHandles] = resize_subplots_for_labels('XLabel', 'Time (s)', 'YLabel', 'Amplitude', 'FigTitle', 'Various Waveforms');
%
% Outputs:
%       supAx       - Handle for the invisible "super" axes holding the labels
%                   specified as a graphics object handle
%       axHandles   - Handles for the original subplot axes that were resized
%                   specified as a vector of graphics object handles
%
% Arguments:
%       varargin    - 'AxesHandles': handles for the subplot axes
%                   must be a vector of axes handles
%                   default == all axes in current figure (excluding legends)
%                   - 'XLabel': overarching x-axis label
%                   must be a string scalar or character vector
%                   default == ''
%                   - 'YLabel': overarching y-axis label
%                   must be a string scalar or character vector
%                   default == ''
%                   - 'FigTitle': overarching figure title
%                   must be a string scalar or character vector
%                   default == ''
%                   - 'XMargin': normalized margin at the bottom for x-label
%                   must be a numeric scalar
%                   default == 0.08
%                   - 'YMargin': normalized margin on the left for y-label
%                   must be a numeric scalar
%                   default == 0.05
%                   - 'TMargin': normalized margin at the top for title
%                   must be a numeric scalar
%                   default == 0.08
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_table_parallel.m
%       cd/plot_traces.m
%       cd/virt_plot_amplitude_correlation.m

% File History:
% 2025-09-13 Created by Gemini using prompt by Adam Lu
%               to extract code from plot_traces.m
% 2025-09-13 Fixed code to create invisible axes at the proper position

%% Hard-coded parameters

%% Default values for optional arguments
axesHandlesDefault = [];        % set later
xLabelDefault = '';
yLabelDefault = '';
figTitleDefault = '';
xMarginDefault = 0.08;
yMarginDefault = 0.05;
tMarginDefault = 0.08;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true; % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AxesHandles', axesHandlesDefault, @(x) isempty(x) || all(isgraphics(x, 'axes')));
addParameter(iP, 'XLabel', xLabelDefault, @(x) ischar(x) || isstring(x));
addParameter(iP, 'YLabel', yLabelDefault, @(x) ischar(x) || isstring(x));
addParameter(iP, 'FigTitle', figTitleDefault, @(x) ischar(x) || isstring(x));
addParameter(iP, 'XMargin', xMarginDefault, @isnumeric);
addParameter(iP, 'YMargin', yMarginDefault, @isnumeric);
addParameter(iP, 'TMargin', tMarginDefault, @isnumeric);

% Read from the Input Parser
parse(iP, varargin{:});
axHandles = iP.Results.AxesHandles;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
figTitle = iP.Results.FigTitle;
xMargin = iP.Results.XMargin;
yMargin = iP.Results.YMargin;
tMargin = iP.Results.TMargin;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% If axes handles are not provided, find them
if isempty(axHandles)
    fig = gcf;
    % Find all axes that are not legends or previously made super axes
    axHandles = findobj(fig, 'type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'super_axis');
end

% If there is only one subplot (or none), do nothing and return
if numel(axHandles) <= 1
    supAx = gobjects(1);
    return;
end

%% Do the job
% Determine whether overarching labels are needed
xLabelNeeded = ~isempty(xLabel);
yLabelNeeded = ~isempty(yLabel);
titleNeeded = ~isempty(figTitle);

% Manually resize all subplots to make room for overarching labels
if xLabelNeeded || yLabelNeeded || titleNeeded
    % Define how much normalized space to reserve for each label
    xMarginToUse = 0;
    if xLabelNeeded, xMarginToUse = xMargin; end
    yMarginToUse = 0;
    if yLabelNeeded, yMarginToUse = yMargin; end
    tMarginToUse = 0;
    if titleNeeded, tMarginToUse = tMargin; end

    % Define the new total area available for plots after margins
    newPlotArea = [yMarginToUse, xMarginToUse, 1 - yMarginToUse, 1 - (xMarginToUse + tMarginToUse)];

    % Get the original outer position of all subplots
    originalPositions = get(axHandles, 'OuterPosition'); % Get all original outer positions

    % Calculate the new outer positions for each subplot
    for i = 1:numel(axHandles)
        % Get the original outer position of the subplot
        pos = originalPositions{i}; % [left, bottom, width, height]

        % Calculate the new outer position by scaling and shifting the 
        % original outer position to fit inside the new plot area.
        newPos = [newPlotArea(1) + pos(1) * newPlotArea(3), ...  % New Left
                  newPlotArea(2) + pos(2) * newPlotArea(4), ...  % New Bottom
                  pos(3) * newPlotArea(3), ...                   % New Width
                  pos(4) * newPlotArea(4)];                      % New Height
        
        % Apply the new outer position to the subplot
        set(axHandles(i), 'OuterPosition', newPos);
    end

    % Create a new, invisible axes that covers the original total plotting area
    supAx = axes('Position', newPlotArea, 'Visible', 'off', ...
                 'Units', 'normalized', 'Tag', 'super_axis');
    set(supAx, 'XTick', [], 'YTick', []); % Remove x and y ticks

    % Add the required labels to this new axes and make them visible
    if xLabelNeeded
        xlabel(supAx, xLabel, 'Visible', 'on');
    end
    if yLabelNeeded
        ylabel(supAx, yLabel, 'Visible', 'on');
    end
    if titleNeeded
        title(supAx, figTitle, 'Visible', 'on');
    end
else
    % If no labels are needed, return an empty graphics object
    supAx = gobjects(1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[allLefts, allBottoms, allWidths, allHeights] = ...
    argfun(@(x) extract_elements(originalPositions, 'specific', 'Index', x), ...
            1, 2, 3, 4);

% Get the Rights and Tops
allRights = allLefts + allWidths;
allTops = allBottoms + allHeights;

% Get the four corners
minLeft = min(allLefts);
minBottom = min(allBottoms);
maxRight = max(allRights);
maxTop = max(allTops);

% Get the original width and height
widthOrig = maxRight - minLeft;
heightOrig = maxTop - minBottom;

% Get the approximate outer position of an axes covering all subplots
originalTotalPosition = [minLeft, minBottom, widthOrig, heightOrig];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
