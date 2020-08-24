function hPlot = minEASE_plot_gapfree_events (tVec, dataRaw, dataLowpass, ...
                                eventInfo, eventClass, outputLabel, varargin)
%% Plots all detected events 
% Usage: hPlot = minEASE_plot_gapfree_events (tVec, dataRaw, dataLowpass, ...
%                               eventInfo, eventClass, outputLabel, varargin)
% Explanation:    
%       TODO
% Outputs:    
%       TODO
%
% Arguments:    
%       TODO
%
% Requires:
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:
%       cd/minEASE_detect_gapfree_events.m
%       cd/minEASE_gui_examine_events.m
%
% File History:
% 2017-06-06 AL - Moved from minEASE_detect_gapfree_events.m
% 2017-06-07 AL - Changed classColorMap to use rgb() instead of colormap
% 2017-06-07 AL - Added Class 8
% 2017-06-09 AL - Added hPlot to manage groups of plots
% 2017-06-09 AL - Made class numbers, colors and labels optional arguments
% 2017-06-09 AL - Fixed time vector when not starting from zero
% 2017-06-13 AL - Account for the case that no events exist
% 2017-10-16 AL - Change Slow Decay -> Wrong Decay
% 2018-01-28 AL - Added noTitleFlag
% 2018-01-28 AL - Added isdeployed
% 2018-01-28 AL - Now checks if eventInfo is empty
% 2018-02-23 AL - Added 'HitTest', 'off' to all plots
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
% 2018-08-03 AL - Updated legend to turn 'AutoUpdate' off for R2017a and beyond

%% Class numbers, colors and labels
classNumbersDefault = 1:8;
classColorsDefault = {'Lime', 'Red', 'Gold', ...
                      'Cyan', 'DarkOrange', 'Purple', ...
                      'Violet', 'Gray'};
classLabelsDefault = {'Class 1: Type I PSC', 'Class 2: Type II PSC', ...
                      'Class 3: Type III PSC', 'Class 4: Wrong Decay', ...
                      'Class 5: Slow Rise', 'Class 6: Too Small', ...
                      'Class 7: In Seal Test', 'Class 8: Removed'};

%% Argument defaults
dataSmoothDefault = [];                     % TODO
legendLocationDefault = 'eastoutside';      % TODO
lineWidthDefault = 0.5;                     % TODO
markerSizeDefault = 3;                      % TODO
isCheckedDefaultTemp = [];                  % TODO
noTitleFlagDefault = false;                 % TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 6
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ClassNumbers', classNumbersDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));
addParameter(iP, 'ClassColors', classColorsDefault, ...
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
                || min(cellfun(@isstring, x))), ...
                ['ClassColors must be a cell array ', ...
                'of strings or character arrays!']));
addParameter(iP, 'ClassLabels', classLabelsDefault, ...
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
                || min(cellfun(@isstring, x))), ...
                ['ClassLabels must be a cell array ', ...
                'of strings or character arrays!']));
addParameter(iP, 'DataSmooth', dataSmoothDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));  
                                    % TODO: validate in a list of all possible legend locations
addParameter(iP, 'MarkerSize', markerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));  
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));  
addParameter(iP, 'IsChecked', isCheckedDefaultTemp, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));  
addParameter(iP, 'NoTitleFlag', noTitleFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));  

% Read from the Input Parser
parse(iP, varargin{:});
classNumbers = iP.Results.ClassNumbers;
classColors = iP.Results.ClassColors;
classLabels = iP.Results.ClassLabels;
dataSmooth = iP.Results.DataSmooth;
legendLocation = iP.Results.LegendLocation;
markerSize = iP.Results.MarkerSize;
lineWidth = iP.Results.LineWidth;
isChecked = iP.Results.IsChecked;
noTitleFlag = iP.Results.NoTitleFlag;

% Check relationships between arguments
% TODO: lengths of classNumbers, classColors, classLabels must be the same

% Set dependent parameter defaults
if isempty(isChecked)
    isChecked = zeros(size(eventInfo, 1), 1);
end

%% Construct color map
nClass = numel(classNumbers);                        % total number of classes
classColorMap = zeros(nClass, 3);
for iClass = 1:nClass
    classColorMap(iClass, :) = rgb(classColors{iClass});
end
%classColorMap = flipud(colormap(jet(nClass)));

%% Plot data
% Hold on
hold on;

% Plot original data in black
hPlot.dataRaw = plot(tVec, dataRaw, '-k', 'LineWidth', 0.5, ...
                     'DisplayName', 'raw data', 'HitTest', 'off');

% Plot lowpass filtered data in blue
hPlot.dataLowpass = plot(tVec, dataLowpass, '-b', 'LineWidth', 0.5, ...
                         'DisplayName', 'lowpass filtered', 'HitTest', 'off');

% Plot moving-average-filtered data in magenta
if ~isempty(dataSmooth)
    hPlot.dataSmooth = ...
        plot(tVec, dataSmooth, '-m', 'LineWidth', 0.5, ...
        'DisplayName', 'moving average filtered', 'HitTest', 'off');
end

%% Plot events
% Plot all event breakpoints as black crosses
hPlot.eventBreaks = gobjects;
if ~isempty(eventInfo)
    h = plot(tVec(eventInfo(:, 1)), eventInfo(:, 3), 'kx', ...
             'MarkerSize', markerSize, 'LineWidth', lineWidth, ...
             'DisplayName', 'Breakpoints', 'HitTest', 'off');
    if ~isempty(h)
        hPlot.eventBreaks = h;
    end
end

% Plot event peaks as circles with color according to class
hPlot.eventPeaks = gobjects(1, nClass);
for iClass = 1:nClass
    if ~isempty(eventInfo)
        % Find peak times of this class of events
        timeThisClassPeak = tVec(eventInfo(eventClass == classNumbers(iClass), 2));

        % Find peak values of this class of events
        valThisClassPeak = eventInfo(eventClass == classNumbers(iClass), 4);

        % Find peak values of this class of events
        h = plot(timeThisClassPeak, valThisClassPeak, ...
                 'o', 'Color', classColorMap(iClass, :), ...
                 'MarkerSize', markerSize, 'LineWidth', lineWidth, ...
                 'DisplayName', classLabels{iClass}, ...
                 'HitTest', 'off');
        if ~isempty(h)
            hPlot.eventPeaks(iClass) = h;
        end
    end
end

% Display legend
hLegend = legend('Location', legendLocation);
set(hLegend, 'AutoUpdate', 'off');              % for R2017a and beyond

% Display ylabel
ylabel('Current (pA)');

% Display xlabel
xlabel('Time (ms)');

% Create title
if ~noTitleFlag
    title(sprintf('All detected events for %s', outputLabel), ...
            'FontSize', 14, 'Interpreter', 'none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%