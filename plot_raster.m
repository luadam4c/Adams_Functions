function [hRaster, eventTimes, yEnds, yTicksTable] = plot_raster (data, varargin)
%% Make a raster plot from a cell array of event time arrays
% Usage: [hRaster, eventTimes, yEnds, yTicksTable] = plot_raster (data, varargin)
% Explanation:
%       Plots a raster plot colored by groups.
%       Each group has a different set of line handles
%
% Example(s):
%       data = {magic(3), 5, (1:5)'};
%       [hRaster, eventTimes, yEnds] = ...
%           plot_raster(data, 'VertBarWidth', 0.6, ...
%                       'LineStyle', '-', 'LineWidth', 2, ...
%                       'ColorMap', {'Blue', 'Red', 'Purple'}, ...
%                       'Labels', {'3 vectors', '1 number', '1 vector'}, ...
%                       'HorzBarWindows', {[1, 5], [1, 4], [1, 3], [1, 2], [1, 1]});
%
% Outputs:
%       hRaster      - handles to the lines for each group
%                   specified as cell array of vectors of primitive line objects
%       eventTimes  - the event times for each group, linearized
%                   specified as a cell array of numeric row vectors
%       yEnds       - the y value endpoints for the bars of each group 
%                       (each column corresponds to an event time)
%                   specified as a cell array of numeric arrays with 2 rows
%       yTicksTable - a table with two fields:
%                       locs    - Y tick values
%                       labels  - Y tick labels
%
% Side Effects:
%       Plots a raster plot colored by groups.
%
% Arguments:    
%       data        - event time arrays
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'HorzBarWindows': horizontal bar windows
%                   must be empty or a cell array of numeric vectors
%                           with the same length as nVectors
%                   default == []
%                   - 'YMid': y value midpoints
%                   must be a numeric vector
%                   default == use trial numbers in reverse order
%                   - 'VertBarWidth': vertical bar width relative to 
%                                   y value increments (0~1)
%                   must be a positive scalar
%                   default == 0.6
%                   - 'ColorMap': a color map for the arrays
%                       Note: By default, each array has a different color
%                   must be a numeric array with 3 columns
%                   default == decide_on_colormap(colorMap, nArrays)
%                   - 'Labels': labels for each array
%                   must be a cell array of character/string arrays
%                   default == array numbers
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses 1 more than max and min of trialNos
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                   default == ['Time (', xUnits, ')']
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                   default == 'suppress'
%                   - 'YTickLocs': locations of Y ticks
%                   must be 'suppress' or a numeric vector
%                   default == ntrials:1
%                   - 'YTickLabels': labels for each raster
%                   must be 'suppress' or a cell array of character/string arrays
%                   default == trial numbers
%                   - 'LegendLocation': location for legend
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'      - use default
%                       'suppress'  - no legend
%                       anything else recognized by the legend() function
%                   default == 'suppress' if nPlots == 1 
%                               'northeast' if nPlots is 2~9
%                               'eastoutside' if nPlots is 10+
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be a positive integer scalar
%                   default == []
%                   - 'FigExpansion': expansion factors for figure position
%                   must be a must be a positive scalar or 2-element vector
%                   default == []
%                   - 'ClearFigure': whether to clear figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == falss
%                   - 'AlwaysNew': whether to always create a new figure even if
%                                   figNumber is not passed in
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       /home/Matlab/Downloaded_Functions/rgb.m
%       cd/apply_iteratively.m
%       cd/create_labels_from_numbers.m
%       cd/create_error_for_nargin.m
%       cd/count_vectors.m
%       cd/extract_common_prefix.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%       cd/set_figure_properties.m
%
% Used by:    
%       cd/parse_multiunit.m
%       cd/plot_relative_swd_raster.m
%       cd/plot_swd_raster.m
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m

% File History:
% 2018-05-16 Created by Adam Lu
% 2018-12-18 Now uses iP.KeepUnmatched
% 2019-02-23 Fixed bugs
% 2019-02-23 Added 'YLimits' as an optional argument
% 2019-02-24 Added maxNYTicks
% 2019-02-25 Added 'HorzBarWindows' as an optional argument
% 2019-03-14 Fixed the case when there is a condition with no spikes
% 2019-09-11 Updated 'Colors' to 'ColorMap'
% 2019-09-11 Updated to not use parfor
% TODO: Distinguish plot_raster.m vs plot_raster_plot.m?
% 

%% Hard-coded parameters
maxNYTicks = 20;             % maximum number of Y ticks

%% Default values for optional arguments
horzBarWindowsDefault = [];     % no horizontal bars by default
yMidDefault = [];               % set later
vertBarWidthDefault = 0.6;      % default bar width relative to y value increments
colorMapDefault = [];           % set later
lineStyleDefault = '-';         % default line style of bars
lineWidthDefault = 2;           % default line width of bars
labelsDefault = {};             % default labels to use for each array
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later
xUnitsDefault = 'unit';         % the default x-axis units
xLabelDefault = '';             % set later
yLabelDefault = '';             % set later
yTickLocsDefault = [];          % set later
yTickLabelsDefault = {};        % set later
legendLocationDefault = 'auto'; % set later
figTitleDefault = '';           % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figExpansionDefault = [];       % no figure expansion by default
clearFigureDefault = false;     % don't clear figure by default
alwaysNewDefault = false;       % don't always create new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If not compiled, add directories to search path for required functions
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for Downloaded_Functions
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
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
addRequired(iP, 'data', ...             % a cell array of event time arrays
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['data must be either empty or a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'HorzBarWindows', horzBarWindowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['HorzBarWindows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'YMid', yMidDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'VertBarWidth', vertBarWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'Labels', labelsDefault, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'Labels must be ''suppress'' or a cell array of character/string arrays!'));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XUnits', xUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'YTickLocs', yTickLocsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || isnumericvector(x), ...
        'YTickLocs must be ''suppress'' or a numeric vector!'));
addParameter(iP, 'YTickLabels', yTickLabelsDefault, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || ...
                iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'YTickLabels must be ''suppress'' or a cell array of character/string arrays!'));
addParameter(iP, 'LegendLocation', legendLocationDefault, ...
    @(x) all(islegendlocation(x, 'ValidateMode', true)));
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
addParameter(iP, 'AlwaysNew', alwaysNewDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
horzBarWindows = iP.Results.HorzBarWindows;
yMidUser = iP.Results.YMid;
vertBarWidth = iP.Results.VertBarWidth;
colorMap = iP.Results.ColorMap;
labels = iP.Results.Labels;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
xUnits = iP.Results.XUnits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
yTickLocs = iP.Results.YTickLocs;
yTickLabels = iP.Results.YTickLabels;
[~, legendLocation] = islegendlocation(iP.Results.LegendLocation, ...
                                        'ValidateMode', true);
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figExpansion = iP.Results.FigExpansion;
clearFigure = iP.Results.ClearFigure;
alwaysNew = iP.Results.AlwaysNew;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Prepare for plotting
% If numeric, force as a column cell array of column vectors
if isnumeric(data)
    data = force_column_vector(data, 'IgnoreNonVectors', false, ...
                                'ForceCellOutput', true);
end

% Force measure bars to be column vectors
if ~isempty(horzBarWindows)
    horzBarWindows = force_column_vector(horzBarWindows, ...
                        'IgnoreNonVectors', false, 'ForceCellOutput', true);
end

% Get the number of event time arrays to plot
nArrays = numel(data);

% If there is nothing to plot, return
if nArrays == 0
    fprintf('Threre is nothing to plot!!\n');
    hRaster = [];
    eventTimes = [];
    yEnds = [];
    yTicksTable = [];
    return
end

% Create a color map for the arrays based on either the rgb function
%   if the colors are provided, or based on the built-in parula map if not
colorMap = decide_on_colormap(colorMap, nArrays);

% Create labels if not provided
if isempty(labels)
    labels = create_labels_from_numbers(nArrays, 'Prefix', 'Group #');
end

% Get the number of vectors in each event time array
nVectors = cellfun(@count_vectors, data);

% If there is an empty vector, consider it a vector
nVectors(nVectors == 0) = 1;

% Get the total number of vectors
nVectorsAll = sum(nVectors);

% Assign trial numbers grouped by each event time array 
trialNos = cell(size(data));
ct = 0;                                 % counts number of trials assigned
for iArray = 1:nArrays
    % Get the number of vectors in this array
    nVectorThis = nVectors(iArray);

    % Assign the trial numbers in ascending order
    trialNosThis = ct + transpose(1:nVectorThis);

    % Store the trial numbers for this event time array
    trialNos{iArray} = trialNosThis;

    % Update number of trials assigned
    ct = ct + nVectorThis;
end

% Ungroup trial numbers
trialNosAll = vertcat(trialNos{:});

% Compute default y midpoints grouped by each event time array
if isempty(yMidUser)
    yMids = cellfun(@(x) nVectorsAll - x + 1, trialNos, ...
                    'UniformOutput', false);
    yMidsAll = vertcat(yMids{:});
else
    yMidsAll = yMidUser;
    yMids = cellfun(@(x) yMidUser(x), trialNos, ...
                    'UniformOutput', false);
end

% Compute y increment
if numel(yMidsAll) == 1
    yIncr = 1;
else
    yIncr = abs(yMidsAll(2) - yMidsAll(1));
end

% Decide on indices for ticks if total number of vectors exceed maxNYTicks
if nVectorsAll > maxNYTicks
    indTicks = create_indices([1, nVectorsAll], 'MaxNum', maxNYTicks);
end

% Decide on Y tick values
if ~ischar(yTickLocs) || ~strcmpi(yTickLocs, 'suppress')
    if ~isempty(yTickLocs)
        % If provided, use custom Y tick values instead
        yTicks.locs = yTickLocs;
    else
        % Set the Y tick values at the midpoints
        if nVectorsAll <= maxNYTicks
            yTicks.locs = yMidsAll;
        else
            yTicks.locs = yMidsAll(indTicks);            
        end
    end
end

% Decide on Y tick labels
if ~ischar(yTickLabels) || ~strcmpi(yTickLabels, 'suppress')
    if ~isempty(yTickLabels)
        % If provided, use custom Y tick labels instead
        yTicks.labels = yTickLabels;
    else
%        if ~isempty(labels)
            % TODO: Not Implemented yet!
%        else
            % Create trial labels from numbers
            trialLabels = create_labels_from_numbers(trialNosAll);

            % Use the trial numbers
            if nVectorsAll <= maxNYTicks
                yTicks.labels = trialLabels;
            else
                yTicks.labels = trialLabels(indTicks);            
            end
%        end
    end
end

% Convert yTicks to a Y tick table
if ~ischar(yTickLabels) || ~strcmpi(yTickLabels, 'suppress')
    % Convert yTicks to a Y tick table
    yTicksTable = struct2table(yTicks);

    % Sort the Y tick table according to Y tick values
    yTicksTable = sortrows(yTicksTable, 'locs');
end

% Get the half bar width in actual coordinates
halfBarWidth = (vertBarWidth / 2) * yIncr;

% Compute the y values for the horizontal lines
if ~isempty(horzBarWindows)
    yHorzBars = yMidsAll + halfBarWidth;
end

% Assign y value endpoints to each event time
[eventTimes, yEnds] = ...
    cellfun(@(x, y) compute_y_endpoints(x, y, halfBarWidth), ...
            data, yMids, 'UniformOutput', false);

% Set the default x axis limits
if isempty(xLimits)
    xLimits = 'suppress';
end

% Set the default y axis limits
if isempty(yLimits)
    maxTrialNo = apply_iteratively(@max, trialNos);
    minTrialNo = apply_iteratively(@min, trialNos);
    yLimits = [minTrialNo - 1, maxTrialNo + 1];
end

% Set the default x-axis label
if isempty(xLabel)
    xLabel = ['Time (', xUnits, ')'];
end

% Set the default x-axis label
if isempty(yLabel)
    yLabel = 'suppress';
end

% Set the default legend location based on number of arrays
if strcmpi(legendLocation, 'auto')
    legendLocation = 'suppress';
end

% Set the default figure title
if isempty(figTitle)
    if ~isempty(labels)
        commonPrefix = extract_common_prefix(labels, 'Delimiter', '_');
        if ~isempty(commonPrefix)
            figTitle = ['Event times for ', commonPrefix];
        else
            figTitle = 'Event times';
        end
    else
        figTitle = 'Event times';
    end
end

%% Plot the event time arrays
% Decide on the figure to plot on
fig = set_figure_properties('FigHandle', figHandle, 'FigNumber', figNumber, ...
                'FigExpansion', figExpansion, 'ClearFigure', clearFigure, ...
                'AlwaysNew', alwaysNew);

% Plot the event time arrays
hRaster = cell(size(data));
for iArray = 1:nArrays
    % Get the event times and y endpoints
    eventTimesThis = eventTimes{iArray};
    yEndsThis = yEnds{iArray};

    % Get the color for this array
    colorThis = colorMap(iArray, :);

    % Get the label for this array
    labelThis = labels{iArray};

    % Plot the event times with the color for this array
    hRaster{iArray} = line(eventTimesThis, yEndsThis, ...
                            'LineStyle', lineStyleDefault, ...
                            'LineWidth', lineWidthDefault, ...
                            'Color', colorThis, ...
                            'DisplayName', labelThis, otherArguments);
end

% Plot horizontal line(s) for duration if provided
if ~isempty(horzBarWindows)
    hold on;
    hBars = cellfun(@(x, y) plot_horizontal_line(x, 'XLimits', y, ...
                        'Color', 'r', 'LineStyle', '-', 'LineWidth', 0.5), ...
                        num2cell(yHorzBars), horzBarWindows, ...
                        'UniformOutput', false);
end

% Change the y tick values and labels
if ~ischar(yTickLabels) || ~strcmpi(yTickLabels, 'suppress')
    set(gca, 'YTick', yTicksTable.locs);
    set(gca, 'YTickLabel', yTicksTable.labels);
end

% Set time axis limits
if ~ischar(xLimits) || ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~ischar(yLimits) || ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Generate an x-axis label
if ~strcmpi(xLabel, 'suppress')
    xlabel(xLabel);
end

% Generate a y-axis label
if ~strcmpi(yLabel, 'suppress')
    ylabel(yLabel);
end

% Generate a legend if there is more than one trace
if ~strcmpi(legendLocation, 'suppress')
%    legend('Location', legendLocation);
% TODO: Use only the first line Object of each set
    error('Not Implemented yet!')
end

% Generate a title
if ~strcmpi(figTitle, 'suppress')
    title(figTitle);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventTimes, yEnds] = ...
                compute_y_endpoints (timesArray, yMids, halfBarWidth);

% Reshape the time values as a single column 
timesVector = timesArray(:);

% Duplicate the time columns, then transpose
eventTimes = transpose([timesVector, timesVector]);

% Get the dimensions of the event time array
sizeThis = size(timesArray);

% Assign y value midpoints to all event times
yMidsArray = ones(sizeThis) * diag(yMids);

% Reshape as a single column
yMidsVector = yMidsArray(:);

% Shift the y midpoints by half the bar width to get the y endpoints,
%   then transpose it so that each column corresponds to a time point
yEnds = transpose(yMidsVector + halfBarWidth * [-1, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
