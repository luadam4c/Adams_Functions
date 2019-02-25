function [hLines, eventTimes, yEnds, yTicksTable] = plot_raster (data, varargin)
%% Make a raster plot from a cell array of event time arrays
% Usage: [hLines, eventTimes, yEnds, yTicksTable] = plot_raster (data, varargin)
% Explanation:
%       Plots a raster plot colored by groups.
%       Each group has a different set of line handles
% Example(s):
%       data = {magic(3), 5, (1:5)'};
%       [hLines, eventTimes, yEnds] = ...
%           plot_raster(data, 'BarWidth', 0.6, ...
%                       'LineStyle', '-', 'LineWidth', 2, ...
%                       'Colors', {'Blue', 'Red', 'Purple'}, ...
%                       'Labels', {'3 vectors', '1 number', '1 vector'});
% Outputs:
%       hLines      - handles to the lines for each group
%                   specified as cell array of vectors of primitive line objects
%       eventTimes  - the event times for each group, linearized
%                   specified as a cell array of numeric row vectors
%       yEnds       - the y value endpoints for the bars of each group 
%                       (each column corresponds to an event time)
%                   specified as a cell array of numeric arrays with 2 rows
%       yTicksTable - a table with two fields:
%                       locs    - Y tick values
%                       labels  - Y tick labels
% Side Effects:
%       Plots a raster plot colored by groups.
% Arguments:    
%       data        - event time arrays
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'YMid': y value midpoints
%                   must be a numeric vector
%                   default == use trial numbers in reverse order
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses 1 more than max and min of trialNos
%                   - 'BarWidth': bar width relative to y value increments (0~1)
%                   must be a positive scalar between 0 and 1
%                   default == 0.6
%                   - 'Colors': colors for each array
%                   must be a cell array of character/string arrays
%                   default == equal samples from the parula map
%                   - 'Labels': labels for each array
%                   must be a cell array of character/string arrays
%                   default == array numbers
%                   - 'YTickLocs': locations of Y ticks
%                   must be 'suppress' or a numeric vector
%                   default == ntrials:1
%                   - 'YTickLabels': labels for each raster
%                   must be 'suppress' or a cell array of character/string arrays
%                   default == trial numbers
%                   - Any other parameter-value pair for the line() function
%
% Requires:
%       /home/Matlab/Downloaded_Functions/rgb.m
%       cd/apply_iteratively.m
%       cd/create_labels_from_numbers.m
%       cd/create_error_for_nargin.m
%       cd/count_vectors.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%
% Used by:    
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m
%
% File History:
% 2018-05-16 Created by Adam Lu
% 2018-12-18 Now uses iP.KeepUnmatched
% 2019-02-23 Fixed bugs
% 2019-02-23 Added 'YLimits' as an optional argument
% 2019-02-24 Added maxYTicks
% 

%% Hard-coded parameters
maxYTicks = 20;             % maximum number of Y ticks

%% Default values for optional arguments
yMidDefault = [];           % set later
yLimitsDefault = [];        % set later
barWidthDefault = 0.6;      % default bar width relative to y value increments
colorsDefault = {};         % default colors to use for each array
labelsDefault = {};         % default labels to use for each array
lineStyleDefault = '-';     % default line style of bars
lineWidthDefault = 2;       % default line width of bars

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
addParameter(iP, 'YMid', yMidDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'BarWidth', barWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'Colors', colorsDefault, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'data must be a cell array of character/string arrays!'));
addParameter(iP, 'Labels', labelsDefault, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'Labels must be ''suppress'' or a cell array of character/string arrays!'));
addParameter(iP, 'YTickLocs', [], ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || isnumericvector(x), ...
        'YTickLocs must be ''suppress'' or a numeric vector!'));
addParameter(iP, 'YTickLabels', {}, ...
    @(x) assert(ischar(x) && strcmpi(x, 'suppress') || ...
                iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'YTickLabels must be ''suppress'' or a cell array of character/string arrays!'));

% Read from the Input Parser
parse(iP, data, varargin{:});
yMidUser = iP.Results.YMid;
yLimits = iP.Results.YLimits;
barWidth = iP.Results.BarWidth;
colors = iP.Results.Colors;
labels = iP.Results.Labels;
yTickLocs = iP.Results.YTickLocs;
yTickLabels = iP.Results.YTickLabels;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Prepare for plotting
% If numeric, force as a column cell array of column vectors
%   Note: if data is already a cell array, don't use this 
%           function or else elements can't be non-vectors!
if isnumeric(data)
    data = force_column_vector(data, 'IgnoreNonVectors', false, ...
                                'ForceCellOutput', true);
end

% Get the number of event time arrays to plot
nArrays = numel(data);

% If there is nothing to plot, return
if nArrays == 0
    fprintf('Threre is nothing to plot!!\n');
    hLines = [];
    eventTimes = [];
    yEnds = [];
    yTicksTable = [];
    return
end

% Create a color map for the arrays based on either the rgb function
%   if the colors are provided, or based on the built-in parula map if not
if isempty(colors)
    % Create a color map based on the built-in parula map
    cm = colormap(parula(nArrays));

    % Place the colors in a cell array
    colors = cell(size(data));
    parfor iArray = 1:nArrays
        colors{iArray} = cm(iArray, :);
    end
else
    % Use the rgb function to convert color-strings into a 3-element array
    colors = cellfun(@rgb, colors, 'UniformOutput', false);
end

% Create labels if not provided
if isempty(labels)
    parfor iArray = 1:nArrays
        labels{iArray} = sprintf('Group #%d', iArray);
    end
end

% Force vectors to be column vectors
parfor iArray = 1:nArrays
    % Get the current array
    dataThis = data{iArray};

    % If it is a vector, force it to be a column
    if isvector(dataThis)
        data{iArray} = dataThis(:);
    end
end

% Get the number of vectors in each event time array
nVectors = cellfun(@count_vectors, data);

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
    yIncr = yMidsAll(2) - yMidsAll(1);
end

% Decide on indices for ticks if total number of vectors exceed maxYTicks
if nVectorsAll > maxYTicks
    indTicks = create_indices([1, nVectorsAll], 'MaxNum', maxYTicks);
end

% Decide on Y tick values
if ~ischar(yTickLocs) || ~strcmpi(yTickLocs, 'suppress')
    if ~isempty(yTickLocs)
        % If provided, use custom Y tick values instead
        yTicks.locs = yTickLocs;
    else
        % Set the Y tick values at the midpoints
        if nVectorsAll <= maxYTicks
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
        % Create trial labels from numbers
        trialLabels = create_labels_from_numbers(trialNosAll);

        % Use the trial numbers
        if nVectorsAll <= maxYTicks
            yTicks.labels = trialLabels;
        else
            yTicks.labels = trialLabels(indTicks);            
        end
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
halfBarWidth = (barWidth / 2) * yIncr;

% Assign y value endpoints to each event time
eventTimes = cell(size(data));
yEnds = cell(size(data));
parfor iArray = 1:nArrays
    % Get the time values for this event time array
    timesArray = data{iArray};

    % Reshape as a single column 
    timesVector = timesArray(:);

    % Duplicate the time columns, then transpose
    eventTimesThis = [timesVector, timesVector]';

    % Get the dimensions of the event time array
    sizeThis = size(timesArray);

    % Get the y value midpoints
    yMidsThis = yMids{iArray};

    % Assign y value midpoints to all event times
    yMidsArray = ones(sizeThis) * diag(yMidsThis);

    % Reshape as a single column 
    yMidsVector = yMidsArray(:);

    % Shift the y midpoints by half the bar width to get the y endpoints,
    %   then transpose it so that each column corresponds to a time point
    yEndsThis = [yMidsVector - halfBarWidth, yMidsVector + halfBarWidth]';

    % Store the event times and y endpoints
    eventTimes{iArray} = eventTimesThis;
    yEnds{iArray} = yEndsThis;
end

% Compute y axis limits
if isempty(yLimits)
    maxTrialNo = apply_iteratively(@max, trialNos);
    minTrialNo = apply_iteratively(@min, trialNos);
    yLimits = [minTrialNo - 1, maxTrialNo + 1];
end

%% Plot the event time arrays
hLines = cell(size(data));
for iArray = 1:nArrays
    % Get the event times and y endpoints
    eventTimesThis = eventTimes{iArray};
    yEndsThis = yEnds{iArray};

    % Get the color for this array
    colorThis = colors{iArray};

    % Get the label for this array
    labelThis = labels{iArray};

    % Plot the event times with the color for this array
    hLines{iArray} = line(eventTimesThis, yEndsThis, ...
                            'LineStyle', lineStyleDefault, ...
                            'LineWidth', lineWidthDefault, ...
                            'Color', colorThis, ...
                            'DisplayName', labelThis, otherArguments);
end

% Change the y tick values and labels
if ~ischar(yTickLabels) || ~strcmpi(yTickLabels, 'suppress')
    set(gca, 'YTick', yTicksTable.locs);
    set(gca, 'YTickLabel', yTicksTable.labels);
end

% Set y axis limits
if ~ischar(yLimits) || ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% The following doesn't work!
data = cellfun(@(x) if isvector(x); x = x(:); end, ...
                data, 'UniformOutput', false); 

%                   - 'LineStyle': line style of bars
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - 'LineWidth': line width of bars
%                   must be a positive scalar
%                   default == 2
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
                            'LineStyle', lineStyle, ...
                            'LineWidth', lineWidth, ...
%       cd/islinestyle.m

yMids = cell(size(data));
% Store the y midpoints for this event time array
yMids{iArray} = yMidsThis;

% Assign trial numbers to each event time array 
%   and save Y tick values and labels
yTicks.locs = zeros(nVectorsAll, 1);
yTicks.labels = cell(nVectorsAll, 1);
for iArray = 1:nArrays
    % Convert to strings for the Y tick labels
    for iTrial = trialNosThis
        yTicks.labels{iTrial} = num2str(iTrial);
    end

    % Set the Y tick values at the midpoints
    yTicks.locs(trialNosThis) = yMidsThis;
end

trialLabels = arrayfun(@num2str, trialNosAll, 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
