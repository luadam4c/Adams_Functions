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
%       data        - a cell array of event time arrays
%                   must be a cell array of numeric arrays
%       varargin    - 'BarWidth': bar width relative to y value increments (0~1)
%                   must be a positive scalar between 0 and 1
%                   default == 0.6
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
%                   - 'Colors': colors for each array
%                   must be a cell array of character/string arrays
%                   default == equal samples from the parula map
%                   - 'Labels': labels for each array
%                   must be a cell array of character/string arrays
%                   default == array numbers
%                   - 'YTickLocs': locations of Y ticks
%                   must be a numeric vector
%                   default == ntrials:1
%                   - 'YTickLabels': labels for each raster
%                   must be a cell array of character/string arrays
%                   default == trial numbers
%
% Requires:
%       /home/Matlab/Downloaded_Functions/rgb.m
%       /home/Matlab/Adams_Functions/islinestyle.m
%
% Used by:    
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m
%
% File History:
% 2018-05-16 Created by Adam Lu
% 

%% Default values for optional arguments
barWidthDefault = 0.6;      % default bar width relative to y value increments
lineStyleDefault = '-';     % default line style of bars
lineWidthDefault = 2;       % default line width of bars
colorsDefault = {};         % default colors to use for each array
labelsDefault = {};         % default labels to use for each array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'data', ...             % a cell array of event time arrays
    @(x) assert(iscell(x) && all(cellfun(@isnumeric, x)), ...
        'data must be a cell array of numeric arrays!'));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BarWidth', barWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Colors', colorsDefault, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'data must be a cell array of character/string arrays!'));
addParameter(iP, 'Labels', labelsDefault, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'data must be a cell array of character/string arrays!'));
addParameter(iP, 'YTickLocs', [], ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'YTickLabels', {}, ...
    @(x) assert(iscell(x) && all(cellfun(@(x) ischar(x) || isstring(x), x)), ...
        'data must be a cell array of character/string arrays!'));

% Read from the Input Parser
parse(iP, data, varargin{:});
barWidth = iP.Results.BarWidth;
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
lineWidth = iP.Results.LineWidth;
colors = iP.Results.Colors;
labels = iP.Results.Labels;
yTickLocs = iP.Results.YTickLocs;
yTickLabels = iP.Results.YTickLabels;

%% Prepare for plotting
% Get the number of event time arrays to plot
nArrays = numel(data);

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
nVectors = cellfun(@(x) size(x, 2), data);

% Get the total number of vectors
nVectorsAll = sum(nVectors);

% Assign trial numbers to each event time array 
%   and save Y tick values and labels
trialNos = cell(size(data));
yMids = cell(size(data));
yTicks.locs = zeros(nVectorsAll, 1);
yTicks.labels = cell(nVectorsAll, 1);
ct = 0;                                 % counts number of trials assigned
for iArray = 1:nArrays
    % Get the number of vectors in this array
    nVectorThis = nVectors(iArray);

    % Assign the trial numbers in ascending order
    trialNosThis = ct + (1:nVectorThis);

    % Convert to strings for the Y tick labels
    for iTrial = trialNosThis
        yTicks.labels{iTrial} = num2str(iTrial);
    end

    % Convert trial numbers to y value midpoints
    yMidsThis = nVectorsAll - trialNosThis + 1;

    % Set the Y tick values at the midpoints
    yTicks.locs(trialNosThis) = yMidsThis;

    % Store the trial numbers and y midpoints for this event time array
    trialNos{iArray} = trialNosThis;
    yMids{iArray} = yMidsThis;

    % Update number of trials assigned
    ct = ct + nVectorThis;
end

% If provided, use custom Y tick values instead
if ~isempty(yTickLocs)
    yTicks.locs = yTickLocs;
end

% If provided, use custom Y tick labels instead
if ~isempty(yTickLabels)
    yTicks.labels = yTickLabels;
end

% Convert yTicks to a Y tick table
yTicksTable = struct2table(yTicks);

% Sort the Y tick table according to Y tick values
yTicksTable = sortrows(yTicksTable, 'locs');

% Get the half bar width
halfBarWidth = barWidth / 2;

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
                            'LineStyle', lineStyle, ...
                            'LineWidth', lineWidth, ...
                            'Color', colorThis, ...
                            'DisplayName', labelThis);
end

% Change the y tick values and labels
set(gca, 'YTick', yTicksTable.locs);
set(gca, 'YTickLabel', yTicksTable.labels);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% The following doesn't work!
data = cellfun(@(x) if isvector(x); x = x(:); end, ...
                data, 'UniformOutput', false); 

%}
