function h = plot_pulse (timeVec, pulseVecs, varargin)
%% Plots pulses, marks the parsed endpoints and displays total number of sweeps
% Usage: h = plot_pulse (timeVec, pulseVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:    
%       timeVec     - time vector(s)
%                   must be a numeric array
%       pulseVecs   - pulse vector(s)
%                   must be a numeric array
%       varargin    - 'PulseParams': a table containing the parsed parameters
%                   must be a table containing at least the fields:
%                       baseValue
%                       pulseValue
%                       idxBeforeStart
%                       idxAfterStart
%                       idxBeforeEnd
%                       idxAfterEnd
%                   default == returned by parse_pulse.m
%                   - 'ToUse': whether each vector will be used
%                   must be a logical/numeric binary vector
%                   default == all will be used
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(timeVec), max(timeVec)]
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == expand by a little bit
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/parse_pulse.m
%
% Used by:    
%       cd/find_passive_params.m

% File History:
% 2018-10-12 Adapted from code in find_passive_params.m
% 2018-12-19 Now uses unmatched varargin parts as parameters for plot()
% TODO: Use improved version of match_vector_numbers.m

%% Hard-coded parameters
yCoverage = 80;             % coverage of y axis (%)

%% Default values for optional arguments
pulseParamsDefault = [];    % set later
toUseDefault = [];          % set later
xLimitsDefault = [];        % set later
yLimitsDefault = [];        % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'timeVec', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'pulseVecs', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PulseParams', pulseParamsDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'ToUse', toUseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, timeVec, pulseVecs, varargin{:});
pulseParams = iP.Results.PulseParams;
toUse = iP.Results.ToUse;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% Count the number of vectors
nVectors = size(pulseVecs, 2);

% Parse pulseParams if not provided
if isempty(pulseParams)
    pulseParams = parse_pulse(pulseVecs);
end

% Use all vectors by default
if isempty(toUse)
    toUse = true(nVectors, 1);
end

% Count the number of vectors to use
nToUse = sum(toUse);

% Set default x limits to the entire time range
if isempty(xLimits)
    xLimits = [min(timeVec), max(timeVec)];
end

% Set default y limits to the entire time range
if isempty(yLimits)
    % Extract all base and pulse values
    baseValue = pulseParams.baseValue;
    pulseValue = pulseParams.pulseValue;

    % Restrict to those to use
    baseValue = baseValue(toUse);
    pulseValue = pulseValue(toUse);

    % Compute the y limits
    yLimits = compute_axis_limits([baseValue; pulseValue], 'y', ...
                                    'Coverage', yCoverage);
end

% Extract from params
idxBeforeStart = pulseParams.idxBeforeStart;
idxAfterStart = pulseParams.idxAfterStart;
idxBeforeEnd = pulseParams.idxBeforeEnd;
idxAfterEnd = pulseParams.idxAfterEnd;

% Force vectors into cell arrays
pulseVecsCell = force_column_cell(pulseVecs);

% Vertically concatenate two sets
pulseVecsCellRepeated = repmat(pulseVecsCell, [2, 1]);

% Generate the indices for left and right triangle markers
indLeftTriangles = [idxBeforeStart; idxBeforeEnd];
indRightTriangles = [idxAfterStart; idxAfterEnd];

% Generate the x and y locations for markers
xLeftTriangles = timeVec(indLeftTriangles);
yLeftTriangles = cellfun(@(x, y) x(y), ...
                            pulseVecsCellRepeated, num2cell(indLeftTriangles));
xRightTriangles = timeVec(indRightTriangles);
yRightTriangles = cellfun(@(x, y) x(y), ...
                            pulseVecsCellRepeated, num2cell(indRightTriangles));

%% Do the job
% Hold on
hold on;

% Plot vectors
if nToUse == nVectors
    % Plot all vectors at once
    plot(timeVec, pulseVecs, '-');
else
    % Plot vectors
    for iVec = 1:nVectors
        if toUse(iVec)
            % Plot vectors as solid line if used
            plot(timeVec, pulseVecs(:, iVec), '-', otherArguments);
        else
            % Plot vectors as dotted line if not used
            plot(timeVec, pulseVecs(:, iVec), '--', otherArguments);
        end
    end
end

% Plot detected endpoints
plot(xLeftTriangles, yLeftTriangles, '<', otherArguments);
plot(xRightTriangles, yRightTriangles, '>', otherArguments);

% Set time axis limits
if ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Get axes ranges
xRange = xLimits(2) - xLimits(1);
yRange = yLimits(2) - yLimits(1);

% Display total number of sweeps
xpos = xLimits(1) + (1/20) * xRange;
ypos = yLimits(1) + (19/20) * yRange;
text(xpos, ypos, ['Total number of sweeps: ', num2str(nVectors)]);

% Label the axes
xlabel('Time (ms)');
ylabel('Current (pA)');

% TODO
h = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract all pulse amplitudes
pulseAmplitude = pulseParams.pulseAmplitude;

% Restrict to those to use
pulseAmplitudeToUse = pulseAmplitude(toUse);

% Take the average
meanPulseAmplitude = mean(pulseAmplitudeToUse);

% Set the y axis limits
yLimits = [meanPulseAmplitude*6/5, -meanPulseAmplitude*1/5];

% Get axes limits
ax = gca;
xLimits = get(ax, 'Xlim');
yLimits = get(ax, 'Ylim');

% Hold on if plotting multiple vectors
if nVectors > 1
    hold on;
end

if ~isempty(firstDipPoint{iVec})
    plot(timeVec(indBaseline(iVec, 1)), pulseVecs(indBaseline(iVec, 1), iVec), '>');
    plot(timeVec(indBaseline(iVec, end)), pulseVecs(indBaseline(iVec, end), iVec), '<');
    plot(timeVec(indSteady(iVec, 1)), pulseVecs(indSteady(iVec, 1), iVec), '>');
    plot(timeVec(indSteady(iVec, end)), pulseVecs(indSteady(iVec, end), iVec), '<');
end

% Find the extremes
maxValue = max([baseValue; pulseValue]);
minValue = min([baseValue; pulseValue]);

% Compute the y limits
yLimits = compute_ylimits(minValue, maxValue, 'Coverage', yCoverage);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%