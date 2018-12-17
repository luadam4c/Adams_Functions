function h = plot_pulse_response (timeVec, responseVecs, varargin)
%% Plots pulse responses, marks the parsed baseline and steady state bounds
% Usage: h = plot_pulse_response (timeVec, responseVecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       h           - figure handle for the created figure
%                   specified as a figure handle
% Arguments:
%       timeVec         - time vector(s)
%                       must be a numeric array
%       responseVecs    - response vector(s)
%                       must be a numeric array
%       varargin    - 'PulseVectors': vector that contains the pulse itself
%                   must be a numeric vector
%                   default == [] (not used)
%                   - 'SameAsPulse': whether always the same as 
%                                       the current pulse endpoints
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MeanValueWindowMs': window in ms for 
%                                           calculating mean values
%                   must be a positive scalar
%                   default == 0.5 ms
%                   - 'MinPeakDelayMs': minimum peak delay (ms)
%                               after the end of the pulse
%                   must be a positive scalar
%                   default == 0 ms
%                   - 'ResponseParams': a table containing the parsed parameters
%                   must be a table containing at least the fields:
%                       baseValue
%                       steadyValue
%                       idxBaseStart
%                       idxBaseEnd
%                       idxSteadyStart
%                       idxSteadyEnd
%                   default == returned by parse_pulse_response.m
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
%
% Requires:
%       cd/compute_ylimits.m
%       cd/parse_pulse_response.m
%
% Used by:    
%       cd/find_passive_params.m

% File History:
% 2018-10-12 Adapted from code in plot_pulse.m
% TODO: Use improved version of match_vector_numbers.m

%% Hard-coded parameters
yCoverage = 80;                 % coverage of y axis (%)

%% Default values for optional arguments
pulseVectorsDefault = [];       % don't use pulse vectors by default
sameAsPulseDefault = true;      % use pulse endpoints by default
meanValueWindowMsDefault = 0.5; % calculating mean values over 0.5 ms by default
minPeakDelayMsDefault = 0;      % no minimum peak delay by default
responseParamsDefault = [];     % set later
toUseDefault = [];              % set later
xLimitsDefault = [];            % set later
yLimitsDefault = [];            % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'timeVec', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addRequired(iP, 'responseVecs', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PulseVectors', pulseVectorsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'SameAsPulse', sameAsPulseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MeanValueWindowMs', meanValueWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'MinPeakDelayMs', minPeakDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'scalar'}));
addParameter(iP, 'ResponseParams', responseParamsDefault, ...
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
parse(iP, timeVec, responseVecs, varargin{:});
pulseVectors = iP.Results.PulseVectors;
sameAsPulse = iP.Results.SameAsPulse;
meanValueWindowMs = iP.Results.MeanValueWindowMs;
minPeakDelayMs = iP.Results.MinPeakDelayMs;
responseParams = iP.Results.ResponseParams;
toUse = iP.Results.ToUse;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;

%% Preparation
% Count the number of vectors
nVectors = size(responseVecs, 2);

% Parse responseParams if not provided
if isempty(responseParams)
    % Compute the sampling interval
    siMs = timeVec(2) - timeVec(1);

    % Parse the pulse response
    responseParams = parse_pulse_response(responseVecs, siMs, ...
                                'PulseVectors', pulseVectors, ...
                                'SameAsPulse', sameAsPulse, ...
                                'MeanValueWindowMs', baselineLengthMs, ...
                                'MinPeakDelayMs', minPeakDelayMs);
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
    baseValue = responseParams.baseValue;
    steadyValue = responseParams.steadyValue;

    % Restrict to those to use
    baseValue = baseValue(toUse);
    steadyValue = steadyValue(toUse);

    % Find the extremes
    maxValue = max([baseValue; steadyValue]);
    minValue = min([baseValue; steadyValue]);

    % Compute the y limits
    yLimits = compute_ylimits(minValue, maxValue, 'Coverage', yCoverage);
end

% Extract from params
idxBaseStart = responseParams.idxBaseStart;
idxBaseEnd = responseParams.idxBaseEnd;
idxSteadyStart = responseParams.idxSteadyStart;
idxSteadyEnd = responseParams.idxSteadyEnd;

% Force vectors into cell arrays
pulseVecsCell = force_column_cell(responseVecs);

% Vertically concatenate two sets
pulseVecsCellRepeated = repmat(pulseVecsCell, [2, 1]);

% Generate the indices for left and right triangle markers
indRightTriangles = [idxBaseStart; idxSteadyStart];
indLeftTriangles = [idxBaseEnd; idxSteadyEnd];

% Generate the x and y locations for markers
xRightTriangles = timeVec(indRightTriangles);
yRightTriangles = cellfun(@(x, y) x(y), ...
                            pulseVecsCellRepeated, num2cell(indRightTriangles));
xLeftTriangles = timeVec(indLeftTriangles);
yLeftTriangles = cellfun(@(x, y) x(y), ...
                            pulseVecsCellRepeated, num2cell(indLeftTriangles));

%% Do the job
% Hold on
hold on;

% Plot vectors
if nToUse == nVectors
    % Plot all vectors at once
    plot(timeVec, responseVecs, '-');
else
    % Plot vectors as dotted line if not used
    parfor iVec = 1:nVectors
        if toUse(iVec)
            plot(timeVec, responseVecs(:, iVec), '-');
        else
            plot(timeVec, responseVecs(:, iVec), '--');
        end
    end
end

% Plot detected endpoints
plot(xLeftTriangles, yLeftTriangles, '<');
plot(xRightTriangles, yRightTriangles, '>');

% Set time axis limits
if ~strcmpi(xLimits, 'suppress')
    xlim(xLimits);
end

% Set y axis limits
if ~strcmpi(yLimits, 'suppress')
    ylim(yLimits);
end

% Label the axes
xlabel('Time (ms)')
ylabel('Voltage (mV)')

% TODO
h = gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract all pulse amplitudes
pulseAmplitude = responseParams.pulseAmplitude;

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
    plot(timeVec(indBaseline(iVec, 1)), responseVecs(indBaseline(iVec, 1), iVec), '>');
    plot(timeVec(indBaseline(iVec, end)), responseVecs(indBaseline(iVec, end), iVec), '<');
    plot(timeVec(indSteady(iVec, 1)), responseVecs(indSteady(iVec, 1), iVec), '>');
    plot(timeVec(indSteady(iVec, end)), responseVecs(indSteady(iVec, end), iVec), '<');
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%