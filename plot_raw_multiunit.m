function handles = plot_raw_multiunit (parsedData, parsedParams, varargin)
%% Plots the raw data from parsed multiunit data
% Usage: handles = plot_raw_multiunit (parsedData, parsedParams, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - a structure with fields:
%                       TODO
%                   specified as a scalar structure
%
% Arguments:
%       parsedData  - parsed data
%                   must be a table with fields:
%                       tVecs
%                       vVecs
%                       vVecsFilt
%       parsedParams- parsed parameters
%                   must be a table with fields:
%                       stimStartSec
%                       phaseVector
%                       figPathBase
%       varargin    - 'PlotFiltered': whether to plot filtered trace
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotStim': whether to plot stimulation time
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotMode': plotting mode for multiple traces
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'parallel'      - in parallel in subPlots
%                       'overlapped'    - overlapped in a single plot
%                       'staggered'     - staggered in a single plot
%                   default == 'staggered'
%                   - 'YAmountToStagger': amount to stagger 
%                                           if 'plotmode' is 'stagger'
%                   must be a positive scalar
%                   default == uses the original y axis range
%                   - 'SweepNumbers' - the sweep numbers to plot
%                   must be a numeric vector or 'all'
%                   default == 'all'
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == uses compute_axis_limits.m
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/compute_axis_limits.m
%       cd/compute_index_boundaries.m
%       cd/extract_common_prefix.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/plot_traces.m
%       cd/plot_vertical_line.m
%       cd/plot_horizontal_line.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-12-02 Moved from parse_multiunit.m
% 2019-12-02 Added 'SweepNumbers' as an optional argument
% 2019-12-02 Added 'YLimits' as an optional argument

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
validPlotModes = {'overlapped', 'parallel', 'staggered'};

%% Default values for optional arguments
plotFilteredDefault = false;    % don't plot filtered traces by default
plotStimDefault = true;         % plot stimulation time by default
plotModeDefault = 'staggered';  % plot traces staggered by default
yAmountToStaggerDefault = [];   % set later  
sweepNumbersDefault = 'all';    % plot all sweeps by default
yLimitsDefault = [];            % set later

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
addRequired(iP, 'parsedData', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'parsedParams', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotFiltered', plotFilteredDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotStim', plotStimDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'YAmountToStagger', yAmountToStaggerDefault, ...
    @(x) assert(isempty(x) || ispositivescalar(x), ...
                ['YAmountToStagger must be either a empty ', ...
                    'or a positive scalar!']));
addParameter(iP, 'SweepNumbers', sweepNumbersDefault);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, parsedData, parsedParams, varargin{:});
plotFiltered = iP.Results.PlotFiltered;
plotStim = iP.Results.PlotStim;
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
yAmountToStagger = iP.Results.YAmountToStagger;
sweepNumbers = iP.Results.SweepNumbers;
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% Restrict to specific sweeps if requested
if isnumeric(sweepNumbers) && ~isempty(sweepNumbers)
    [parsedData, parsedParams] = ...
        argfun(@(x) x(sweepNumbers, :), parsedData, parsedParams);
end

% Extract parameters
tVecs = parsedData.tVec;
vVecs = parsedData.vVec;
vVecsFilt = parsedData.vVecFilt;

% Extract parameters
stimStartSec = parsedParams.stimStartSec;
phaseVector = parsedParams.phaseNumber;
figPathBase = parsedParams.figPathBase;

% Compute the phase boundaries
phaseBoundaries = compute_index_boundaries('Grouping', phaseVector, ...
                                            'TreatNaNsAsGroup', false);

% Count the number of sweeps
nSweeps = height(parsedParams);

% Convert time vector to seconds
% TODO: Use convert_units.m instead
tVecsSec = transform_vectors(tVecs, MS_PER_S, 'divide');

% Decide on x axis label
xLabel = 'Time (s)';

% Decide on y axis label
switch plotMode
    case {'parallel', 'overlapped'}
        yLabel = 'Voltage (uV)';
    case 'staggered'
        yLabel = 'Trace #';
end

% Decide on title prefix
if plotFiltered
    titlePrefix = 'Raw and Filtered traces';
else
    titlePrefix = 'Raw traces';
end

% Decide on figure title(s)
switch plotMode
    case 'parallel'
        % Create figure titles
        titleBase = replace(figPathBase, '_', '\_');
        figTitle = 'suppress';
        figSubTitles = strcat(titlePrefix, {' for '}, titleBase);
    case  {'staggered', 'overlapped'}
        % Extract the original file base
        fileBase = extract_common_prefix(figPathBase);

        % Create figure titles
        titleBase = replace(fileBase, '_', '\_');
        figTitle = strcat(titlePrefix, {' for '}, titleBase);
        figSubTitles = {};
    otherwise
        error('plotMode unrecognized!');
end

% Compute the original y limits from data
if isempty(yLimits)
    yLimits = compute_axis_limits(vVecs, 'y', 'AutoZoom', true);
end

% Compute a default amount of y to stagger if not provided
if isempty(yAmountToStagger)
    yAmountToStagger = range(yLimits);
end

% Decide on what to plot
if plotFiltered
    dataToPlot = vVecsFilt;
    dataToCompare = vVecs;
else
    dataToPlot = vVecs;
    dataToCompare = [];
end

% Decide on the color map
if plotFiltered
    colorMap = 'b';
else
    colorMap = 'k';
end

%% Do the job
% Hold on
wasHold = hold_on;

% Plot traces
handles = plot_traces(tVecsSec, dataToPlot, 'DataToCompare', dataToCompare, ...
            'Verbose', false, 'PlotMode', plotMode, ...
            'YAmountToStagger', yAmountToStagger, ...
            'YLimits', yLimits, 'XLabel', xLabel, ...
            'YLabel', yLabel, 'TraceLabels', 'suppress', ...
            'FigTitle', figTitle, 'FigSubTitles', figSubTitles, ...
            'ColorMap', colorMap, 'LegendLocation', 'suppress', ...
            otherArguments);

% Plot stimulation start
if plotStim
    vertLine = plot_vertical_line(mean(stimStartSec), 'Color', 'g', ...
                                    'LineStyle', '--', 'LineWidth', 0.5);
end

% Plot phase boundaries
%   TODO: Use plot_window_boundaries.m instead
if ~isempty(phaseBoundaries) && strcmp(plotMode, 'staggered')
    yBoundaries = (nSweeps - phaseBoundaries + 1) * yAmountToStagger;
    horzLine = plot_horizontal_line(yBoundaries, 'Color', 'g', ...
                                    'LineStyle', '--', 'LineWidth', 2);
end

% Hold off
hold_off(wasHold);

%% Output results
if plotStim
    handles.vertLine = vertLine;
end
if exist('horzLine', 'var')
    handles.horzLine = horzLine;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%