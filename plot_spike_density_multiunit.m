function handles = plot_spike_density_multiunit (parsedData, parsedParams, varargin)
%% Plots a spike density plot from parsed multiunit data
% Usage: handles = plot_spike_density_multiunit (parsedData, parsedParams, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - a structure with fields:
%                       im
%                       boundaryLine
%                       stimStartLine
%                   specified as a scalar structure
%
% Arguments:
%       parsedData  - parsed data
%                   must be a table with fields:
%                       spikeDensityHz
%       parsedParams- parsed parameters
%                   must be a table with fields:
%                       siSeconds
%                       minTimeSec
%                       maxTimeSec
%                       stimStartSec
%                       phaseNumber
%                       figPathBase
%       varargin    - 'FigPosition':
%                   must be a TODO
%                   default == TODO
%                   - 'XLimits':
%                   must be a TODO
%                   default == TODO
%                   - 'PlotStimStart':
%                   must be a TODO
%                   default == TODO
%                   - 'BoundaryType':
%                   must be a TODO
%                   default == TODO
%                   - 'MaxNYTicks':
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for imagesc()
%
% Requires:
%       cd/compute_index_boundaries.m
%       cd/create_indices.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/decide_on_colormap.m
%       cd/extract_common_prefix.m
%       cd/force_matrix.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/plot_vertical_line.m
%       cd/plot_window_boundaries.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-11-30 Moved from parse_multiunit.m
% TODO: Pull out code as function plot_heat_map(spikeDensityHz);

%% Hard-coded parameters
validBoundaryTypes = {'horizontalLines', 'verticalBars'};

% TODO: Make optional argument
colorMap = [];

% Default values for optional arguments
figPositionDefault = [];
xLimitsDefault = [];
plotStimStartDefault = true;
boundaryTypeDefault = 'horizontalLines';
maxNYTicksDefault = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to an input Parser
addRequired(iP, 'parsedData', ..
    @(x) validateattributes(x, {'table'}, {'2d'}));
addRequired(iP, 'parsedParams', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigPosition', figPositionDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PlotStimStart', plotStimStartDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'BoundaryType', boundaryTypeDefault, ...
    @(x) any(validatestring(x, validBoundaryTypes)));
addParameter(iP, 'MaxNYTicks', maxNYTicksDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, parsedData, parsedParams, varargin{:});
figPosition = iP.Results.FigPosition;
xLimits = iP.Results.XLimits;
plotStimStart = iP.Results.PlotStimStart;
boundaryType = validatestring(iP.Results.BoundaryType, validBoundaryTypes);
maxNYTicks = iP.Results.MaxNYTicks;

% Keep unmatched arguments for the imagesc() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Retrieve data for plotting
spikeDensityHz = parsedData.spikeDensityHz;

siSeconds = parsedParams.siSeconds;
minTimeSec = parsedParams.minTimeSec;
maxTimeSec = parsedParams.maxTimeSec;
stimStartSec = parsedParams.stimStartSec;
phaseVector = parsedParams.phaseNumber;
figPathBase = parsedParams.figPathBase;

% Extract the original file base
fileBase = extract_common_prefix(figPathBase);

% Create a figure title base
titleBase = replace(fileBase, '_', '\_');

% Compute the phase boundaries
phaseBoundaries = compute_index_boundaries('Grouping', phaseVector, ...
                                            'TreatNaNsAsGroup', false);

%% Plot as a heatmap
% Hold on
wasHold = hold_on;

% Compute stimulation start
meanStimStartSec = mean(stimStartSec);

% Count traces
nSweeps = numel(spikeDensityHz);

% Get the average sampling interval in seconds
siSeconds = mean(siSeconds);

% Set x and y end points
xEnds = [min(minTimeSec); max(maxTimeSec)];
yEnds = [1; nSweeps];

% Set x and y limits
if isempty(xLimits)
    xLimits = [xEnds(1) - 0.5 * siSeconds; xEnds(2) + 0.5 * siSeconds];
end
yLimits = [yEnds(1) - 0.5; yEnds(2) + 0.5];

% Compute bar X position if it is used
barXPosition = mean([xLimits(1), meanStimStartSec]);

% Decide on y ticks and labels
yTicks = create_indices('IndexEnd', nSweeps, 'MaxNum', maxNYTicks, ...
                        'AlignMethod', 'left');
yTickLabels = create_labels_from_numbers(nSweeps - yTicks + 1);

% Force as a matrix and transpose it so that
%   each trace is a row
spikeDensityMatrix = transpose(force_matrix(spikeDensityHz));

% Set a gray-scale color map
% colormap(flipud(gray));
% colormap(jet);
cm = decide_on_colormap(colorMap, 'ColorMapFunc', @gray, ...
                        'ReverseOrder', true, 'HighContrast', true);
colormap(cm);

% Generate plot
im = imagesc(xEnds, flipud(yEnds), spikeDensityMatrix, otherArguments{:});

% Set the y ticks and labels
yticks(yTicks);
yticklabels(yTickLabels);

% Plot stimulation start
if plotStimStart
    stimStartLine = plot_vertical_line(meanStimStartSec, 'Color', rgb('Green'), ...
                                    'LineStyle', '--', 'LineWidth', 0.5, ...
                                    'YLimits', yLimits);
end

% Plot phase boundaries
if ~isempty(phaseBoundaries)
    % Compute y values for phase boundaries
    yBoundaries = nSweeps - phaseBoundaries + 1;

    % Plot phase boundaries
    boundaryLine = ...
        plot_window_boundaries(yBoundaries, 'BoundaryType', boundaryType, ...
                                'BarValue', barXPosition);
end

% Update x limits and labels
xlim(xLimits);
ylim(yLimits);
xlabel('Time (s)');
ylabel('Trace #');
title(['Spike density (Hz) for ', titleBase]);

% Show a scale bar
colorbar;

% Set figure position if requested
if ~isempty(figPosition)
    set(gcf, 'Position', figPosition);
end

% Hold off
hold_off(wasHold);

%% Save output
handles.im = im;
if exist('boundaryLine', 'var')
    handles.boundaryLine = boundaryLine;
end
if exist('stimStartLine', 'var')
    handles.stimStartLine = stimStartLine;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%