function handles = plot_spike_histogram (spHistData, spHistParams, varargin)
%% Plots a spike histogram from the results of compute_spike_histogram.m
% Usage: handles = plot_spike_histogram (spHistData, spHistParams, varargin)
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
%       spHistData      - spike histogram data
%                       must be a scalar structure with fields:
%                           spikeCounts
%                           edgesSec
%                           timeBurstStartsSec
%                           timeBurstEndsSec
%       spHistParams    - spike histogram parameters
%                       must be a scalar structure with fields:
%                           oscDurationSec
%                           nSpikesInOsc
%                           figTitleBase
%       varargin    - 'XLimits': x-axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == minimum and maximum edges of bins
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - Any other parameter-value pair for plot_histogram()
%
% Requires:
%       cd/alternate_elements.m
%       cd/apply_iteratively.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/plot_histogram.m
%       cd/plot_horizontal_line.m
%
% Used by:
%       cd/m3ha_oscillations_analyze.m
%       cd/parse_multiunit.m

% File History:
% 2019-12-02 Moved from parse_multiunit.m
% 

%% Hard-coded parameters
minNXTicks = 5;

%% Default values for optional arguments
xLimitsDefault = [];            % set later
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
addRequired(iP, 'spHistData', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addRequired(iP, 'spHistParams', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);

% Read from the Input Parser
parse(iP, spHistData, spHistParams, varargin{:});
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;

% Keep unmatched arguments for the plot_histogram() function
otherArguments = iP.Unmatched;

%% Preparation
% Extract from spHistData
spikeCounts = spHistData.spikeCounts;
edgesSec = spHistData.edgesSec;
timeBurstStartsSec = spHistData.timeBurstStartsSec;
timeBurstEndsSec = spHistData.timeBurstEndsSec;

% Extract from spHistParams
oscDurationSec = spHistParams.oscDurationSec;
nSpikesInOsc = spHistParams.nSpikesInOsc;
figTitleBase = spHistParams.figTitleBase;

% Compute burst windows
burstWindows = alternate_elements(timeBurstStartsSec, timeBurstEndsSec, ...
                                    'ReturnNaNIfEmpty', true);

% Compute default x limits
if isempty(xLimits)
    % Extract parameters
    histLeftSec = spHistParams.histLeftSec;
    timeOscEndSec = spHistParams.timeOscEndSec;
    maxInterBurstIntervalSec = spHistParams.maxInterBurstIntervalSec;

    % Compute left and right ends of histogram to show
    histLeft = min(histLeftSec);
    histRight = timeOscEndSec + 1.5 * max(maxInterBurstIntervalSec);
    % histRight = 10;

    % Find appropriate x limits
    xLimits = [histLeft, histRight];
end

% Compute x tick locations
xTickLocs = linspace(xLimits(1), xLimits(2), minNXTicks);

% Compute default y limits
if isempty(yLimits)
    % Extract parameters
    binWidthSec = spHistParams.binWidthSec;

    % Find the last bin to show for all traces
    lastBinToShow = floor(range(xLimits) ./ binWidthSec) + 1;

    % Find appropriate y limits
    spikeCountsOfInterest = extract_subvectors(spikeCounts, ...
                            'IndexEnd', lastBinToShow);
    largestSpikeCount = apply_iteratively(@max, spikeCountsOfInterest);
    yLimits = [0, largestSpikeCount * 1.2];
end

%% Plot histogram
% Hold on
wasHold = hold_on;

% Plot histogram
[histBars, histFig] = ...
    plot_histogram([], 'Counts', spikeCounts, 'Edges', edgesSec, ...
                    'XLimits', xLimits, 'YLimits', yLimits, ...
                    'XTickLocs', xTickLocs, 'XLabel', 'Time (seconds)', ...
                    'YLabel', 'Spike Count per 10 ms', ...
                    'FigTitle', ['Spike histogram for ', figTitleBase], ...
                    otherArguments);

% Show oscillation duration
texts(1) = text(0.2, 0.95, sprintf('Oscillation Duration = %.2g seconds', ...
                oscDurationSec), 'Units', 'normalized');

% Show number of spikes in oscillation
texts(2) = text(0.2, 0.9, sprintf('Number of spikes in oscillation = %d', ...
                nSpikesInOsc), 'Units', 'normalized');

% Plot burst windows
horzLines = plot_horizontal_line(0, 'XLimits', burstWindows, ...
                            'Color', 'g', 'LineStyle', '-', 'LineWidth', 3);

% Hold off
hold_off(wasHold);

%% Output results
handles.fig = histFig;
handles.bars = histBars;
handles.texts = texts;
if exist('horzLines', 'var')
    handles.horzLines = horzLines;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%