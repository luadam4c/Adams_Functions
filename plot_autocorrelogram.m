function handles = plot_autocorrelogram (autoCorrData, autoCorrParams, varargin)
%% Plots an autocorrelation function from the results of compute_autocorrelogram.m
% Usage: handles = plot_autocorrelogram (autoCorrData, autoCorrParams, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       plot_autocorrelogram(autoCorrData, autoCorrParams)
%       plot_autocorrelogram(autoCorrData, autoCorrParams, 'PlotType', 'acfFiltered')
%
% Outputs:
%       handles     - a structure with fields:
%                       TODO
%                   specified as a scalar structure
%
% Arguments:
%       autoCorrData    - autocorrelogram data
%                       must be a scalar structure with fields:
%                       TODO
%       autoCorrParams  - autocorrelogram parameters
%                       must be a scalar structure with fields:
%                       TODO
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'autocorrelogram' - two-sided autocorr
%                       'acfFiltered'     - one-sided autocorrelation function
%                   default == 'autocorrelogram'
%                   - 'PlotFiltered': whether to plot filtered acf
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotPeaks': whether to plot detected peaks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotTroughs': whether to plot detected troughs
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotDuration': whether to plot oscillation duration
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotText': whether to plot detected results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'XLimits': x-axis limits
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == minimum and maximum edges of bins
%                   - 'YLimits': limits of y axis, 
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == 'suppress'
%                   - 'BarYValue': TODO
%                   must be TODO
%                   default == -(yLimits(2) * 0.025)
%                   - Any other parameter-value pair for plot()
%
% Requires:
%       cd/compute_axis_limits.m
%       cd/compute_stats.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/hold_off.m
%       cd/hold_on.m
%       cd/plot_horizontal_line.m
%       cd/set_figure_properties.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_oscillations_analyze.m
%       cd/parse_multiunit.m

% File History:
% 2019-12-02 Moved from parse_multiunit.m
% 

%% Hard-coded parameters
validePlotTypes = {'autocorrelogram', 'acfFiltered'};
% TODO:
figHandle = [];

%% Default values for optional arguments
plotTypeDefault = 'autocorrelogram';
plotFilteredDefault = true;     % plot filtered traces by default
plotPeaksDefault = true;        % plot detected peaks by default
plotTroughsDefault = true;      % plot detected troughs by default
plotDurationDefault = true;     % plot duration by default
plotTextDefault = true;         % plot detected results by default
xLimitsDefault = [];                    % set later
yLimitsDefault = [];                    % set later
barYValueDefault = [];                  % set later

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
addRequired(iP, 'autoCorrData', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));
addRequired(iP, 'autoCorrParams', ...
    @(x) validateattributes(x, {'struct'}, {'scalar'}));


% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validePlotTypes)));
addParameter(iP, 'PlotFiltered', plotFilteredDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotPeaks', plotPeaksDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotTroughs', plotTroughsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotDuration', plotDurationDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotText', plotTextDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isvector(x) && length(x) == 2);
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'BarYValue', barYValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, autoCorrData, autoCorrParams, varargin{:});
plotType = validatestring(iP.Results.PlotType, validePlotTypes);
plotFiltered = iP.Results.PlotFiltered;
plotPeaks = iP.Results.PlotPeaks;
plotTroughs = iP.Results.PlotTroughs;
plotDuration = iP.Results.PlotDuration;
plotText = iP.Results.PlotText;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
barYValue = iP.Results.BarYValue;

% Keep unmatched arguments for the plot() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Retrieve data for plotting
autoCorr = autoCorrData.autoCorr;
acf = autoCorrData.acf;
acfFiltered = autoCorrData.acfFiltered;
indPeaks = autoCorrData.indPeaks;
indTroughs = autoCorrData.indTroughs;
ampPeaks = autoCorrData.ampPeaks;
ampTroughs = autoCorrData.ampTroughs;

binWidthSec = autoCorrParams.binWidthSec;
nBins = autoCorrParams.nBins;
halfNBins = autoCorrParams.halfNBins;
oscIndex1 = autoCorrParams.oscIndex1;
oscIndex2 = autoCorrParams.oscIndex2;
oscIndex3 = autoCorrParams.oscIndex3;
oscIndex4 = autoCorrParams.oscIndex4;
oscPeriod2Ms = autoCorrParams.oscPeriod2Ms;
oscPeriod1Ms = autoCorrParams.oscPeriod1Ms;
oscDurationSec = autoCorrParams.oscDurationSec;
nSpikesInOsc = autoCorrParams.nSpikesInOsc;

% Return if nothing to plot
if isempty(autoCorr)
    handles = struct;
    return;
end

% Decide on figure title base
if isfield(autoCorrParams, 'figTitleBase')
    figTitleBase = autoCorrParams.figTitleBase;
else
    figTitleBase = 'unnamed data';
end

% Create time values 
if nBins > 1
    tAcfTemp = create_time_vectors(nBins - 1, 'SamplingIntervalSec', binWidthSec, ...
                                    'TimeUnits', 's');
    tAcf = [0; tAcfTemp(1:halfNBins)];
    tAutoCorr = [-flipud(tAcfTemp); 0; tAcfTemp];
    timePeaksSec = (indPeaks - indPeaks(1)) * binWidthSec;
    timeTroughsSec = (indTroughs - indPeaks(1)) * binWidthSec;
else
    tAcf = NaN(size(acf));
    tAutoCorr = NaN(size(autoCorr));
    timePeaksSec = NaN(size(ampPeaks));
    timeTroughsSec = NaN(size(ampTroughs));
end

% Compute the x limits for the oscillation duration line
xLimitsOscDur = [0, oscDurationSec];

% Compute default x limits
if isempty(xLimits)
    allLastPeaksBins = extract_elements(indPeaks, 'last');
    allLastPeaksSec = allLastPeaksBins .* binWidthSec;
    allOscDur = oscDurationSec;
    bestRightForAll = max([allOscDur, allLastPeaksSec], [], 2) + 1;
    acfFilteredRight = compute_stats(bestRightForAll, 'upper95', ...
                                    'RemoveOutliers', true);

    switch plotType
    case 'autocorrelogram'
        xLimits = [-acfFilteredRight, acfFilteredRight];
    case 'acfFiltered'
        xLimits = [0, acfFilteredRight];
    end
end

% Compute default y limits
if isempty(yLimits)
    % Find the best upper limits
    lastIndexToShow = floor(xLimits(2) ./ binWidthSec) + 1;
    acfOfInterest = extract_subvectors(acf, 'IndexEnd', lastIndexToShow);
    largestAcfValues = extract_elements(acfOfInterest, 'max');
    bestUpperLimit = compute_stats(largestAcfValues, 'upper95', ...
                                    'RemoveOutliers', true);

    % Compute appropriate y limits
    switch plotType
    case 'autocorrelogram'
        yLimits = compute_axis_limits([0, bestUpperLimit], ...
                                        'y', 'Coverage', 95);
    case 'acfFiltered'
        yLimits = compute_axis_limits([0, bestUpperLimit], ...
                                        'y', 'Coverage', 90);
    end
end

% Compute default oscillation duration bar y value
if isempty(barYValue) && numel(yLimits) == 2
    barYValue = -(yLimits(2) * 0.025);
end

% Decide on the figure title
switch plotType
    case 'autocorrelogram'
        figTitle = ['Autocorrelation for ', figTitleBase];
    case 'acfFiltered'
        figTitle = ['Autocorrelation function for ', figTitleBase];
end

% Decide on the figure handle
fig = set_figure_properties('FigHandle', figHandle);

%% Plot
% Hold on
wasHold = hold_on;

% Plot stuff
switch plotType
    case 'autocorrelogram'
        % Just plot the autocorrelogram
        lines = plot(tAutoCorr, autoCorr, otherArguments{:});
    case 'acfFiltered'
        % Plot the autocorrelation function
        lines(1) = plot(tAcf, acf, 'k', otherArguments{:});

        % Plot filtered acf if requested
        if plotFiltered
            lines(2) = plot(tAcf, acfFiltered, 'r', 'LineWidth', 1);
        end

        % Plot detected peaks if requested
        if plotPeaks
            plot(timePeaksSec, ampPeaks, 'go', 'LineWidth', 2);
        end

        % Plot detected troughs if requested
        if plotTroughs
            plot(timeTroughsSec, ampTroughs, 'bx', 'LineWidth', 2);
        end

        % Plot oscillation duration if requested
        if plotDuration && ~isempty(barYValue)
            plot_horizontal_line(barYValue, 'XLimits', xLimitsOscDur, ...
                                'Color', 'g', 'LineStyle', '-', 'LineWidth', 2);
        end

        % Plot analyzed numbers as text if requested
        if plotText
            text(0.3, 0.98, sprintf('Oscillatory Index 4 = %.2g', ...
                oscIndex4), 'Units', 'normalized');
            text(0.3, 0.94, sprintf('Oscillatory Index 3 = %.2g', ...
                oscIndex3), 'Units', 'normalized');
            text(0.3, 0.90, sprintf('Oscillatory Index 2 = %.2g', ...
                oscIndex2), 'Units', 'normalized');
            text(0.3, 0.86, sprintf('Oscillatory Index 1 = %.2g', ...
                oscIndex1), 'Units', 'normalized');
            text(0.3, 0.82, sprintf('Oscillation Period 2 = %.3g ms', ...
                oscPeriod2Ms), 'Units', 'normalized');
            text(0.3, 0.78, sprintf('Oscillation Period 1 = %.3g ms', ...
                oscPeriod1Ms), 'Units', 'normalized');
            text(0.3, 0.74, sprintf('Total spike count = %g', ...
                nSpikesInOsc), 'Units', 'normalized');
            text(0.3, 0.70, sprintf('Oscillation Duration = %.2g seconds', ...
                oscDurationSec), 'Units', 'normalized');
        end
    otherwise
        error('plotType unrecognized!');
end

% Set x axis limits
xlim(xLimits);

% Set y axis limits
if numel(yLimits) >= 2 && yLimits(1) < yLimits(2)
    ylim(yLimits);
end

% Set x axis label
xlabel('Lag (s)');

% Set y axis label
ylabel('Spike rate squared (Hz^2)');

% Set figure title
title(figTitle);

% Hold off
hold_off(wasHold);

%% Save handles in a structure
handles.fig = fig;
handles.lines = lines;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%