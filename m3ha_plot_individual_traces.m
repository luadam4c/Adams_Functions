function handles = m3ha_plot_individual_traces (tVecs, data, varargin)
%% Plots individual voltage traces
% Usage: handles = m3ha_plot_individual_traces (tVecs, data, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       handles     - a handles structure with fields:
%                       endPointsToPlot
%                       baseWindow
%                       fitWindow
%                       baseNoise
%                       sweepWeights
%                       fig
%                       subPlots
%                       plotsData
%                       plotsDataToCompare
%                       subTitles
%                       boundaries
%                   specified as a scalar structure
% Arguments:
%       tVecs       - time vector(s) for plotting
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       data        - data vectors(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'DataToCompare': data vector(s) to compare against
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == []
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(tVec), max(tVec)]
%                   - 'XUnits': x-axis units
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == 'unit'
%                   - 'XLabel': label for the time axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector 
%                       or a cell array of strings or character vectors
%                   default == ['Time (', xUnits, ')']
%                   - 'YLabel': label(s) for the y axis, 
%                               suppress by setting value to 'suppress'
%                   must be a string scalar or a character vector
%                   default == 'Data' if plotMode is 'overlapped'
%                               {'Trace #1', 'Trace #2', ...}
%                                   if plotMode is 'parallel'
%                   - 'ColorMap': a color map that also groups traces
%                                   each set of traces will be on the same row
%                                   if plot mode is 'parallel'
%                   must be a numeric array with 3 columns
%                   default == colormap(jet(nTraces))
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigHandle': figure handle for created figure
%                   must be a empty or a figure object handle
%                   default == []
%                   - 'FigNumber': figure number for creating figure
%                   must be empty or a positive integer scalar
%                   default == []
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == first half of the trace
%                   - 'FitWindow': time window to fit for each trace
%                   must be a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == second half of the trace
%                   - 'BaseNoise': baseline noise value(s)
%                   must be a numeric vector
%                   default == apply compute_default_sweep_info.m
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 1 ./ baseNoise
%                   - 'BaseErrors': sweep errors in baseline window
%                   must be a numeric vector
%                   default == apply compute_sweep_errors.m
%                   - 'SweepErrors': sweep errors in fitting window
%                   must be a numeric vector
%                   default == apply compute_sweep_errors.m
%                   - 'PlotSwpWeightsFlag': whether to plot sweep weights
%                   must be numeric/logical 1 (true) or 0 (false) or 'auto'
%                   default == 'auto'
%                   - Any other parameter-value pair for the plot() function
%
% Requires:
%       ~/Downloaded_Functions.m/rgb.m
%       cd/argfun.m
%       cd/compute_baseline_noise.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sweep_errors.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_column_numeric.m
%       cd/isbinaryscalar.m
%       cd/iscellnumeric.m
%       cd/isemptycell.m
%       cd/isnumericvector.m
%       cd/ispositiveintegerscalar.m
%       cd/match_format_vector_sets.m
%       cd/match_row_count.m
%       cd/plot_traces.m
%       cd/plot_window_boundaries.m
%       cd/save_all_figtypes.m
%
% Used by:    
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-10-29 Created by Adam Lu
% 2018-12-19 Now uses unmatched varargin parts as parameters for plot()
% 2018-12-19 Now returns all object handles in a structure
% 2018-12-19 Now does not create new figure if figNumber is not provided
% 2018-12-19 Now restricts vectors to x limits first
% 2018-12-20 Added 'XUnits', 'XLabel', 'YLabel' as parameters
% 2018-12-20 Now returns subTitles in handles
% 2018-12-20 Now computes baseErrors and sweepErrors
% 

%% Hard-coded parameters
maxNTracesForAnnotations = 8;
nSigFig = 3;
fontSize = 8;
plotMode = 'parallel';
% linkAxesOption = 'xy';
linkAxesOption = 'x';

%% Default values for optional arguments
dataToCompareDefault = [];      % no data to compare against by default
xLimitsDefault = [];            % set later
xUnitsDefault = 'ms';           % time in ms by default
xLabelDefault = '';             % set later
yLabelDefault = 'suppress';     % no y axis labels by default
colorMapDefault = [];           % set later
figTitleDefault = '';           % set later
figHandleDefault = [];          % no existing figure by default
figNumberDefault = [];          % no figure number by default
figNameDefault = '';            % don't save figure by default
baseWindowDefault = [];         % set later
fitWindowDefault = [];          % set later
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set later
baseErrorsDefault = [];         % set later
sweepErrorsDefault = [];        % set later
plotSwpWeightsFlagDefault = 'auto'; % set later

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'tVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'data', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataToCompare', dataToCompareDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'XUnits', xUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLabel', xLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'YLabel', yLabelDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault, ...
    @(x) isempty(x) || isnumeric(x) && size(x, 2) == 3);
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigHandle', figHandleDefault);
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'FigNumber must be a empty or a positive integer scalar!'));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['BaseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'FitWindow', fitWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['FitWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'BaseNoise', baseNoiseDefault, ...
    @(x) assert(isnumericvector(x), 'BaseNoise must be a numeric vector!'));
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeights must be a numeric vector!'));
addParameter(iP, 'BaseErrors', baseErrorsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepErrors must be a numeric vector!'));
addParameter(iP, 'SweepErrors', sweepErrorsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepErrors must be a numeric vector!'));
addParameter(iP, 'PlotSwpWeightsFlag', plotSwpWeightsFlagDefault, ...
    @(x) assert(isbinaryscalar(x) || ischar(x) && strcmpi(x, 'auto'), ...
                'PlotSwpWeightsFlag must be a binary scalar or ''auto''!'));

% Read from the Input Parser
parse(iP, tVecs, data, varargin{:});
dataToCompare = iP.Results.DataToCompare;
xLimits = iP.Results.XLimits;
xUnits = iP.Results.XUnits;
xLabel = iP.Results.XLabel;
yLabel = iP.Results.YLabel;
colorMap = iP.Results.ColorMap;
figTitle = iP.Results.FigTitle;
figHandle = iP.Results.FigHandle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
baseWindow = iP.Results.BaseWindow;
fitWindow = iP.Results.FitWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;
baseErrors = iP.Results.BaseErrors;
sweepErrors = iP.Results.SweepErrors;
plotSwpWeightsFlag = iP.Results.PlotSwpWeightsFlag;

% Keep unmatched arguments for the plot() function
otherArguments = iP.Unmatched;

%% Preparation
% Initialize output structure
handles = struct;

% If data is empty, return
if isempty(data) || iscell(data) && all(isemptycell(data))
    fprintf('Nothing to plot!\n');
    return
end

% Force time and data vectors as column cell arrays of column vectors
[tVecs, data] = argfun(@force_column_cell, tVecs, data);

% Count the number of sweeps
nSweeps = numel(data);

% Decide whether to annotate at all
if nSweeps <= maxNTracesForAnnotations
    toAnnotate = true;
else
    toAnnotate = false;
end

% Decide whether to plot sweep weights
if plotSwpWeightsFlag == 'auto'
    if nSweeps > 1 && toAnnotate
        plotSwpWeightsFlag = true;
    else
        plotSwpWeightsFlag = false;
    end
end

% Decide on the data for baseline noise and sweep weights
if isempty(dataToCompare)
    % Compute default windows, noise and weights from data
    dataForWeights = data;
else
    % Compute default windows, noise and weights from dataToCompare
    dataForWeights = dataToCompare;
end

% Compute default windows, noise and weights from dataForWeights
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tVecs, dataForWeights, ...
            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
            'BaseNoise', baseNoise, 'SweepWeights', sweepWeights);

% TODO: Merge below into a function and use it in m3ha_xolotl_plot.m as well
% Compute baseline errors if not provided
if isempty(baseErrors)
    if isempty(dataToCompare)
        % Just use the baseline noise
        baseErrors = baseNoise;
    else
        % Compute sweep errors over the baseline window
        errorStructTemp = compute_sweep_errors(data, dataToCompare, ...
                            'TimeVecs', tVecs, 'FitWindow', baseWindow, ...
                            'NormalizeError', false);

        % Extract baseline errors for each trace
        baseErrors = errorStructTemp.swpErrors;
    end
end

% Compute sweep errors if not provided
if isempty(sweepErrors)
    if isempty(dataToCompare)
        % Compute the baseline noise over the fitting window
        sweepErrors = compute_baseline_noise(dataForWeights, tVecs, fitWindow);
    else
        % Compute sweep errors over the fitting window
        errorStructTemp = compute_sweep_errors(data, dataToCompare, ...
                            'TimeVecs', tVecs, 'FitWindow', fitWindow, ...
                            'NormalizeError', false);

        % Extract sweep errors for each trace
        sweepErrors = errorStructTemp.swpErrors;
    end
end

% Match vectors format and numbers of sweep-dependent vectors with data
[tVecs, dataToCompare, baseWindow, fitWindow] = ...
    argfun(@(x) match_format_vector_sets(x, data), ...
            tVecs, dataToCompare, baseWindow, fitWindow);

% Make sure vectors are columns
[baseNoise, sweepErrors] = ...
    argfun(@force_column_numeric, baseNoise, sweepErrors);

% Match numbers of sweep-dependent scalars with data
[baseNoise, sweepErrors] = ...
    argfun(@(x) match_row_count(x, nSweeps), baseNoise, sweepErrors);

% Determine the number of rows and the number of traces per row
nRows = size(colorMap, 1);
nTracesPerRow = ceil(nSweeps / nRows);

% Initialize graphics object arrays for boundaries
subTitles = gobjects(nSweeps, 1);
boundaries = gobjects(nSweeps, 2);

%% Do the job
% Restrict to x limits for faster processing
%   Note: this should only be done after the errors are computed
if ~isempty(xLimits) && isnumeric(xLimits)
    % Find the end points to plot
    endPointsToPlot = find_window_endpoints(xLimits, tVecs);

    % Restrict to these end points
    [tVecs, data, dataToCompare] = ...
        argfun(@(x) extract_subvectors(x, 'EndPoints', endPointsToPlot), ...
                tVecs, data, dataToCompare);
else
    % Use the first and last indices
    endPointsToPlot = find_window_endpoints([], tVec);
end

% Plot traces
[fig, subPlots, plotsData, plotsDataToCompare] = ...
    plot_traces(tVecs, data, 'DataToCompare', dataToCompare, ...
                'ColorMap', colorMap, 'XLimits', xLimits, ...
                'XUnits', xUnits, 'XLabel', xLabel, 'YLabel', yLabel, ...
                'LegendLocation', 'suppress', ...
                'PlotMode', plotMode, 'LinkAxesOption', linkAxesOption, ...
                'FigHandle', figHandle, 'FigNumber', figNumber, ...
                otherArguments);

% Show sweep info and error as a subplot title 
%   TODO: make a function plot_title(subPlots, titles)
if toAnnotate
    for iSwp = 1:nSweeps
        % Generate the error string
        errorString = ['Noise = ', num2str(baseErrors(iSwp), nSigFig), '; ', ...
                        'RMSE = ', num2str(sweepErrors(iSwp), nSigFig)];

        % Get the subplot of interest
        subplot(subPlots(iSwp));

        % Show the error string as the subplot title
        subTitles(iSwp) = title(errorString);
    end
end

% Plot sweep weights
if plotSwpWeightsFlag
    for iSwp = 1:nSweeps
        % Get the subplot of interest
        subplot(subPlots(iSwp));

        % Hold on
        hold on

        % Get the current sweep weight
        sweepWeight = sweepWeights(iSwp);

        % Decide on the text color
        if sweepWeight ~= 0
            colorText = rgb('DarkGreen');
        else
            colorText = rgb('Gray');
        end

        % Show sweep weight
        text('String', ['w: ', num2str(sweepWeight, nSigFig)], ...
            'Color', colorText, 'FontSize', fontSize, ...
            'Position', [0.1 0.9], 'Units', 'normalized');
    end
end

% Plot fit windows as vertical lines
if toAnnotate
    for iSwp = 1:nSweeps
        % Get the subplot of interest
        subplot(subPlots(iSwp));

        % Hold on
        hold on

        % Plot fit window
        boundaries(iSwp, :) = ...
            plot_window_boundaries(fitWindow{iSwp}, ...
                                    'Color', 'g', 'LineStyle', '--');
    end
end

% Change settings for each subplot
%   Note: must do this after plotting vertical lines
for iSwp = 1:nSweeps
    % Get the subplot of interest
    subplot(subPlots(iSwp));

    % Make y axis automatically update
    axis 'auto y'
end

% Extract y axis limits
yLimits = zeros(nSweeps, 2);
for iSwp = 1:nSweeps
    yLimits(iSwp, :) = get(subPlots(iSwp), 'YLim');
end

% Create a title
suptitle(figTitle);

%% Output results
% Save figure
if ~isempty(figName)
    save_all_figtypes(fig, figName);
end

% Store in handles structure
handles.endPointsToPlot = endPointsToPlot;
handles.baseWindow = baseWindow;
handles.fitWindow = fitWindow;
handles.baseNoise = baseNoise;
handles.sweepWeights = sweepWeights;
handles.fig = fig;
handles.subPlots = subPlots;
handles.plotsData = plotsData;
handles.plotsDataToCompare = plotsDataToCompare;
handles.subTitles = subTitles;
handles.boundaries = boundaries;
handles.yLimits = yLimits;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if plotSwpWeightsFlag && nSweeps < 20

% Hold off
hold off

if sweepWeight ~= 0
    colorText = 'green';
else
    colorText = 'gray';
end

text('String', ['\color{', colorText, '} \bf ', ...
                num2str(sweepWeight, 2)], ...
        'Units', 'normalized', 'Position', [0.1 0.9]);

subplot(nRows, nTracesPerRow, iSwp); hold on;

'FigTitle', figTitle, 

plot_window_boundaries(fitWindow{iSwp}, ...
                        'LineColor', 'g', 'LineStyle', '--');

fig = figure('Visible', 'off');

fig = [];

figNumberDefault = 104;         % figure 104 by default

% Create and clear figure
if ~isempty(figNumber)
    fig = figure(figNumber);
else
    fig = gcf;
end
set(fig, 'Name', 'All individual voltage traces');
clf(fig);

% Extract y axis limits and store in individual structure
yLimits = zeros(nSweeps, 2);
for iSwp = 1:nSweeps
    % Use the limits from the left boundary line
    yLimits(iTrace, :) = boundaries(iSwp, 1).YData;
end
individual.yLimits = yLimits;

% Compute sweep errors over the fitting window
errorStructTemp = compute_sweep_errors(data, dataToCompare, ...
                    'TimeVecs', tVecs, 'FitWindow', fitWindow, ...
                    'SweepWeights', sweepWeights, 'NormalizeError', false);

% Extract sweep errors for each trace
sweepErrors = errorStructTemp.swpErrors;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
