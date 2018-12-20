function xolotlObject = m3ha_xolotl_plot (xolotlObject, varargin)
%% Plots the simulation results from a xolotl object
% Usage: xolotlObject = m3ha_xolotl_plot (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%       varargin    - 'DataToCompare': data vector(s) to compare against
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == []
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == [min(tVec), max(tVec)]
%                   - 'ColorMap': a color map that also groups traces
%                                   each set of traces will be on the same row
%                                   if plot mode is 'parallel'
%                   must be a numeric array with 3 columns
%                   default == colormap(jet(nTraces))
%                   - 'FigTitle': title for the figure
%                   must be a string scalar or a character vector
%                   default == ['Traces for ', figName]
%                               or [yLabel, ' over time']
%                   - 'FigNumber': figure number for creating figure
%                   must be empty or a positive integer scalar
%                   default == 104
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
%                   - 'SweepErrors': sweep errors
%                   must be a numeric vector
%                   default == apply compute_sweep_errors.m
%                   - 'PlotSwpWeightsFlag': whether to plot sweep weights
%                   must be numeric/logical 1 (true) or 0 (false) or 'auto'
%                   default == 'auto'
%
% Requires:
%       cd/argfun.m
%       cd/compute_axis_limits.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/create_time_vectors.m
%       cd/find_ind_str_in_cell.m
%       cd/m3ha_plot_individual_traces.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2018-12-19 Now uses Children to find the soma column
% 2018-12-19 Now updates plots by changing the data values
% 2018-12-19 Now restricts vectors to x limits first
% 2018-12-19 Now updates the y limits for the boundaries
% 2018-12-20 Added x and y labels
% 2018-12-20 Now updates error strings
% 

%% Hard-coded parameters
compsToLabel = {'dend2', 'dend1', 'soma'};
yLabels = {'Dendrite 2 (mV)', 'Dendrite 1 (mV)', 'Soma (mV)', 'Stim (nA)'};
nSigFig = 3;

%% Default values for optional arguments
dataToCompareDefault = [];      % no data to compare against by default
xLimitsDefault = [];            % for current pulse
colorMapDefault = [];           % set in m3ha_plot_individual_traces.m
figTitleDefault = 'Simulation by xolotl';
figNumberDefault = [];          % set in m3ha_plot_individual_traces.m
figNameDefault = '';            % don't save figure by default
baseWindowDefault = [];         % set later
fitWindowDefault = [];          % set later
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set later
baseErrorsDefault = [];         % set later
sweepErrorsDefault = [];        % set later
plotSwpWeightsFlagDefault = false;      % will not update yet if true

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
addRequired(iP, 'xolotlObject');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataToCompare', dataToCompareDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'ColorMap', colorMapDefault, ...
    @(x) isempty(x) || isnumeric(x) && size(x, 2) == 3);
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
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
parse(iP, xolotlObject, varargin{:});
dataToCompare = iP.Results.DataToCompare;
xLimits = iP.Results.XLimits;
colorMap = iP.Results.ColorMap;
figTitle = iP.Results.FigTitle;
figNumber = iP.Results.FigNumber;
figName = iP.Results.FigName;
baseWindow = iP.Results.BaseWindow;
fitWindow = iP.Results.FitWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;
baseErrors = iP.Results.BaseErrors;
sweepErrors = iP.Results.SweepErrors;
plotSwpWeightsFlag = iP.Results.PlotSwpWeightsFlag;

%% Preparation
% Extract the compartment names
compNames = xolotlObject.Children;

% Find the indices for compartments to label in xolotl
idxInXolotl = ...
    cellfun(@(x) find_ind_str_in_cell(x, compNames, 'IgnoreCase', true, ...
                                'SearchMode', 'substrings', 'MaxNum', 1), ...
            compsToLabel);

%% Initialize a plot or retrieve plot
% Retrieve handles from xolotl object
xHandles = xolotlObject.handles;

% TODO: Make this a function?
if isempty(xHandles) || ~isfield(xHandles, 'individual')
    % Extract the time step in ms
    timeStep = xolotlObject.dt;

    % Extract the external current injection protocol
    currentProtocol = xolotlObject.I_ext;

    % Find the idx for soma
    idxSoma = find_ind_str_in_cell('soma', compNames, 'IgnoreCase', true, ...
                                    'SearchMode', 'substrings', 'MaxNum', 1);

    % Save the stimulation protocol
    iStim = currentProtocol(:, idxSoma);

    % Get the number of samples for the stimulation protocol
    nSamples = count_samples(iStim);

    % Create a time vector in milliseconds
    tVec = create_time_vectors(nSamples, 'SamplingIntervalMs', timeStep, ...
                                'TimeUnits', 'ms');

    % Restrict to x limits for faster processing
    if ~isempty(xLimits) && isnumeric(xLimits)
        % Find the end points
        endPointsToPlot = find_window_endpoints(xLimits, tVec);

        % Restrict to these end points
        [tVec, iStim, dataToCompare] = ...
            argfun(@(x) extract_subvectors(x, 'EndPoints', endPointsToPlot), ...
                    tVec, iStim, dataToCompare);
    else
        % Use the first and last indices
        endPointsToPlot = find_window_endpoints([], tVec);
    end

    % Create NaN data for the initial plot
    nanVoltageData = NaN * tVec;
    nanData = [nanVoltageData, nanVoltageData, nanVoltageData, iStim];

    % Make a figure
    figHandle = ...
        figure('Outerposition', [100, 100, 800, 800], ...
            'PaperUnits', 'points', 'PaperSize', [1200, 600]); hold on

    % Initialize the plot by using NaN data
    individual = m3ha_plot_individual_traces(tVec, nanData, ...
                                    'DataToCompare', dataToCompare, ...
                                    'XLimits', xLimits, ...
                                    'XUnits', 'ms', ...
                                    'YLabel', yLabels, ...
                                    'ColorMap', colorMap, ...
                                    'FigTitle', figTitle, ...
                                    'FigHandle', figHandle, ...
                                    'FigNumber', figNumber, ...
                                    'FigName', figName, ...
                                    'BaseWindow', baseWindow, ...
                                    'FitWindow', fitWindow, ...
                                    'BaseNoise', baseNoise, ...
                                    'SweepWeights', sweepWeights, ...
                                    'SweepErrors', sweepErrors, ...
                                    'PlotSwpWeightsFlag', plotSwpWeightsFlag);

    % Extract y axis limits
    yLimits = zeros(3, 2);
    for iTrace = 1:4
        % Use the limits from the left boundary line
        yLimits(iTrace, :) = individual.boundaries(iTrace, 1).YData;
    end

    % If using x.manipulate
    if isfield(xHandles, 'puppeteer_object')
        % Attach figure to puppeteer so that 
        %   "the figure can be closed automatically"
        xHandles.puppeteer_object.attachFigure(individual.fig);
    end

    % Store information in individual structure
    individual.endPointsToPlot = endPointsToPlot;
    individual.baseWindow = baseWindow;
    individual.fitWindow = fitWindow;
    individual.yLimits = yLimits;

    % Store the figure handle in the xolotl object
    xHandles.individual = individual;
else
    % Extract the individual structure from the xolotl object
    individual = xHandles.individual;

    % Extract information from the individual structure
    endPointsToPlot = individual.endPointsToPlot;
    baseWindow = individual.baseWindow;
    fitWindow = individual.fitWindow;
    yLimits = individual.yLimits;

    % Extract tVec from the first plot
    tVec = individual.plotsData(1).XData;

    % Recompute the number of samples
    nSamples = count_samples(tVec);

    % Extract dataToCompare
    dataToCompare = zeros(nSamples, 3);
    for iTrace = 1:3    
        dataToCompare(:, iTrace) = individual.plotsDataToCompare(iTrace).YData;
    end
end

%% Simulate
% Get voltage traces for all compartments
vVecs = xolotlObject.integrate;

% Restrict to same end points to be consistent with the default plot
vVecs = extract_subvectors(vVecs, 'EndPoints', endPointsToPlot);

%% Computed updated data for plots
% Reorganize voltage traces to match plots
vVecPlots = vVecs(:, idxInXolotl);

% Compute new y limits
newYLimits = zeros(3, 2);
parfor iTrace = 1:3
    % Put the new data and old y limits together
    newDataForLimits = {vVecPlots(:, iTrace); yLimits(iTrace, :)};
    
    % Update y limits
    newYLimits(iTrace, :) = compute_axis_limits(newDataForLimits, 'y');
end

% Re-compute default windows, noise and weights
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tVec, vVecPlots, ...
            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
            'BaseNoise', baseNoise, 'SweepWeights', sweepWeights);

% Re-compute baseline errors if not provided
if isempty(baseErrors)
    if isempty(dataToCompare)
        % Just use the baseline noise
        baseErrors = baseNoise;
    else
        % Compute sweep errors over the baseline window
        errorStructTemp = compute_sweep_errors(vVecPlots, dataToCompare, ...
                            'TimeVecs', tVec, 'FitWindow', baseWindow, ...
                            'SweepWeights', sweepWeights, 'NormalizeError', false);

        % Extract baseline errors for each trace
        baseErrors = errorStructTemp.swpErrors;
    end
end

% compute sweep errors if not provided
if isempty(sweepErrors)
    if isempty(dataToCompare)
        % Compute the baseline noise over the fitting window
        sweepErrors = compute_baseline_noise(dataForWeights, tVecs, fitWindow);
    else
        % Compute sweep errors over the fitting window
        errorStructTemp = compute_sweep_errors(vVecPlots, dataToCompare, ...
                            'TimeVecs', tVec, 'FitWindow', fitWindow, ...
                            'SweepWeights', sweepWeights, 'NormalizeError', false);

        % Extract sweep errors for each trace
        sweepErrors = errorStructTemp.swpErrors;
    end
end

% Re-generate the error strings
errorStrings = ...
    arrayfun(@(x) ['Noise = ', num2str(baseErrors(x), nSigFig), '; ', ...
                    'RMSE = ', num2str(sweepErrors(x), nSigFig)], 1:3, ...
                'UniformOutput', false);

%% Update plots
% Update data in the chart line objects
for iTrace = 1:3    
    individual.plotsData(iTrace).YData = vVecPlots(:, iTrace);
end

% Update data in the primitive line objects
for iTrace = 1:3
    individual.boundaries(iTrace, 1).YData = newYLimits(iTrace, :);
    individual.boundaries(iTrace, 2).YData = newYLimits(iTrace, :);
end

% Update subplot titles
for iTrace = 1:3
    individual.subTitles(iTrace).String = errorStrings{iTrace};
end


%% Save handles in xolotl object
xolotlObject.handles = xHandles;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

tVec = (1:length(vVecs)) * xolotlObject.dt * 1e-3;

set(gca, 'YLim', [-80 50])

nSamples = count_samples(vVecs);

% Create figure
figure('Outerposition', [300, 300, 1200, 600], ...
        'PaperUnits', 'points', 'PaperSize', [1200, 600]);
hold on

% Create first subplot
subplot(4, 1, 1); hold on
plot(tVec, vVecDendrite2, 'r')
ylabel('Voltage in dendrite 2 (mV)')
set(gca, 'XLim', xLimits)

% Create second subplot
subplot(4, 1, 2); hold on
plot(tVec, vVecDendrite1, 'r')
ylabel('Voltage in dendrite 1 (mV)')
set(gca, 'XLim', xLimits)

% Create third subplot
subplot(4, 1, 3); hold on
plot(tVec, vVecSoma, 'r')
ylabel('Voltage in soma (mV)')
set(gca, 'XLim', xLimits)

% Create fourth subplot
subplot(4, 1, 4); hold on
plot(tVec, iVec, 'r')
xlabel('Time (s)')
ylabel('Stimulation Pulse (nA)')
set(gca, 'XLim', xLimits)

xolotlObject.plot;

subplot(3, 1, 1); hold on
xlim(xLimits)

subplot(3, 1, 2); hold on
xlim(xLimits)

subplot(3, 1, 3); hold on
xlim(xLimits)

'FigNumber', figNumber, ...
'FigName', figName, ...

% Plot the traces
    m3ha_plot_individual_traces(tVec, data, ...
                                'DataToCompare', dataToCompare, ...
                                'XLimits', xLimits, ...
                                'ColorMap', colorMap, ...
                                'FigTitle', figTitle, ...
                                'FigNumber', figNumber, ...
                                'FigName', figName, ...
                                'BaseWindow', baseWindow, ...
                                'FitWindow', fitWindow, ...
                                'BaseNoise', baseNoise, ...
                                'SweepWeights', sweepWeights, ...
                                'SweepErrors', sweepErrors, ...
                                'PlotSwpWeightsFlag', plotSwpWeightsFlag);

% Retrieve all parameter values in the object tree
paramValues = xolotlObject.serialize;

% Restore the xolotl object to its original state
xolotlObject.deserialize(paramValues);

% Place all traces into a data array
data = [vVecDendrite2, vVecDendrite1, vVecSoma, iStim];

% Update figure
drawnow;

xolotlObject.handles.puppeteer_object.attachFigure(individual.fig);
xolotlObject.handles.individual = individual;
xolotlObject.handles.individual.plotsData(1).YData = vVecDendrite2;
xolotlObject.handles.individual.plotsData(2).YData = vVecDendrite1;
xolotlObject.handles.individual.plotsData(3).YData = vVecSoma;

% Extract the voltage traces for each compartment
vVecDendrite2 = vVecs(:, idxDend2);
vVecDendrite1 = vVecs(:, idxDend1);
vVecSoma = vVecs(:, idxSoma);

individual.plotsData(1).YData = vVecDendrite2;
individual.plotsData(2).YData = vVecDendrite1;
individual.plotsData(3).YData = vVecSoma;

% Extract the old minimums and maximums
[oldMinY, oldMaxY] = argfun(@(x) oldYLimits(iTrace, x), 1, 2);

% Compute the minimum of new data and old y limits
newMinY = apply_iteratively(@min, {vVecPlots(:, iTrace); oldMinY});
newMaxY = apply_iteratively(@max, {vVecPlots(:, iTrace); oldMaxY});

% Update y limits
newYLimits(iTrace, :) = compute_ylimits(newMinY, newMaxY, 'Coverage', 80);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%