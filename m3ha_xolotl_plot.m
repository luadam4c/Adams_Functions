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
%                   - 'HoldingPotential': holding potential
%                   must be a numeric vector
%                   default == baseline of dataToCompare
%                   - 'CompToPatch': compartment name to patch
%                   must be a string scalar or a character vector
%                   default == set in xolotl_compartment_index.m
%                   - 'TimeToStabilize': time to stabilize in ms
%                   must be a nonnegative scalar
%                   default == 1000 ms
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
%       cd/compute_baseline_noise.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/create_time_vectors.m
%       cd/find_in_strings.m
%       cd/find_window_endpoints.m
%       cd/force_row_vector.m
%       cd/m3ha_plot_individual_traces.
%       cd/parse_xolotl_object.m
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
% 2018-12-20 Now updates holding current here
% 2019-08-15 Fixed bugs for computing new y limits
% 

%% Hard-coded parameters
compsToLabel = {'dend2', 'dend1', 'soma'};
yLabels = {'Dendrite 2 (mV)', 'Dendrite 1 (mV)', 'Soma (mV)', 'Stim (nA)'};
nSigFig = 3;
xUnits = 'ms';

%% Default values for optional arguments
dataToCompareDefault = [];      % no data to compare against by default
compToPatchDefault = '';        % set later
timeToStabilizeDefault = 1000;  % simulate for 1000 ms by default
holdingPotentialDefault = [];   % set later
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
addParameter(iP, 'HoldingPotential', holdingPotentialDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'CompToPatch', compToPatchDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeToStabilize', timeToStabilizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
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
holdingPotential = iP.Results.HoldingPotential;
compToPatch = iP.Results.CompToPatch;
timeToStabilize = iP.Results.TimeToStabilize;
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
% Extract all compartments
compartments = xolotlObject.Children;

% Count the number of compartments
nCompartments = numel(compartments);

% Count the number of vectors
nTraces = nCompartments + 1;

% Find the indices for compartments to label in xolotl
idxInXolotl = ...
    cellfun(@(x) find_in_strings(x, compartments, 'IgnoreCase', true, ...
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

    % Find the idx for the compartment to patch
    idxCompToPatch = ...
        find_in_strings(compToPatch, compartments, 'MaxNum', 1, ...
                            'SearchMode', 'substrings', 'IgnoreCase', true);

    % Save the stimulation protocol
    iStim = currentProtocol(:, idxCompToPatch);

    % Get the number of samples for the stimulation protocol
    nSamples = count_samples(iStim);

    % Create a time vector in milliseconds
    tVec = create_time_vectors(nSamples, 'SamplingIntervalMs', timeStep, ...
                                'TimeUnits', 'ms');

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
                                    'XUnits', xUnits, ...
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

    % If using x.manipulate
    if isfield(xHandles, 'puppeteer_object')
        % Attach figure to puppeteer so that 
        %   "the figure can be closed automatically"
        xHandles.puppeteer_object.attachFigure(individual.fig);
    end

    % If not provided, compute the holding potential(s) needed for simulations
    if isempty(holdingPotential)
        % Find the baseline window endpoints
        baseEndPoints = find_window_endpoints(individual.baseWindow, tVec);
        
        % Compute the holding potential from dataToCompare
        holdingPotential = ...
            compute_stats(dataToCompare, 'mean', 'EndPoints', baseEndPoints, ...
                            'IgnoreNaN', true);
                        
        % Force as a row vector
        holdingPotential = force_row_vector(holdingPotential);
    end

    % Retrieve the original external current
    externalCurrentOrig = xolotlObject.I_ext;
    
    % Store arguments in handles structure 
    %   (so that they can be used by x.manipulate)
    xHandles.tVec = tVec;
    xHandles.dataToCompare = dataToCompare;
    xHandles.compToPatch = compToPatch;
    xHandles.idxCompToPatch = idxCompToPatch;
    xHandles.timeToStabilize = timeToStabilize;
    xHandles.holdingPotential = holdingPotential;
    xHandles.externalCurrentOrig = externalCurrentOrig;

    % Store the figure handle in the xolotl object
    xHandles.individual = individual;
else
    % Extract the individual structure from the xolotl object
    individual = xHandles.individual;

    % Retrieve from the handles structure
    %   Note: tVec & dataToCompare must be retrieved in full 
    %           (not just what's plotted)
    tVec = xHandles.tVec;
    dataToCompare = xHandles.dataToCompare;
    compToPatch = xHandles.compToPatch;
    idxCompToPatch = xHandles.idxCompToPatch;
    timeToStabilize = xHandles.timeToStabilize;
    holdingPotential = xHandles.holdingPotential;
    externalCurrentOrig = xHandles.externalCurrentOrig;
end

% Extract information from the individual structure
yLimits = individual.yLimits;
baseWindow = individual.baseWindow;
fitWindow = individual.fitWindow;
endPointsToPlot = individual.endPointsToPlot;
subPlots = individual.subPlots;

%% Set up simulation
% Reset to the original external current
xolotlObject.I_ext = externalCurrentOrig;

% Find the holding current (nA) necessary to match the holding potential
holdingCurrent = ...
    xolotl_estimate_holding_current(xolotlObject, holdingPotential, ...
                                    'CompToPatch', compToPatch, ...
                                    'TimeToStabilize', timeToStabilize);

% Add the holding current (nA)
xolotl_add_holding_current(xolotlObject, 'Compartment', compToPatch, ...
                            'Amplitude', holdingCurrent);

% Retrieve the new external current
externalCurrent = xolotlObject.I_ext;

%% Simulate
% Get voltage traces for all compartments
vVecs = xolotlObject.integrate;

%% Computed updated data for plots
% Reorganize voltage traces to match plots
vVecPlots = vVecs(:, idxInXolotl);

% Put together into a data array
dataPlots = [vVecPlots, externalCurrent(:, idxCompToPatch)];

% Re-compute baseline errors if not provided
if isempty(baseErrors)
    if isempty(dataToCompare)
        % Just use the baseline noise
        baseErrors = baseNoise;
    else
        % Compute sweep errors over the baseline window
        errorStructTemp = compute_sweep_errors(dataPlots, dataToCompare, ...
                                'TimeVecs', tVec, 'FitWindow', baseWindow, ...
                                'NormalizeError', false);

        % Extract baseline errors for each trace
        baseErrors = errorStructTemp.swpErrors;
    end
end

% Compute sweep errors if not provided
if isempty(sweepErrors)
    if isempty(dataToCompare)
        % Compute the baseline noise over the fitting window
        sweepErrors = compute_baseline_noise(dataPlots, tVec, fitWindow);
    else
        % Compute sweep errors over the fitting window
        errorStructTemp = compute_sweep_errors(dataPlots, dataToCompare, ...
                                'TimeVecs', tVec, 'FitWindow', fitWindow, ...
                                'NormalizeError', false);

        % Extract sweep errors for each trace
        sweepErrors = errorStructTemp.swpErrors;
    end
end

% Re-generate the error strings
if ~isempty(baseErrors)
    errorStrings = ...
        arrayfun(@(x) ['Noise = ', num2str(baseErrors(x), nSigFig), '; ', ...
                        'RMSE = ', num2str(sweepErrors(x), nSigFig)], ...
                    1:nTraces, 'UniformOutput', false);
else
    errorStrings = ...
        arrayfun(@(x) ['Noise = NaN; ', ...
                        'RMSE = ', num2str(sweepErrors(x), nSigFig)], ...
                    1:nTraces, 'UniformOutput', false);
end

%% Update plots
% Restrict to same end points to be consistent with the default plot
%   Note: this should only be done after the errors are computed
dataPlots = extract_subvectors(dataPlots, 'EndPoints', endPointsToPlot);

% Compute new y limits
%   Note: this should only be done after dataPlots is restricted
newYLimits = zeros(nTraces, 2);
for iTrace = 1:nTraces
    % Put the new data and old y limits together
    newDataForLimits = {dataPlots(:, iTrace); yLimits(iTrace, :)};

    % Update y limits
    newYLimits(iTrace, :) = compute_axis_limits(newDataForLimits, 'y');
end

% Update data in the chart line objects
for iTrace = 1:nTraces
    individual.plotsData(iTrace).YData = dataPlots(:, iTrace);
end

% Update data in the primitive line objects and set new y limits
for iTrace = 1:nTraces
    newYLimitsThis = newYLimits(iTrace, :);
    individual.boundaries(iTrace, 1).YData = newYLimitsThis;
    individual.boundaries(iTrace, 2).YData = newYLimitsThis;
    if ~any(isnan(newYLimitsThis))
        set(subPlots(iTrace), 'YLim', newYLimitsThis);
    end
end

% Update subplot titles
for iTrace = 1:nTraces
    individual.subTitles(iTrace).String = errorStrings{iTrace};
end

%% Save handles in xolotl object if under x.manipulate
if isstruct(xolotlObject.handles)
    xolotlObject.handles = xHandles;
end

%% Reset to the original external current
% Note: this is important for x.manipulate to work!
xolotlObject.I_ext = externalCurrentOrig;

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
vVecSoma = vVecs(:, idxCompToPatch);

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

% Extract dataToCompare
dataToCompare = zeros(nSamples, nCompartments);
for iTrace = 1:nCompartments   
    dataToCompare(:, iTrace) = individual.plotsDataToCompare(iTrace).YData;
end

% Restrict to x limits for faster processing
if ~isempty(xLimits) && isnumeric(xLimits)
    % Find the end points
    endPointsToPlot = find_window_endpoints(xLimits, tVec);

    % Restrict to these end points
    [tVec, iStim, dataToCompare] = ...
        argfun(@(x) extract_subvectors(x, 'EndPoints', endPointsToPlot), ...
                tVindividual.plotsData(iTrace).YData = dataPlots{iTrace};
ec, iStim, dataToCompare);
else
    % Use the first and last indices
    endPointsToPlot = find_window_endpoints([], tVec);
end

sweepWeights = individual.sweepWeights;

% holdingCurrent = 0.0743;
% holdingCurrent = 0.1125;
% holdingCurrent = 0;
% holdingCurrent = 0.095;
% holdingCurrent = 0.8;

%                   - 'CpDelay': current pulse delay (ms)
%                   must be a numeric scalar
%                   default == 
%                   - 'CpDuration': current pulse duration (ms)
%                   must be a numeric scalar
%                   default == 
%                   - 'CpAmplitude': current pulse amplitude (nA)
%                   must be a numeric scalar
%                   default == 
%                   - 'TimeEndCpr': time end of current pulse response (ms)
%                   must be a numeric scalar
%                   default == 
cpDelayDefault = [];            % set later
cpDurationDefault = [];         % set later
cpAmplitudeDefault = [];        % set later
timeEndCprDefault = [];         % set later
addParameter(iP, 'CpDelay', cpDelayDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CpDuration', cpDurationDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CpAmplitude', cpAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'TimeEndCpr', timeEndCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
cpDelay = iP.Results.CpDelay;
cpDuration = iP.Results.CpDuration;
cpAmplitude = iP.Results.CpAmplitude;
timeEndCpr = iP.Results.TimeEndCpr;

% Recompute the number of samples
nSamples = count_samples(tVec);

% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract the compartments
compartments = parsedParams.compartments;
nCompartments = parsedParams.nCompartments;

% Use just the first nCompartments vectors
vVecToCompare = dataToCompare(:, 1:nCompartments);

% Use just the first nCompartments windows or end points
[baseWindow, fitWindow, endPointsToPlot] = ...
    argfun(@(x) x(1:nCompartments), baseWindow, fitWindow, endPointsToPlot);

% Extract tVec from the first plot
tVec = transpose(individual.plotsData(1).XData);

individual.plotsData(iTrace).YData = dataPlots{iTrace};

if ~isempty(holdingPotential)
else
    holdingCurrent = 0;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
