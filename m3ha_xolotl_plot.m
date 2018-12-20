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
% 

%% Hard-coded parameters

%% Default values for optional arguments
dataToCompareDefault = [];      % no data to compare against by default
xLimitsDefault = [];            % for current pulse
colorMapDefault = [];           % set in m3ha_plot_individual_traces.m
figTitleDefault = 'Simulation by xolotl';
figNumberDefault = 104;         % set in m3ha_plot_individual_traces.m
figNameDefault = '';            % don't save figure by default
baseWindowDefault = [];         % set in m3ha_plot_individual_traces.m
fitWindowDefault = [];          % set in m3ha_plot_individual_traces.m
baseNoiseDefault = [];          % set in m3ha_plot_individual_traces.m
sweepWeightsDefault = [];       % set in m3ha_plot_individual_traces.m
sweepErrorsDefault = [];        % set in m3ha_plot_individual_traces.m
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
sweepErrors = iP.Results.SweepErrors;
plotSwpWeightsFlag = iP.Results.PlotSwpWeightsFlag;

%% Preparation
% Extract the compartment names
compNames = xolotlObject.Children;

% Find the idx for soma
[idxSoma, idxDend1, idxDend2] = ...
    argfun(@(x) find_ind_str_in_cell(x, compNames, 'IgnoreCase', true, ...
                                'SearchMode', 'substrings', 'MaxNum', 1), ...
            'soma', 'dend1', 'dend2');

%% Retrieve handles from xolotl object
xHandles = xolotlObject.handles;

%% Initialize a plot if not already done
% TODO: Make this a function?
if isempty(xHandles) || ~isfield(xHandles, 'individual')
    % Extract the time step in ms
    timeStep = xolotlObject.dt;

    % Extract the external current injection protocol
    currentProtocol = xolotlObject.I_ext;

    % Save the stimulation pulse
    iStim = currentProtocol(:, idxSoma);

    % Get the number of samples for each trace
    nSamples = count_samples(iStim);

    % Create a time vector in milliseconds
    tVecs = create_time_vectors(nSamples, 'SamplingIntervalMs', timeStep, ...
                                'TimeUnits', 'ms');

    % Restrict to x limits for faster processing
    if ~isempty(xLimits) && isnumeric(xLimits)
        % Find the end points
        endPoints = find_window_endpoints(xLimits, tVecs);

        % Restrict to these end points
        [tVecs, iStim, dataToCompare] = ...
            argfun(@(x) extract_subvectors(x, 'EndPoints', endPoints), ...
                    tVecs, iStim, dataToCompare);
    else
        % Use the first and last end point
        endPoints = find_window_endpoints([], tVecs);
    end

    % Create NaN data for the initial plot
    nanVoltageData = NaN * tVecs;
    nanData = [nanVoltageData, nanVoltageData, nanVoltageData, iStim];

    % Make a figure
    figHandle = ...
        figure('Outerposition', [100, 100, 600, 600], ...
            'PaperUnits', 'points', 'PaperSize', [1200, 600]); hold on

    % Initialize the plot by using NaNs
    individual = m3ha_plot_individual_traces(tVecs, nanData, ...
                                    'DataToCompare', dataToCompare, ...
                                    'XLimits', xLimits, ...
                                    'ColorMap', colorMap, ...
                                    'FigTitle', figTitle, ...
                                    'FigHandle', figHandle, ...
                                    'BaseWindow', baseWindow, ...
                                    'FitWindow', fitWindow, ...
                                    'BaseNoise', baseNoise, ...
                                    'SweepWeights', sweepWeights, ...
                                    'SweepErrors', sweepErrors, ...
                                    'PlotSwpWeightsFlag', plotSwpWeightsFlag);

    % If using x.manipulate, do the following
    if isfield(xHandles, 'puppeteer_object')
        % Attach figure to puppeteer so that 
        %   "the figure can be closed automatically"
        xHandles.puppeteer_object.attachFigure(individual.fig);

        % Prevent the time axis limits from automatically updating

    end

    % Store endpoints
    individual.endPoints = endPoints;

    % Store the figure handle in the xolotl object
    xHandles.individual = individual;
end

%% Simulate
% Get voltage traces for all compartments
vVecs = xolotlObject.integrate;

% Restrict to same end points to be consistent with the default plot
vVecs = extract_subvectors(vVecs, 'EndPoints', xHandles.individual.endPoints);

%% Computed updated data for plots
% Extract the voltage traces for each compartment
vVecSoma = vVecs(:, idxSoma);
vVecDendrite1 = vVecs(:, idxDend1);
vVecDendrite2 = vVecs(:, idxDend2);

%% Update plots
% Update data in the corresponding line object
xHandles.individual.plotsData(1).YData = vVecDendrite2;
xHandles.individual.plotsData(2).YData = vVecDendrite1;
xHandles.individual.plotsData(3).YData = vVecSoma;

% Update the 
xHandles.individual.boundaries(1, :).YData = 
xHandles.individual.boundaries(2, :).YData = 
xHandles.individual.boundaries(3, :).YData = 

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
    m3ha_plot_individual_traces(tVecs, data, ...
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%