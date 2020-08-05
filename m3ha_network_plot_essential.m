function handles = m3ha_network_plot_essential (varargin)
%% Plots essential traces from special cells in the network
% Usage: handles = m3ha_network_plot_essential (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - handles to plotted objects
%                   specified as a scalar structure
%
% Arguments:
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'essential'     - essential plots
%                       'm2h'           - m2h plot
%                   default == 'essential'
%                   - 'InFolder': directory containing the .singsp files
%                   must be a string scalar or a character vector
%                   default == pwds
%                   - 'AmpScaleFactor': amplitude scaling factor
%                   must be a numeric scalar
%                   default == 200%
%                   - 'PharmCondition': pharmacological condition
%                   must be a numeric scalar
%                   default == 1
%                   - 'XLimits': x value limits
%                   must be empty or a numeric vector of 2 elements
%                   default == []
%                   - 'OutFolder': output folder
%                   must be a string scalar or a character vector
%                   default == inFolder
%                   - 'FigTitle': figure title
%                   must be a string scalar or a character vector
%                   default == ['Essential traces for ', commonSuffix]
%                   - 'FigName': figure path for saving
%                   must be a string scalar or a character vector
%                   default == [commonSuffix, '_essential.png']
%                   - 'SaveNewFlag': whether to create and save new figure
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/compute_total_current.m
%       cd/convert_units.m
%       cd/extract_columns.m
%       cd/extract_fileparts.m
%       cd/read_neuron_outputs.m
%       cd/plot_traces.m
%       cd/plot_window_boundaries.m
%       cd/set_default_flag.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/m3ha_plot_figure08.m

% File History:
% 2020-01-30 Modified from m3ha_network_plot_gabab.m
% 2020-02-06 Added 'XLimits' as an optional argument
% 2020-02-06 Now downsamples vectors
% 2020-02-06 Added 'PlotType' as an optional argument

%% Hard-coded parameters
validPlotTypes = {'essential', 'm2h'};
specialExtension = 'singsp';
cellIdRT = 0;
cellIdTC = 0;

%   Note: Must be consistent with m3ha_network_launch.m
timeToStabilize = 2000;         % time for everything to stabilize

% Column numbers for simulated data
%   Note: Must be consistent with m3ha_net.hoc
RT_TIME = 1;
RT_VOLT = 2;
RT_INA = 3;
RT_IK = 4;
RT_ICA = 5;
RT_IAMPA = 6;
RT_IGABAA = 7;
RT_CAI = 8;
RT_CLI = 9;

TC_TIME = 1;
TC_VOLT = 2;
TC_IN = 3;
TC_IK = 4;
TC_ICA_SOMA = 5;
TC_IGABAA = 6;
TC_IGABAB = 7;
TC_CAI = 8;
TC_GGABAB = 9;
TC_IT_M_DEND2 = 10;
TC_IT_MINF_DEND2 = 11;
TC_IT_H_DEND2 = 12;
TC_IT_HINF_DEND2 = 13;
TC_ICA_DEND1 = 14;
TC_ICA_DEND2 = 15;
TC_GGABAA = 16;

% Plot parameters
xLabel = 'Time (ms)';
pharmLabels = {'{\it s}-Control', '{\it s}-GAT1 Block', ...
                '{\it s}-GAT3 Block', '{\it s}-Dual Block'};

tcParamsPrefix = 'TCparams';
simParamsPrefix = 'sim_params';
itm2hDiffLowerLimit = 1e-9;
itm2hDiffThreshold = 1e-2;

% TODO: Make optional arguments
figTypes = 'png';
tcParamsTable = table.empty;
simParamsTable = table.empty;

%% Default values for optional arguments
plotTypeDefault = 'essential';
inFolderDefault = pwd;      % use current directory by default
ampScaleFactorDefault = []; % set later
pharmConditionDefault = []; % set later
xLimitsDefault = timeToStabilize + [0, 8000];
outFolderDefault = '';      % set later
figTitleDefault = '';           % set later
figNameDefault = '';        % no figure name by default
saveNewFlagDefault = true;  % create and save new figure by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'AmpScaleFactor', ampScaleFactorDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) && isscalar(x), ...
                ['AmpScaleFactor must be either empty ', ...
                    'or a numeric scalar!']));
addParameter(iP, 'PharmCondition', pharmConditionDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) && isscalar(x), ...
                ['PharmCondition must be either empty ', ...
                    'or a numeric scalar!']));
addParameter(iP, 'XLimits', xLimitsDefault);
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigTitle', figTitleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveNewFlag', saveNewFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
inFolder = iP.Results.InFolder;
ampScaleFactor = iP.Results.AmpScaleFactor;
pharmCondition = iP.Results.PharmCondition;
xLimits = iP.Results.XLimits;
outFolder = iP.Results.OutFolder;
figTitle = iP.Results.FigTitle;
figName = iP.Results.FigName;
saveNewFlag = iP.Results.SaveNewFlag;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% Set default parameters
if isempty(ampScaleFactor)
    ampScaleFactor = 200;
end
if isempty(pharmCondition)
    pharmCondition = 1;
end

% Decide whether RT data is needed
loadRT = set_default_flag([], strcmp(plotType, 'essential'));

% Set default output folder
if isempty(outFolder)
    outFolder = inFolder;
end

% Set up the special file prefixes
spPrefixRT = sprintf('RE[%d]', cellIdRT);
spPrefixTC = sprintf('TC[%d]', cellIdTC);

% Find the appropriate file keyword
ampScaleFactorNetwork = ampScaleFactor / 12;
spKeyword = ['pCond_', num2str(pharmCondition), '_', ...
            'gIncr_', num2str(ampScaleFactorNetwork)];

% Locate the RT neuron data
if loadRT
    [~, dataPathRT] = all_files('Directory', inFolder, ...
                                'Prefix', spPrefixRT, 'Keyword', spKeyword, ...
                                'Extension', specialExtension, 'MaxNum', 1);
end

% Locate the TC neuron data
[~, dataPathTC] = all_files('Directory', inFolder, ...
                            'Prefix', spPrefixTC, 'Keyword', spKeyword, ...
                            'Extension', specialExtension, 'MaxNum', 1);

% Decide on the TC parameters table
if isempty(tcParamsTable)
    % Find the corresponding parameters file
    [~, tcParamsPath] = ...
        all_files('Directory', inFolder, 'Keyword', spKeyword, ...
                    'Prefix', tcParamsPrefix, 'Extension', 'csv', ...
                    'MaxNum', 1);

    % Load the simulation parameters table
    tcParamsTable = readtable(tcParamsPath);
end

% Decide on the simulation parameters table
if isempty(simParamsTable)
    % Find the corresponding parameters file
    [~, simParamsPath] = ...
        all_files('Directory', inFolder, 'Keyword', spKeyword, ...
                    'Prefix', simParamsPrefix, 'Extension', 'csv', ...
                    'MaxNum', 1);

    % Load the simulation parameters table
    simParamsTable = readtable(simParamsPath, 'ReadRowNames', true);
end

% Decide on figure name
if isempty(figName) && saveNewFlag
    commonSuffix = extract_fileparts({dataPathTC, dataPathRT}, 'commonsuffix');
    figName = fullfile(outFolder, [commonSuffix, '_', plotType, '.png']);
end

% Decide on figure title
if isempty(figTitle)
    switch plotType
    case 'essential'
        figTitle = ['Essential traces for ', commonSuffix];
    case 'm2h'
        figTitle = ['m2h for ', commonSuffix];            
    end
    figTitle = replace(figTitle, '_', '\_');
end

%% Do the job
% Extract stimulation start and duration in ms
stimStartMs = simParamsTable{'stimStart', 'Value'};
stimDurMs = simParamsTable{'stimDur', 'Value'};

% Construct stimulation window
stimWindow = [stimStartMs, stimStartMs + stimDurMs];

% Load simulated data
if loadRT
    [simDataRT, ~, nColumnsRT] = read_neuron_outputs('FileNames', dataPathRT);
end
[simDataTC, ~, nColumnsTC] = read_neuron_outputs('FileNames', dataPathTC);

% Decide whether to plot gGABAA
plotGGABAA = set_default_flag([], nColumnsTC >= TC_GGABAA);

% Convert the table to a structure array
tcParamsStructArray = table2struct(tcParamsTable);

% Extract vectors from simulated data
if loadRT
    vVecRT = extract_columns(simDataRT, RT_VOLT);
end
switch plotType
case 'essential'
    [tVecsMs, vVecTC, gGababTCUs, itSomaTC, itDend1TC, itDend2TC, ...
            itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
        extract_columns(simDataTC, [TC_TIME, TC_VOLT, TC_GGABAB, ...
                                    TC_ICA_SOMA, TC_ICA_DEND1, TC_ICA_DEND2, ...
                                    TC_IT_M_DEND2, TC_IT_MINF_DEND2, ...
                                    TC_IT_H_DEND2, TC_IT_HINF_DEND2]);
    if plotGGABAA
        gGabaaTCUs = extract_columns(simDataTC, TC_GGABAA);
    end
case 'm2h'
    [tVecsMs, itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
        extract_columns(simDataTC, [TC_TIME, TC_IT_M_DEND2, TC_IT_MINF_DEND2, ...
                                    TC_IT_H_DEND2, TC_IT_HINF_DEND2]);
end

% Clear simData to release memory
clear simDataRT simDataTC

% Downsample by 10;
switch plotType
case 'essential'
    [tVecsMs, vVecRT, vVecTC, gGababTCUs, ...
            itSomaTC, itDend1TC, itDend2TC, ...
            itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
        argfun(@(x) downsample(x, 10), ...
                tVecsMs, vVecRT, vVecTC, gGababTCUs, ...
                itSomaTC, itDend1TC, itDend2TC, ...
                itmDend2, itminfDend2, ithDend2, ithinfDend2);
    if plotGGABAA
        gGabaaTCUs = downsample(gGabaaTCUs, 10);
    end
case 'm2h'
    [tVecsMs, itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
        argfun(@(x) downsample(x, 10), ...
                tVecsMs, itmDend2, itminfDend2, ithDend2, ithinfDend2);
end

% Compute m2h and its steady state
itm2h = (itmDend2 .^ 2) .* ithDend2;
itminf2hinf = (itminfDend2 .^ 2) .* ithinfDend2;
itm2hDiff = itm2h - itminf2hinf;
itm2hDiff(itm2hDiff < itm2hDiffLowerLimit) = itm2hDiffLowerLimit;

switch plotType
case 'essential'
    % Convert conductance from uS to nS
    gGababTCNs = convert_units(gGababTCUs, 'uS', 'nS');
    if plotGGABAA
        gGabaaTCNs = convert_units(gGabaaTCUs, 'uS', 'nS');
    end

    % Compute total T current
    itTotalTC = compute_total_current([itSomaTC, itDend1TC, itDend2TC], ...
                                        'GeomParams', tcParamsStructArray);

    % List all possible items to plot
    vecsAll = {vVecRT; vVecTC; gGababTCNs; itTotalTC; ...
                itm2h; itminf2hinf; itm2hDiff};

    % List corresponding labels
    labelsAll = {'V_{RT} (mV)'; 'V_{TC,soma} (mV)'; ...
                'g_{GABA_B} (nS)'; ...
                'I_{T} (nA)'; 'm_{T,dend2}^2h_{T,dend2}'; ...
                'm_{\infty,T,dend2}^2h_{\infty,T,dend2}'; ...
                'm2hDiff'};

    if plotGGABAA
        % List all possible items to plot
        vecsAll = vertcat(vecsAll, {gGabaaTCNs});

        % List corresponding labels
        labelsAll = vertcat(labelsAll, {'g_{GABA_A} (nS)'});
    end
case 'm2h'
    % List all possible items to plot
    vecsAll = {itm2h; itminf2hinf};

    % List corresponding labels
    labelsAll = {'m_{T,dend2}^2h_{T,dend2}'; ...
                'm_{\infty,T,dend2}^2h_{\infty,T,dend2}'};    
end

% Create a figure
if saveNewFlag
    fig = set_figure_properties('AlwaysNew', true);
end

switch plotType
case 'essential'
    % Plot traces
    handles = plot_traces(tVecsMs, vecsAll, ...
                            'PlotMode', 'parallel', 'XLimits', xLimits, ...
                            'XLabel', 'suppress', 'YLabel', labelsAll, ...
                            'FigTitle', figTitle, 'LegendLocation', 'suppress', ...
                            'FigName', figName, 'FigTypes', figTypes, ...
                            otherArguments);

    % Extract the axes handles for the subplots
    subPlots = handles.subPlots;

    % Make the 5-7th subplots log-scaled
    arrayfun(@(x) set(subPlots(x), 'YScale', 'log'), 5:7);

    % Add a threshold line
    subplot(subPlots(7));
    plot_horizontal_line(itm2hDiffThreshold, 'ColorMap', 'DarkGreen', ...
                            'LineStyle', ':', 'LineWidth', 1);

    % Plot stimulation boundaries
    for i = 1:numel(subPlots)
        subplot(subPlots(i));
        plot_window_boundaries(stimWindow, 'BoundaryType', 'verticalShades', ...
                                'Color', 'PaleGreen');
    end
case 'm2h'
    handlesInstantaneous = ...
        plot_traces(tVecsMs, itm2h, ...
                    'LineStyle', '-', 'PlotMode', 'overlapped', ...
                    'LegendLocation', 'suppress', ...
                    'XLimits', xLimits, ...
                    'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                    'YLabel', 'm^2h', 'FigTitle', figTitle, otherArguments);
    hold on;
    handlesSteadyState = ...
        plot_traces(tVecsMs, itminf2hinf, 'PlotOnly', true, ...
                    'LineStyle', ':', 'PlotMode', 'overlapped', ...
                    'XLimits', xLimits, otherArguments);

    % set(gca, 'YLim', [1e-6, 1]);

    % Set the y axis to be log-scaled
    set(gca, 'YScale', 'log');

    % Plot stimulation boundaries
    plot_window_boundaries(stimWindow, 'BoundaryType', 'verticalShades', ...
                            'Color', 'PaleGreen');

    handles.handlesInstantaneous = handlesInstantaneous;
    handles.handlesSteadyState = handlesSteadyState;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
