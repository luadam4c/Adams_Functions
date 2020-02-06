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
%       varargin    - 'InFolder': directory containing the .singsp files
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
%       cd/load_neuron_outputs.m
%       cd/plot_traces.m
%       cd/plot_window_boundaries.m
%       cd/set_figure_properties.m
%
% Used by:
%       cd/m3ha_plot_figure07.m

% File History:
% 2020-01-30 Modified from m3ha_network_plot_gabab.m
% 2020-02-06 Added 'XLimits' as an optional argument

%% Hard-coded parameters
spExtension = 'singsp';
cellIdRT = 0;
cellIdTC = 0;
ipscStartMs = 3000;

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

% Plot parameters
xLabel = 'Time (ms)';
pharmLabels = {'{\it s}-Control', '{\it s}-GAT1 Block', ...
                '{\it s}-GAT3 Block', '{\it s}-Dual Block'};

tcParamsPrefix = 'TCparams';
simParamsPrefix = 'sim_params';

% TODO: Make optional arguments
figTypes = 'png';
tcParamsTable = table.empty;
simParamsTable = table.empty;

%% Default values for optional arguments
inFolderDefault = pwd;      % use current directory by default
ampScaleFactorDefault = []; % set later
pharmConditionDefault = []; % set later
xLimitsDefault = [2000, 10000];
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
[~, dataPathRT] = all_files('Directory', inFolder, ...
                            'Prefix', spPrefixRT, 'Keyword', spKeyword, ...
                            'Extension', spExtension, 'MaxNum', 1);

% Locate the TC neuron data
[~, dataPathTC] = all_files('Directory', inFolder, ...
                            'Prefix', spPrefixTC, 'Keyword', spKeyword, ...
                            'Extension', spExtension, 'MaxNum', 1);

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
    figName = fullfile(outFolder, [commonSuffix, '_essential.png']);
end

% Decide on figure title
if isempty(figTitle)
    figTitle = ['Essential traces for ', commonSuffix];
    figTitle = replace(figTitle, '_', '\_');
end

%% Do the job
% Extract stimulation start and duration in ms
stimStartMs = simParamsTable{'stimStart', 'Value'};
stimDurMs = simParamsTable{'stimDur', 'Value'};

% Construct stimulation window
stimWindow = [stimStartMs, stimStartMs + stimDurMs];

% Load simulated data
[simDataRT, simDataTC] = ...
    argfun(@(x) load_neuron_outputs('FileNames', x), dataPathRT, dataPathTC);

% Convert the table to a structure array
tcParamsStructArray = table2struct(tcParamsTable);

% Extract vectors from simulated data
[tVecsMs, vVecRT] = ...
    extract_columns(simDataRT, [RT_TIME, RT_VOLT]);
[vVecTC, gCmdTCUs, itSomaTC, itDend1TC, itDend2TC, ...
        itmDend2, itminfDend2, ithDend2, ithinfDend2] = ...
    extract_columns(simDataTC, [TC_VOLT, TC_GGABAB, ...
                                TC_ICA_SOMA, TC_ICA_DEND1, TC_ICA_DEND2, ...
                                TC_IT_M_DEND2, TC_IT_MINF_DEND2, ...
                                TC_IT_H_DEND2, TC_IT_HINF_DEND2]);

% Convert conductance from uS to nS
gCmdTCNs = convert_units(gCmdTCUs, 'uS', 'nS');

% Compute total T current
itTotalTC = ...
    compute_total_current([itSomaTC, itDend1TC, itDend2TC], ...
                            'GeomParams', tcParamsStructArray);

% Compute m2h
itm2hDend2 = (itmDend2 .^ 2) .* ithDend2;
itminf2hinfDend2 = (itminfDend2 .^ 2) .* ithinfDend2;

% Clear simData to release memory
clear simDataRT simDataTC

% List all possible items to plot
vecsAll = {vVecRT; vVecTC; gCmdTCNs; itTotalTC; itm2hDend2; itminf2hinfDend2};

% List corresponding labels
labelsAll = {'V_{RT} (mV)'; 'V_{TC,soma} (mV)'; 'g_{GABA_B} (nS)'; ...
            'I_{T} (nA)'; 'm_{T,dend2}^2h_{T,dend2}'; ...
            'm_{\infty,T,dend2}^2h_{\infty,T,dend2}'};

% Create a figure
if saveNewFlag
    fig = set_figure_properties('AlwaysNew', true);
end

% Plot traces
handles = plot_traces(tVecsMs, vecsAll, ...
                        'PlotMode', 'parallel', 'XLimits', xLimits, ...
                        'XLabel', 'suppress', 'YLabel', labelsAll, ...
                        'FigTitle', figTitle, 'LegendLocation', 'suppress', ...
                        'FigName', figName, 'FigTypes', figTypes, ...
                        otherArguments);

% Extract the axes handles for the subplots
subPlots = handles.subPlots;

% Make the 5th and 6th subplot log-scaled
arrayfun(@(x) set(subPlots(x), 'YScale', 'log'), 5:6);

% Plot stimulation boundaries
for i = 1:numel(subPlots)
    subplot(subPlots(i));
    plot_window_boundaries(stimWindow, 'BoundaryType', 'verticalShades', ...
                            'Color', 'PaleGreen');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
