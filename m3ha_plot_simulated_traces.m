function handles = m3ha_plot_simulated_traces (varargin)
%% Plots simulated traces from NEURON output files
% Usage: handles = m3ha_plot_simulated_traces (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'individual'    - voltage traces on separate subplots
%                       'residual'      - residual traces 
%                                               between simulated and recorded
%                       'overlapped'    - all traces of interest overlapped
%                       'm2h'           - m2h plot
%                   default == 'individual'
%                   - 'BuildMode': TC neuron build mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - inserted leak channels only
%                       'active'  - inserted both passive and active channels
%                   default == detected
%                   - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulated a current pulse response
%                       'active'  - simulated an IPSC response
%                   default == detected
%                   - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == set in all_files.m
%                   - 'FileNames': paths to simulated data
%                   must be a string array or a cell array of character vectors
%                   default == detected from Directory
%                   - 'Extension': data file extension
%                   must be a string scalar or a character vector
%                   default == 'out'
%                   - 'ColorMap': color map
%                   must be TODO
%                   default == TODO
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == no restrictions
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == same as Directory
%                   - 'ExpStr': experiment string file names
%                   must be a character array
%                   default == extract_common_prefix(fileNames)
%                   - 'tVecs': time vectors to match
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'vVecsRec': recorded voltage vectors
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'iVecsRec': recorded current vectors
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'gVecsRec': recorded conductance vectors
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'Residuals': voltage residuals
%                   must be a numeric array or a cell array of numeric arrays
%                   default == [] (none provided)
%                   - 'LineWidth': line width of plots
%                   must be empty or a positive scalar
%                   default == TODO
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/construct_fullpath.m
%       cd/decide_on_colormap.m
%       cd/extract_columns.m
%       cd/extract_common_prefix.m
%       cd/isemptycell.m
%       cd/load_neuron_outputs.m
%       cd/m3ha_extract_sweep_names.m
%       cd/m3ha_import_raw_traces.m
%       cd/plot_traces.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-10-14 Created by Adam Lu
% 2019-12-22 Added 'PlotType' as an optional argument

%% Hard-coded parameters
validPlotTypes = {'individual', 'residual', 'overlapped', 'm2h'};
validBuildModes = {'', 'active', 'passive'};
validSimModes = {'', 'active', 'passive'};
maxRowsWithOneOnly = 8;
lineWidthParallel = 1;
lineWidthIndividual = 0.5;

% Note: Must be consistent with m3ha_neuron_run_and_analyze.m
importedSuffix = 'imported_files';

%% Column numbers for recorded data
%   Note: Must be consistent with m3ha_resave_sweeps.m
TIME_COL_REC = 1;
VOLT_COL_REC = 2;
CURR_COL_REC = 3;
COND_COL_REC = 4;

% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
DEND1_COL_SIM = 3;
DEND2_COL_SIM = 4;
IDCLAMP_COL_SIM = 5;
GGABAB_COL_SIM = 6;
iCP_COL_SIM = 7;
IEXT_COL_SIM = 8;

%% Default values for optional arguments
plotTypeDefault = 'individual';
buildModeDefault = '';          % set later
simModeDefault = '';            % set later
directoryDefault = '';          % set in all_files.m
fileNamesDefault = {};
extensionDefault = 'out';       % 
colorMapDefault = [];
% xLimitsDefault = [];          % set later
xLimitsDefault = [2800, 4500];
outFolderDefault = '';          % set later
expStrDefault = '';             % set later
tVecsDefault = [];
vVecsRecDefault = [];
iVecsRecDefault = [];
gVecsRecDefault = [];
residualsDefault = [];
lineWidthDefault = [];      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'BuildMode', buildModeDefault, ...
    @(x) any(validatestring(x, validBuildModes)));
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fileNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ExpStr', expStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'tVecs', tVecsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['tVecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'vVecsRec', vVecsRecDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['vVecsRec must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'iVecsRec', iVecsRecDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['iVecsRec must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'gVecsRec', gVecsRecDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['gVecsRec must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'Residuals', residualsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['Residuals must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'LineWidth', lineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
buildMode = validatestring(iP.Results.BuildMode, validBuildModes);
simMode = validatestring(iP.Results.SimMode, validSimModes);
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
extension = iP.Results.Extension;
colorMap = iP.Results.ColorMap;
xLimits = iP.Results.XLimits;
outFolder = iP.Results.OutFolder;
expStr = iP.Results.ExpStr;
tVecs = iP.Results.tVecs;
vVecsRec = iP.Results.vVecsRec;
iVecsRec = iP.Results.iVecsRec;
gVecsRec = iP.Results.gVecsRec;
residuals = iP.Results.Residuals;
lineWidth = iP.Results.LineWidth;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% Use the present working directory for both inputs and output by default
if isempty(directory)
    directory = pwd;
end

% Set default output directory
if isempty(outFolder)
    outFolder = directory;
end

% Decide on input paths
if isempty(fileNames)
    [~, fileNames] = ...
        all_files('Directory', directory, 'Extension', extension);
else
    % Make sure they are full paths
    fileNames = construct_fullpath(fileNames, 'Directory', directory);
end

% Use the common expStr as the experiment string
if isempty(expStr)
    expStr = extract_common_prefix(fileNames);
end

% Decide on the build mode
% TODO: Detect from output structure instead
if isempty(buildMode)
    if all(contains(fileNames, 'cpr'))
        buildMode = 'passive';
    else
        buildMode = 'active';
    end
end

% Decide on the simulation mode
% TODO: Detect from output structure instead
if isempty(simMode)
    if all(contains(fileNames, 'cpr'))
        simMode = 'passive';
    else
        simMode = 'active';
    end
end

% Create an experiment identifier for title
expStrForTitle = replace(expStr, '_', '\_');

% Count the number of files
nFiles = numel(fileNames);

% Decide on nRows
nRows = decide_on_nrows(nFiles, maxRowsWithOneOnly);

% Decide on the color map if not provided
if isempty(colorMap)
    switch plotType
        case {'individual', 'residual'}
            % Decide on the color map for individual and residual plots
            colorMap = decide_on_colormap('r', nRows);
        case {'overlapped', 'm2h'}
            % Decide on the colors for parallel plots
            colorMap = decide_on_colormap([], 4);
            if nSweeps > nRows
                nColumns = ceil(nSweeps / nRows);
                nSlots = nColumns * nRows;
                colorMap = reshape(repmat(reshape(colorMap, 1, []), ...
                                    nColumns, 1), nSlots, 3);
            end
        otherwise
            % Decide on the colors for each row in the plots
            colorMap = decide_on_colormap(colorMap, nFiles);
    end
end

% Decide on the plot line width
if isempty(lineWidth)
    switch plotType
        case {'individual', 'residual'}
            lineWidth = lineWidthIndividual;
        case {'overlapped', 'm2h'}
            lineWidth = lineWidthParallel;
        otherwise
            error('plotType unrecognized!');
    end
end

%% Data
% Look for matching recorded data
[~, importedPath] = all_files('Prefix', expStr, 'Suffix', importedSuffix, ...
                                'MaxNum', 1);

% Look for matching recorded sweep names
if ~isempty(importedPath)
    % TODO FOR SHINSHIN: read_log.m
    % Read the log file
    logTable = readtable(importedPath, 'ReadVariableNames', false, ...
                            'Delimiter', ' ');

    % Extract sweep names
    sweepNames = logTable.Var1;
else
    sweepNames = m3ha_extract_sweep_names(fileNames);
end

% Import and extract from recorded data
if ~all(isemptycell(sweepNames))
    % Import recorded traces
    realData = m3ha_import_raw_traces(sweepNames, 'ImportMode', simMode, ...
                                    'Verbose', true, 'OutFolder', outFolder);

    % Extract vectors from recorded data
    %   Note: these will be empty if realData not provided
    [tVecs, vVecsRec, iVecsRec, gVecsRec] = ...
        extract_columns(realData, [TIME_COL_REC, VOLT_COL_REC, ...
                                    CURR_COL_REC, COND_COL_REC]);
end

% Load simulated data
% If recorded data provided (tVecs not empty at this point),
%   interpolate simulated data to match the time points of recorded data
% Note: This is necessary because CVODE (variable time step method) 
%       is applied in NEURON
simData = load_neuron_outputs('FileNames', fileNames, 'tVecs', tVecs);

% Extract vectors from simulated data
[tVecs, vVecsSim, gVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
    extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                    GGABAB_COL_SIM, IEXT_COL_SIM, ...
                    DEND1_COL_SIM, DEND2_COL_SIM]);

%% Plots
% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
[tVecs, vVecsRec, iVecsRec, gVecsRec, ...
    vVecsSim, iVecsSim, gVecsSim, vVecsDend1, vVecsDend2] = ...
    argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
            tVecs, vVecsRec, iVecsRec, gVecsRec, ...
            vVecsSim, iVecsSim, gVecsSim, vVecsDend1, vVecsDend2);

% Plot according to plot type
switch plotType
    case 'individual'
        handles = m3ha_plot_individual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    simMode, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    case 'residual'
        handles = m3ha_plot_residual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    residuals, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    case 'overlapped'
        handles = m3ha_plot_overlapped_traces(simData, vVecsRec, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    case 'm2h'
        handles = m3ha_plot_m2h(simData, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    otherwise
        error('plotType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nRows = decide_on_nrows(nFiles, maxRowsWithOneOnly)
%% Decide on the number of rows

% Decide on the number of rows
if mod(nFiles, 4) == 0
    nRows = 4;
elseif nFiles <= maxRowsWithOneOnly
    nRows = nFiles;
else
    nRows = floor(sqrt(nSweeps));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_individual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    simMode, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

% TODO
plotSwpWeightsFlag = false;

%% Preparation
% Decide on figure title
figTitle = sprintf('All traces for Experiment %s', expStrForTitle);

% Decide on the axes to be linked
if strcmp(simMode, 'passive')
    linkAxesOption = 'x';
else
    linkAxesOption = 'xy';
end

%% Do the job
% Print to standard output
fprintf('Plotting figure of individual voltage traces for %s ...\n', expStr);

% Plot the individual traces
handles = plot_fitted_traces(tVecs, vVecsSim, ...
            'DataToCompare', vVecsRec, 'PlotMode', 'parallel', ...
            'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
            'ColorMap', colorMap, 'XLimits', xLimits, ...
            'LineWidth', lineWidth, 'LinkAxesOption', linkAxesOption, ...
            'FigTitle', figTitle, 'PlotSwpWeightsFlag', plotSwpWeightsFlag, ...
            otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_residual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    residuals, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

% TODO
plotSwpWeightsFlag = false;

% Calculate voltage residuals (simulated - recorded) if necessary
if isempty(residuals) && ~isempty(vVecsRec)
    residuals = compute_residuals(vVecsSim, vVecsRec);
end

% Decide on figure title
figTitle = sprintf('Residuals for Experiment %s', expStrForTitle);

%% Do the job
% Print to standard output
fprintf('Plotting figure of residual traces for %s ...\n', expStr);

% Plot the individual traces
handles = plot_fitted_traces(tVecs, residuals, ...
            'PlotMode', 'residuals', ...
            'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
            'ColorMap', colorMap, 'XLimits', xLimits, ...
            'LineWidth', lineWidth, 'LinkAxesOption', 'xy', ...
            'FigTitle', figTitle, 'PlotSwpWeightsFlag', plotSwpWeightsFlag, ...
            otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_overlapped_traces (simData, vVecsRec, buildMode, ...
                                        xLimits, colorMap, lineWidth, ...
                                        expStr, expStrForTitle, otherArguments)

%% Hard-coded parameters
% TODO: Use this
labelsAll = {'Time (ms)'; 'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
        'V_{dend2} (mV)'; 'I_{GABA_B} (nA)'; 'g_{GABA_B} (uS)'; ...
        'I_{cp} (nA)'; 'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
        'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
        'I_{h} (mA/cm^2)'; 'm_{h}'; 'I_{A} (mA/cm^2)'; ...
        'm_{1,A}'; 'h_{1,A}'; 'm_{2,A}'; 'h_{2,A}'; ...
        'I_{Kir} (mA/cm^2)'; 'm_{\infty,Kir}'; ...
        'I_{NaP} (mA/cm^2)'; 'm_{\infty,NaP}'; 'h_{NaP}'};

% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
DEND1_COL_SIM = 3;
DEND2_COL_SIM = 4;
IDCLAMP_COL_SIM = 5;
GGABAB_COL_SIM = 6;
iCP_COL_SIM = 7;
IEXT_COL_SIM = 8;
ICA_COL_SIM = 9;
ITM_COL_SIM = 10;
ITMINF_COL_SIM = 11;
ITH_COL_SIM = 12;
ITHINF_COL_SIM = 13;
IH_COL_SIM = 14;
IHM_COL_SIM = 15;
IKA_COL_SIM = 16;
IAM1_COL_SIM = 17;
IAH1_COL_SIM = 18;
IAM2_COL_SIM = 19;
IAH2_COL_SIM = 20;
IKKIR_COL_SIM = 21;
IKIRM_COL_SIM = 22;
INAPNA_COL_SIM = 23;
INAPM_COL_SIM = 24;
INAPH_COL_SIM = 25;

%% Preparation
% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
if strcmpi(buildMode, 'passive')
    [tVecs, vVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        IEXT_COL_SIM, DEND1_COL_SIM, DEND2_COL_SIM]);
elseif strcmpi(buildMode, 'active')
    [tVecs, vVecsSim, gVecsSim, iVecsSim, icaVecsSim, ...
            itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
            ihVecsSim, ihmVecsSim, ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
            iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
            inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        GGABAB_COL_SIM, IEXT_COL_SIM, ...
                        ICA_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                        ITH_COL_SIM, ITHINF_COL_SIM, ...
                        IH_COL_SIM, IHM_COL_SIM, ...
                        IKA_COL_SIM, IAM1_COL_SIM, IAH1_COL_SIM, ...
                        IAM2_COL_SIM, IAH2_COL_SIM, ...
                        IKKIR_COL_SIM, IKIRM_COL_SIM, ...
                        INAPNA_COL_SIM, INAPM_COL_SIM, INAPH_COL_SIM]);
end

% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
if strcmpi(buildMode, 'passive')
    [tVecs, vVecsRec, vVecsSim, vVecsDend1, vVecsDend2, iVecsSim] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, vVecsRec, vVecsSim, vVecsDend1, vVecsDend2, iVecsSim);
elseif strcmpi(buildMode, 'active')
    [tVecs, vVecsRec, vVecsSim, gVecsSim, iVecsSim, ...
        icaVecsSim, itmVecsSim, itminfVecsSim, ...
        ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
        ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
        iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
        inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, vVecsRec, vVecsSim, gVecsSim, iVecsSim, ...
                icaVecsSim, itmVecsSim, itminfVecsSim, ...
                ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
                ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
                iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
                inapnaVecsSim, inapmVecsSim, inaphVecsSim);
end

% Select data to plot
if strcmpi(buildMode, 'passive')
    dataForOverlapped = {vVecsSim; vVecsDend1; vVecsDend2; iVecsSim};
elseif strcmpi(buildMode, 'active')
    dataForOverlapped = {vVecsSim; gVecsSim; iVecsSim; ...
            icaVecsSim; itm2hVecsSim; itminf2hinfVecsSim; ...
            itmVecsSim; itminfVecsSim; ithVecsSim; ithinfVecsSim; ...
            ihVecsSim; ikaVecsSim; ikkirVecsSim; inapnaVecsSim};
end

% Construct matching y labels
if strcmpi(buildMode, 'passive')
    yLabelsOverlapped = {'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                        'V_{dend2} (mV)'; 'I_{stim} (nA)'};
elseif strcmpi(buildMode, 'active')
    yLabelsOverlapped = {'V_{soma} (mV)'; 'g_{GABA_B} (uS)'; ...
            'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
            'm^2h_{T}'; 'm_{\infty}^2h_{\infty,T}'; ...
            'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
            'I_{h} (mA/cm^2)'; 'I_{A} (mA/cm^2)'; ...
            'I_{Kir} (mA/cm^2)'; 'I_{NaP} (mA/cm^2)'};
end

% Add recorded voltage on the top if exists
if ~isempty(vVecsRec)
    dataForOverlapped = [{vVecsRec}; dataForOverlapped];
    yLabelsOverlapped = [{'V_{rec} (mV)'}; yLabelsOverlapped];
end
        
% Construct matching time vectors
tVecsForOverlapped = repmat({tVecs}, size(dataForOverlapped));

% Decide on figure title and file name
figTitle = sprintf('Simulated traces for Experiment %s', expStrForTitle);

%% Plots
% Print to standard output
fprintf('Plotting figure of overlapped traces for %s ...\n', expStr);

% Plot overlapped traces
handles = ...
    plot_traces(tVecsForOverlapped, dataForOverlapped, ...
                'Verbose', false, 'PlotMode', 'parallel', ...
                'SubplotOrder', 'list', 'ColorMode', 'byTraceInPlot', ...
                'LegendLocation', 'suppress', ...
                'ColorMap', colorMap, 'XLimits', xLimits, ...
                'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                'YLabel', yLabelsOverlapped, ...
                'FigTitle', figTitle, 'LineWidth', lineWidth, ...
                otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_m2h (simData, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

%% Hard-coded parameters
% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
ITM_COL_SIM = 10;
ITMINF_COL_SIM = 11;
ITH_COL_SIM = 12;
ITHINF_COL_SIM = 13;

% Only do this for active mode
if strcmpi(buildMode, 'passive')
    handles = struct;
    return
end

%% Process data
% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
[tVecs, itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim] = ...
    extract_columns(simData, [TIME_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                    ITH_COL_SIM, ITHINF_COL_SIM]);

% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
[tVecs, itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim] = ...
    argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
            tVecs, itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim);

% Compute m2h
itm2hVecsSim = (itmVecsSim .^ 2) .* ithVecsSim;
itminf2hinfVecsSim = (itminfVecsSim .^ 2) .* ithinfVecsSim;

% Decide on figure title and file name
figTitle = sprintf('m2h for Experiment %s', expStrForTitle);

%% Plots
% Print to standard output
fprintf('Plotting figure of m2h for %s ...\n', expStr);

handlesInstantaneous = ...
    plot_traces(tVecs, itm2hVecsSim, ...
                'LineStyle', '-', 'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'LegendLocation', 'suppress', ...
                'ColorMap', colorMap, 'XLimits', xLimits, ...
                'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                'YLabel', 'm^2h', 'FigTitle', figTitle, otherArguments);
hold on;
handlesSteadyState = ...
    plot_traces(tVecs, itminf2hinfVecsSim, ...
                'PlotOnly', true, ...
                'LineStyle', '--', 'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'ColorMap', colorMap, 'XLimits', xLimits, otherArguments);

set(gca, 'YLim', [1e-6, 1]);
set(gca, 'YScale', 'log');

handles.handlesInstantaneous = handlesInstantaneous;
handles.handlesSteadyState = handlesSteadyState;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecs = prepare_for_plotting(vecs, endPointsForPlots)
%% Prepare vectors for plotting

% Restrict vectors to xLimits to save time on plotting
vecs = extract_subvectors(vecs, 'Endpoints', endPointsForPlots);

% Combine vectors into matrices
vecs = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

dataForOverlapped = {vVecsSim; gVecsSim; iVecsSim; ...
        icaVecsSim; itm2hVecsSim; itminf2hinfVecsSim; ...
        itmVecsSim; itminfVecsSim; ithVecsSim; ithinfVecsSim; ...
        ihVecsSim; ihmVecsSim; ...
        ikaVecsSim; iam1VecsSim; iah1VecsSim; ...
        iam2VecsSim; iah2VecsSim; ikkirVecsSim; ikirmVecsSim; ...
        inapnaVecsSim; inapmVecsSim; inaphVecsSim};
yLabelsOverlapped = {'V_{soma} (mV)'; 'g_{GABA_B} (uS)'; ...
        'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
        'm^2h_{T}'; 'm_{\infty}^2h_{\infty,T}'; ...
        'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
        'I_{h} (mA/cm^2)'; 'm_{h}'; 'I_{A} (mA/cm^2)'; ...
        'm_{1,A}'; 'h_{1,A}'; 'm_{2,A}'; 'h_{2,A}'; ...
        'I_{Kir} (mA/cm^2)'; 'm_{\infty,Kir}'; ...
        'I_{NaP} (mA/cm^2)'; 'm_{\infty,NaP}'; 'h_{NaP}'};

figName = fullfile(outFolder, [expStr, '_simulated.png']);
% Count the number of subplots
nSubPlots = numel(yLabelsOverlapped);
% Create figure
figOverlapped = set_figure_properties('AlwaysNew', true, ...
                'FigExpansion', [1, nSubPlots/4]);
figM2h = set_figure_properties('AlwaysNew', true, ...
                'FigExpansion', [1, 1/2]);
figNameM2h = fullfile(outFolder, [expStr, '_simulated_m2h.png']);
save_all_figtypes(figM2h, figNameM2h, figTypes);
handles.handlesOverlapped = handlesOverlapped;
handles.handlesM2h = handlesM2h;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%