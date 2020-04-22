function handles = m3ha_plot_simulated_traces (varargin)
%% Plots simulated traces from single neuron output files
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
%                       'essential'
%                       'somaVoltage'
%                       'allVoltages'
%                       'allTotalCurrents'
%                       'allComponentCurrents'
%                       'allITproperties'
%                       'dend2ITproperties'
%                       'm2h'           - m2h plot
%                       'voltageVsOpd'
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
%                   - 'CompareWithRecorded': whether to compare with recorded
%                                               data when available
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == set in all_files.m
%                   - 'FileNames': paths to simulated data
%                   must be a string array or a cell array of character vectors
%                   default == detected from Directory
%                   - 'SimParamsTable': simulation parameters table
%                   must be a table
%                   default == detected from Directory
%                   - 'Extension': data file extension
%                   must be a string scalar or a character vector
%                   default == 'out'
%                   - 'ColorMap': color map
%                   must be TODO
%                   default == TODO
%                   - 'TimeLimits': limits of time axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == no restrictions
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == same as timeLimits
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
%       cd/compute_derivative_trace.m
%       cd/compute_stats.m
%       cd/compute_total_current.m
%       cd/construct_fullpath.m
%       cd/convert_units.m
%       cd/count_vectors.m
%       cd/decide_on_colormap.m
%       cd/find_zeros.m
%       cd/force_column_cell.m
%       cd/extract_columns.m
%       cd/extract_common_prefix.m
%       cd/extract_elements.m
%       cd/extract_fields.m
%       cd/extract_subvectors.m
%       cd/extract_vars.m
%       cd/find_first_match.m
%       cd/force_matrix.m
%       cd/isemptycell.m
%       cd/load_neuron_outputs.m
%       cd/movingaveragefilter.m
%       cd/m3ha_extract_sweep_name.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_figure05.m
%       cd/parse_peaks.m
%       cd/plot_fitted_traces.m
%       cd/plot_selected.m
%       cd/plot_traces.m
%       cd/plot_window_boundaries.m
%       cd/read_lines_from_file.m
%       cd/set_default_flag.m
%       cd/sscanf_full.m
%       cd/unique_custom.m
%
% Used by:
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_simulate_population.m

% File History:
% 2019-10-14 Created by Adam Lu
% 2019-12-22 Added 'PlotType' as an optional argument
% 2019-12-29 Added 'allVoltages', 'allTotalCurrents', 'allITproperties', 
%               and 'dend2ITproperties' 
% 2019-12-29 Reordered simulated ouptut columns to include ipas
% 2020-01-06 - Now makes the individual plot figure size proportional to the 
%               number of rows and columns
% 2020-01-30 Added 'somaVoltage'
% 2020-02-08 Added m2h difference
% 2020-02-09 Now plots m2h subplots in log scale
% 2020-02-10 Added m2h ratio
% 2020-04-09 The default expStr is just the base name of the common prefix
% 2020-04-09 Now defaults outFolder to common directory of files provided
% 2020-04-12 Removed absolute value from itm2hDiffDend2
% 2020-04-12 Now plots IDX_M2HDIFF_DEND2 in essential
% 2020-04-13 Added 'voltageVsOpd' as a valid plot type
% 2020-04-13 Added 'TimeLimits' as an optional argument
% 2020-04-17 Now plots voltage traces in 'voltageVsOpd' plot
% 2020-04-20 Added 'SimParamsTable' as an optional argument
% 2020-04-20 Now makes sure that simParamsTable, sweepName and fileNames 
%               match each other
% 2020-04-20 Now highlights time limits for voltage traces 
%               in 'voltageVsOpd' plot

%% Hard-coded parameters
validPlotTypes = {'individual', 'residual', 'overlapped', ...
                    'essential', 'somaVoltage', ...
                    'allVoltages', 'allTotalCurrents', ...
                    'allComponentCurrents', 'allITproperties', ...
                    'dend2ITproperties', 'm2h', 'voltageVsOpd'};
validBuildModes = {'', 'active', 'passive'};
validSimModes = {'', 'active', 'passive'};
maxRowsWithOneOnly = 8;
lineWidthParallel = 1;
lineWidthIndividual = 0.5;

% Note: The following must be consistent with m3ha_neuron_run_and_analyze.m
importedSuffix = 'imported_files';
paramsSuffix = 'simulation_parameters';
simOutPathStr = 'outFilePath';

% Note: The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 2000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized


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
compareWithRecordedDefault = true;
directoryDefault = '';          % set in all_files.m
fileNamesDefault = {};
simParamsTableDefault = table.empty;
extensionDefault = 'out';       % 
colorMapDefault = [];
timeLimitsDefault = [];         % set later
xLimitsDefault = [];            % set later
outFolderDefault = '';          % set later
expStrDefault = '';             % set later
tVecsDefault = [];
vVecsRecDefault = [];
iVecsRecDefault = [];
gVecsRecDefault = [];
residualsDefault = [];
lineWidthDefault = [];          % set later

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
addParameter(iP, 'CompareWithRecorded', compareWithRecordedDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fileNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'SimParamsTable', simParamsTableDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColorMap', colorMapDefault);
addParameter(iP, 'TimeLimits', timeLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
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
compareWithRecorded = iP.Results.CompareWithRecorded;
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
simParamsTable = iP.Results.SimParamsTable;
extension = iP.Results.Extension;
colorMap = iP.Results.ColorMap;
timeLimits = iP.Results.TimeLimits;
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
% Determine whether recorded traces needs to be imported
switch plotType
    case {'individual', 'residual', 'overlapped', 'allVoltages'}
        toImportRecorded = set_default_flag([], compareWithRecorded);
    case {'essential', 'somaVoltage',...
            'allTotalCurrents', 'allComponentCurrents', ...
            'allITproperties', 'dend2ITproperties', 'm2h', 'voltageVsOpd'}
        toImportRecorded = false;
end

% Use the present working directory for both inputs and output by default
if isempty(directory) && isempty(fileNames)
    directory = pwd;
end

% Decide on input paths
if isempty(fileNames)
    [~, fileNames] = ...
        all_files('Directory', directory, 'Extension', extension, ...
                    'Keyword', 'sim', 'ForceCellOutput', true);
else
    % Extract common directory
    directory = extract_fileparts(fileNames, 'commondirectory');

    % Make sure they are full paths
    fileNames = construct_fullpath(fileNames, 'Directory', directory);
end

% Set default output directory
if isempty(outFolder)
    outFolder = directory;
end

% Reorder the input paths correctly
fileNames = reorder_simulation_output_files(fileNames);

% Use the common expStr as the experiment string
if isempty(expStr)
    expStr = extract_common_prefix(fileNames);
    expStr = extract_fileparts(expStr, 'base');
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

% Decide on timeLimits
if isempty(timeLimits)
    if strcmp(simMode, 'active')
%        timeLimits = [2800, 4500]; 
        timeLimits = [2800, 4800];
    else
        timeLimits = [timeToStabilize, Inf];
    end
end

% Decide on xLimits
if isempty(xLimits) && ~strcmp(plotType, 'voltageVsOpd')
    xLimits = timeLimits;
end

% Count the number of files
nFiles = numel(fileNames);

% Decide on nRows
nRows = decide_on_nrows(nFiles, simMode, maxRowsWithOneOnly);

% Decide on the color map if not provided
if isempty(colorMap)
    switch plotType
        case {'individual', 'residual'}
            % Decide on the color map for individual and residual plots
            colorMap = decide_on_colormap('r', nRows);
        case {'overlapped', 'essential', 'somaVoltage', ...
                'allVoltages', 'allTotalCurrents', ...
                'allComponentCurrents', 'allITproperties', ...
                'dend2ITproperties', 'm2h', 'voltageVsOpd'}
            colorMap = decide_on_colormap([], 4);
            if nFiles > nRows
                nColumns = ceil(nFiles / nRows);
                nSlots = nColumns * nRows;
                colorMap = reshape(repmat(reshape(colorMap, 1, []), ...
                                    nColumns, 1), nSlots, 3);
            end

            % Make sure the color map matches the number of files
            colorMap = decide_on_colormap(colorMap, nFiles);
        otherwise
            % Use default
    end
end

% Decide on the plot line width
if isempty(lineWidth)
    switch plotType
        case {'individual', 'residual'}
            lineWidth = lineWidthIndividual;
        case {'overlapped', 'essential', 'somaVoltage', ...
                'allVoltages', 'allTotalCurrents', ...
                'allComponentCurrents', 'allITproperties', ...
                'dend2ITproperties', 'm2h', 'voltageVsOpd'}
            lineWidth = lineWidthParallel;
        otherwise
            error('plotType unrecognized!');
    end
end

% Decide on the simulation parameters table
if isempty(simParamsTable)
    if contains(expStr, 'sim')
        % Find the corresponding parameters file
        [~, simParamsPath] = all_files('Directory', directory, 'MaxNum', 1, ...
                                'Suffix', paramsSuffix, 'Extension', 'csv');
    else
        % Find the corresponding parameters file
        [~, simParamsPath] = ...
            all_files('Directory', directory, 'Keyword', expStr, ...
                    'MaxNum', 1, 'Suffix', paramsSuffix, 'Extension', 'csv');
    end

    % Load the simulation parameters table
    simParamsTable = readtable(simParamsPath);

    % Restrict table to match with file names
    if contains(expStr, 'sim')
        % Restrict to the simulation number
        simStr = extract_substrings(expStr, 'RegExp', 'sim[\d]*');
        rowsToUse = sscanf_full(simStr, '%d');
    else
        % Restrict to the output file names
        simOutPaths = simParamsTable.(simOutPathStr);
        rowsToUse = find_first_match(fileNames, simOutPaths);
    end
    simParamsTable = simParamsTable(rowsToUse, :);
end

%% Data
if toImportRecorded
    % Extract sweep names from simulated file names
    sweepNames = m3ha_extract_sweep_name(fileNames);

    % Look for matching recorded sweep names
    if any(isemptycell(sweepNames))
        % Look for the imported files log
        [~, importedPath] = ...
            all_files('Directory', directory, 'Prefix', expStr, ...
                                        'Suffix', importedSuffix, 'MaxNum', 1);

        % Read from the imported file log
        if ~isempty(importedPath)
            % Extract sweep names
            sweepNames = read_lines_from_file(importedPath);
        else
            error('Recorded sweep names unrecognized!'); 
        end
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
end

% Load simulated data
% If recorded data provided (tVecs not empty at this point),
%   interpolate simulated data to match the time points of recorded data
% Note: This is necessary because CVODE (variable time step method) 
%       is applied in NEURON
simData = load_neuron_outputs('FileNames', fileNames, 'tVecs', tVecs, ...
                                'ForceCellOutput', true);

% Extract vectors from simulated data
[tVecs, vVecsSim] = extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM]);

%% Plots
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
    case {'overlapped', 'essential', 'somaVoltage', ...
            'allVoltages', 'allTotalCurrents', ...
            'allComponentCurrents', 'allITproperties', 'dend2ITproperties'}
        handles = m3ha_plot_overlapped_traces(simData, vVecsRec, ...
                                    simParamsTable, plotType, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    case 'm2h'
        handles = m3ha_plot_m2h(simData, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments);
    case 'voltageVsOpd'
        handles = m3ha_plot_voltage_vs_opd(simData, buildMode, ...
                                timeLimits, xLimits, colorMap, lineWidth, ...
                                expStr, expStrForTitle, otherArguments);
    otherwise
        error('plotType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fileNames = reorder_simulation_output_files(fileNames)

% Return if there is only one file name
if ischar(fileNames) || numel(fileNames) == 1
    return
end

% Extract just the file base
fileBases = extract_fileparts(fileNames, 'base');

% Extract the simulation number strings with 'sim'
simStrs = extract_substrings(fileBases, 'Regexp', 'sim[\d]*');

% Extract the simulation numbers
simNums = sscanf_full(simStrs, '%d');

% Sort the numbers
[~, origIndex] = sort(simNums);

% Reorder the file names
fileNames = fileNames(origIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nRows = decide_on_nrows(nFiles, simMode, maxRowsWithOneOnly)
%% Decide on the number of rows

% Decide on the number of rows
if nFiles > 1 && strcmp(simMode, 'active')
    nRows = 4;
elseif nFiles <= 3 && strcmp(simMode, 'passive')
    nRows = 3;
elseif nFiles <= maxRowsWithOneOnly
    nRows = nFiles;
else
    nRows = floor(sqrt(nFiles));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_individual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    simMode, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

% TODO
plotSwpWeightsFlag = false;
visibleStatus = 'on';

%% Preparation
% Decide on figure title
figTitle = sprintf('All traces for Experiment %s', expStrForTitle);

% Decide on the axes to be linked
if strcmp(simMode, 'passive')
    linkAxesOption = 'x';
else
    linkAxesOption = 'xy';
end

% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
[tVecs, vVecsSim, vVecsRec] = ...
    argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
            tVecs, vVecsSim, vVecsRec);

%% Do the job
% Print to standard output
fprintf('Plotting figure of individual voltage traces for %s ...\n', expStr);

% Decide on the figure width and height
nSweeps = count_vectors(vVecsSim);
nRows = 4;
nColumns = ceil(nSweeps / nRows);
figExpansion = [nColumns / 3, nRows / 4];

% Plot the individual traces
figHandle = set_figure_properties('Visible', visibleStatus, ...
                'FigExpansion', figExpansion, 'Name', 'All traces');

% Plot the individual traces
handles = plot_fitted_traces(tVecs, vVecsSim, 'ToAnnotate', false, ...
            'DataToCompare', vVecsRec, 'PlotMode', 'parallel', ...
            'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
            'ColorMap', colorMap, 'XLimits', xLimits, ...
            'LineWidth', lineWidth, 'LinkAxesOption', linkAxesOption, ...
            'FigTitle', figTitle, 'PlotSwpWeightsFlag', plotSwpWeightsFlag, ...
            'FigHandle', figHandle, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_residual_traces(tVecs, vVecsSim, vVecsRec, ...
                                    residuals, xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

% TODO
plotSwpWeightsFlag = false;

%% Preparation
% Calculate voltage residuals (simulated - recorded) if necessary
if isempty(residuals) && ~isempty(vVecsRec)
    residuals = compute_residuals(vVecsSim, vVecsRec);
end

% Decide on figure title
figTitle = sprintf('Residuals for Experiment %s', expStrForTitle);

% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
[tVecs, residuals] = ...
    argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
            tVecs, residuals);

%% Do the job
% Print to standard output
fprintf('Plotting figure of residual traces for %s ...\n', expStr);

% Decide on the figure width and height
nSweeps = count_vectors(residuals);
nRows = 4;
nColumns = ceil(nSweeps / nRows);
figExpansion = [nColumns / 3, nRows / 4];

% Plot the individual traces
figHandle = set_figure_properties('Visible', visibleStatus, ...
                'FigExpansion', figExpansion, 'Name', 'All traces');

% Plot the individual traces
handles = plot_fitted_traces(tVecs, residuals, 'ToAnnotate', false, ...
            'PlotMode', 'residuals', ...
            'SubplotOrder', 'bycolor', 'ColorMode', 'byRow', ...
            'ColorMap', colorMap, 'XLimits', xLimits, ...
            'LineWidth', lineWidth, 'LinkAxesOption', 'xy', ...
            'FigTitle', figTitle, 'PlotSwpWeightsFlag', plotSwpWeightsFlag, ...
            'FigHandle', figHandle, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_overlapped_traces (simData, vVecsRec, ...
                                        simParamsTable, plotType, buildMode, ...
                                        xLimits, colorMap, lineWidth, ...
                                        expStr, expStrForTitle, otherArguments)

%% Hard-coded parameters
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

IPAS_SOMA = 9;
IPAS_DEND1 = 10;
IPAS_DEND2 = 11;

IT_SOMA = 12;
IT_M_SOMA = 13;
IT_MINF_SOMA = 14;
IT_H_SOMA = 15;
IT_HINF_SOMA = 16;
IH_SOMA = 17;
IH_M_SOMA = 18;
IA_SOMA = 19;
IA_M1_SOMA = 20;
IA_H1_SOMA = 21;
IA_M2_SOMA = 22;
IA_H2_SOMA = 23;
IKIR_SOMA = 24;
IKIR_M_SOMA = 25;
INAP_SOMA = 26;
INAP_M_SOMA = 27;
INAP_H_SOMA = 28;

IT_DEND1 = 29;
IT_M_DEND1 = 30;
IT_MINF_DEND1 = 31;
IT_H_DEND1 = 32;
IT_HINF_DEND1 = 33;
IH_DEND1 = 34;
IH_M_DEND1 = 35;
IA_DEND1 = 36;
IA_M1_DEND1 = 37;
IA_H1_DEND1 = 38;
IA_M2_DEND1 = 39;
IA_H2_DEND1 = 40;
IKIR_DEND1 = 41;
IKIR_M_DEND1 = 42;
INAP_DEND1 = 43;
INAP_M_DEND1 = 44;
INAP_H_DEND1 = 45;

IT_DEND2 = 46;
IT_M_DEND2 = 47;
IT_MINF_DEND2 = 48;
IT_H_DEND2 = 49;
IT_HINF_DEND2 = 50;
IH_DEND2 = 51;
IH_M_DEND2 = 52;
IA_DEND2 = 53;
IA_M1_DEND2 = 54;
IA_H1_DEND2 = 55;
IA_M2_DEND2 = 56;
IA_H2_DEND2 = 57;
IKIR_DEND2 = 58;
IKIR_M_DEND2 = 59;
INAP_DEND2 = 60;
INAP_M_DEND2 = 61;
INAP_H_DEND2 = 62;

itm2hDiffLowerLimit = 1e-8;

%% Preparation
% Initialize handles
handles = struct;

% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
if strcmpi(buildMode, 'passive')
    [tVecs, vVecsSim, vVecsDend1, vVecsDend2, iExtSim, ...
            iPasSoma, iPasDend1, iPasDend2] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        DEND1_COL_SIM, DEND2_COL_SIM, IEXT_COL_SIM, ...
                        IPAS_SOMA, IPAS_DEND1, IPAS_DEND2]);
else
    [tVecs, vVecsSim, vVecsDend1, ...
            vVecsDend2, gCmdSimUs, iExtSim, ...
            iPasSoma, iPasDend1, iPasDend2, ...
            itSoma, itmSoma, itminfSoma, ithSoma, ithinfSoma, ...
            ihSoma, ihmSoma, ikirSoma, ikirmSoma, ...
            iaSoma, iam1Soma, iah1Soma, iam2Soma, iah2Soma, ...
            inapSoma, inapmSoma, inaphSoma, ...
            itDend1, itmDend1, itminfDend1, ithDend1, ithinfDend1, ...
            ihDend1, ihmDend1, ikirDend1, ikirmDend1, ...
            iaDend1, iam1Dend1, iah1Dend1, iam2Dend1, iah2Dend1, ...
            inapDend1, inapmDend1, inaphDend1, ...
            itDend2, itmDend2, itminfDend2, ithDend2, ithinfDend2, ...
            ihDend2, ihmDend2, ikirDend2, ikirmDend2, ...
            iaDend2, iam1Dend2, iah1Dend2, iam2Dend2, iah2Dend2, ...
            inapDend2, inapmDend2, inaphDend2] = ...
        extract_columns(simData, ...
            [TIME_COL_SIM, VOLT_COL_SIM, DEND1_COL_SIM, ...
            DEND2_COL_SIM, GGABAB_COL_SIM, IEXT_COL_SIM, ...
            IPAS_SOMA, IPAS_DEND1, IPAS_DEND2, ...
            IT_SOMA, IT_M_SOMA, IT_MINF_SOMA, IT_H_SOMA, IT_HINF_SOMA, ...
            IH_SOMA, IH_M_SOMA, IKIR_SOMA, IKIR_M_SOMA, ...
            IA_SOMA, IA_M1_SOMA, IA_H1_SOMA, IA_M2_SOMA, IA_H2_SOMA, ...
            INAP_SOMA, INAP_M_SOMA, INAP_H_SOMA, ...
            IT_DEND1, IT_M_DEND1, IT_MINF_DEND1, IT_H_DEND1, IT_HINF_DEND1, ...
            IH_DEND1, IH_M_DEND1, IKIR_DEND1, IKIR_M_DEND1, ...
            IA_DEND1, IA_M1_DEND1, IA_H1_DEND1, IA_M2_DEND1, IA_H2_DEND1, ...
            INAP_DEND1, INAP_M_DEND1, INAP_H_DEND1, ...
            IT_DEND2, IT_M_DEND2, IT_MINF_DEND2, IT_H_DEND2, IT_HINF_DEND2, ...
            IH_DEND2, IH_M_DEND2, IKIR_DEND2, IKIR_M_DEND2, ...
            IA_DEND2, IA_M1_DEND2, IA_H1_DEND2, IA_M2_DEND2, IA_H2_DEND2, ...
            INAP_DEND2, INAP_M_DEND2, INAP_H_DEND2]);
end

% Convert the table to a structure array
simParamsStructArray = table2struct(simParamsTable);

% Calculate total currents from current densities
compute_current_across_cells = @(x, y, z) ...
    cellfun(@(a, b, c, d) compute_total_current([a, b, c], 'GeomParams', d), ...
            x, y, z, num2cell(simParamsStructArray), 'UniformOutput', false);
iPasTotal = compute_current_across_cells(iPasSoma, iPasDend1, iPasDend2);
if strcmpi(buildMode, 'active')
    % Compute total and component currents (nA)
    [itTotal, itTotalEachCompartment] = ...
        compute_current_across_cells(itSoma, itDend1, itDend2);
    [ihTotal, ihTotalEachCompartment] = ...
        compute_current_across_cells(ihSoma, ihDend1, ihDend2);
    [iaTotal, iaTotalEachCompartment] = ...
        compute_current_across_cells(iaSoma, iaDend1, iaDend2);
    [ikirTotal, ikirTotalEachCompartment] = ...
        compute_current_across_cells(ikirSoma, ikirDend1, ikirDend2);
    [inapTotal, inapTotalEachCompartment] = ...
        compute_current_across_cells(inapSoma, inapDend1, inapDend2);


    % Extract component currents (nA) for each compartment
    [itTotalSoma, itTotalDend1, itTotalDend2] = ...
        extract_columns(itTotalEachCompartment, 1:3);
    [iaTotalSoma, iaTotalDend1, iaTotalDend2] = ...
        extract_columns(iaTotalEachCompartment, 1:3);

    % Compute the total current
    itaTotal = cellfun(@(a, b) a + b, itTotal, iaTotal, 'UniformOutput', false);
end

% Compute the total intrinsic current
if strcmpi(buildMode, 'passive')
    iIntTotal = iPasTotal;
else
    iIntTotal = cellfun(@(a, b, c, d, e, f) a + b + c + d + e + f, ...
                        iPasTotal, itTotal, ihTotal, ...
                        iaTotal, ikirTotal, inapTotal, ...
                        'UniformOutput', false);
end

% Compute the total current
% TODO: Use combine_traces.m?
iTotal = cellfun(@(a, b) a + b, iExtSim, iIntTotal, 'UniformOutput', false);

% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Extract region of interest and force as a matrix
if strcmpi(buildMode, 'passive')
    [tVecs, vVecsRec, vVecsSim, vVecsDend1, vVecsDend2, ...
            iTotal, iExtSim, iIntTotal, iPasTotal] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, vVecsRec, vVecsSim, vVecsDend1, vVecsDend2, ...
                iTotal, iExtSim, iIntTotal, iPasTotal);
elseif strcmpi(buildMode, 'active')
    [tVecs, vVecsRec, vVecsSim, vVecsDend1, ...
            vVecsDend2, gCmdSimUs, iTotal, iExtSim, iIntTotal, iPasTotal, ...
            itTotal, ihTotal, iaTotal, ikirTotal, inapTotal, itaTotal, ...
            itTotalSoma, itTotalDend1, itTotalDend2, ...
            iaTotalSoma, iaTotalDend1, iaTotalDend2, ...
            itSoma, itmSoma, itminfSoma, ithSoma, ithinfSoma, ...
            ihSoma, ihmSoma, ikirSoma, ikirmSoma, ...
            iaSoma, iam1Soma, iah1Soma, iam2Soma, iah2Soma, ...
            inapSoma, inapmSoma, inaphSoma, ...
            itDend1, itmDend1, itminfDend1, ithDend1, ithinfDend1, ...
            ihDend1, ihmDend1, ikirDend1, ikirmDend1, ...
            iaDend1, iam1Dend1, iah1Dend1, iam2Dend1, iah2Dend1, ...
            inapDend1, inapmDend1, inaphDend1, ...
            itDend2, itmDend2, itminfDend2, ithDend2, ithinfDend2, ...
            ihDend2, ihmDend2, ikirDend2, ikirmDend2, ...
            iaDend2, iam1Dend2, iah1Dend2, iam2Dend2, iah2Dend2, ...
            inapDend2, inapmDend2, inaphDend2] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, vVecsRec, vVecsSim, vVecsDend1, ...
                vVecsDend2, gCmdSimUs, iTotal, iExtSim, iIntTotal, iPasTotal, ...
                itTotal, ihTotal, iaTotal, ikirTotal, inapTotal, itaTotal, ...
                itTotalSoma, itTotalDend1, itTotalDend2, ...
                iaTotalSoma, iaTotalDend1, iaTotalDend2, ...
                itSoma, itmSoma, itminfSoma, ithSoma, ithinfSoma, ...
                ihSoma, ihmSoma, ikirSoma, ikirmSoma, ...
                iaSoma, iam1Soma, iah1Soma, iam2Soma, iah2Soma, ...
                inapSoma, inapmSoma, inaphSoma, ...
                itDend1, itmDend1, itminfDend1, ithDend1, ithinfDend1, ...
                ihDend1, ihmDend1, ikirDend1, ikirmDend1, ...
                iaDend1, iam1Dend1, iah1Dend1, iam2Dend1, iah2Dend1, ...
                inapDend1, inapmDend1, inaphDend1, ...
                itDend2, itmDend2, itminfDend2, ithDend2, ithinfDend2, ...
                ihDend2, ihmDend2, ikirDend2, ikirmDend2, ...
                iaDend2, iam1Dend2, iah1Dend2, iam2Dend2, iah2Dend2, ...
                inapDend2, inapmDend2, inaphDend2);
end

% Convert conductance from uS to nS
gCmdSimNs = convert_units(gCmdSimUs, 'uS', 'nS');

% Compute m2h, minf2hinf and m2h difference
itm2hDend2 = (itmDend2 .^ 2) .* ithDend2;
itminf2hinfDend2 = (itminfDend2 .^ 2) .* ithinfDend2;
itm2hDiffDend2 = itm2hDend2 - itminf2hinfDend2;
itm2hDiffDend2(itm2hDiffDend2 < itm2hDiffLowerLimit) = itm2hDiffLowerLimit;
itm2hAbsDiffDend2 = abs(itm2hDend2 - itminf2hinfDend2);
itm2hRatioDend2 = itm2hDend2 ./ itminf2hinfDend2;

% List all possible items to plot
if strcmpi(buildMode, 'passive')
    vecsAll = {vVecsRec; vVecsSim; vVecsDend1; ...
                vVecsDend2; iTotal; iExtSim; ...
                gCmdSimNs; iIntTotal; iPasTotal};
else
    vecsAll = {vVecsRec; vVecsSim; vVecsDend1; ...
                vVecsDend2; iTotal; iExtSim; ...
                gCmdSimNs; iIntTotal; iPasTotal; ...
                itTotal; ihTotal; iaTotal; ikirTotal; inapTotal; itaTotal; ...
                itTotalSoma; itTotalDend1; itTotalDend2; ...
                iaTotalSoma; iaTotalDend1; iaTotalDend2; ...
                itmSoma; itminfSoma; ithSoma; ithinfSoma; ...
                itmDend1; itminfDend1; ithDend1; ithinfDend1; ...
                itmDend2; itminfDend2; ithDend2; ithinfDend2; ...
                itm2hDend2; itminf2hinfDend2; itm2hAbsDiffDend2; ...
                itm2hDiffDend2; itm2hRatioDend2};
end

% List corresponding labels
if strcmpi(buildMode, 'passive')
    labelsAll = {'V_{rec} (mV)'; 'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                'V_{dend2} (mV)'; 'I_{total} (nA)'; 'I_{stim} (nA)'; ...
                'g_{GABA_B} (uS)'; 'I_{int} (nA)'; 'I_{pas} (nA)'};
else
    labelsAll = {'V_{rec} (mV)'; 'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                'V_{dend2} (mV)'; 'I_{total} (nA)'; 'I_{stim} (nA)'; ...
                'g_{GABA_B} (nS)'; 'I_{int} (nA)'; 'I_{pas} (nA)'; ...
                'I_{T} (nA)'; 'I_{h} (nA)'; 'I_{A} (nA)'; ...
                'I_{Kir} (nA)'; 'I_{NaP} (nA)'; 'I_{T} + I_{A} (nA)'; ...
                'I_{T,soma} (nA)'; 'I_{T,dend1} (nA)'; 'I_{T,dend2} (nA)'; ...
                'I_{A,soma} (nA)'; 'I_{A,dend1} (nA)'; 'I_{A,dend2} (nA)'; ...
                'm_{T,soma}'; 'm_{\infty,T,soma}'; ...
                'h_{T,soma}'; 'h_{\infty,T,soma}'; ...
                'm_{T,dend1}'; 'm_{\infty,T,dend1}'; ...
                'h_{T,dend1}'; 'h_{\infty,T,dend1}'; ...
                'm_{T,dend2}'; 'm_{\infty,T,dend2}'; ...
                'h_{T,dend2}'; 'h_{\infty,T,dend2}'; ...
                'm_{T,dend2}^2h_{T,dend2}'; ...
                'm_{\infty,T,dend2}^2h_{\infty,T,dend2}'; ...
                '|m_{T}^2h_{T} - m_{\infty,T}^2h_{\infty,T}|'; ...
                'm_{T}^2h_{T} - m_{\infty,T}^2h_{\infty,T}'; ...
                'm_{T}^2h_{T} / m_{\infty,T}^2h_{\infty,T}'};
end

% List whether y axis should be log scaled
if strcmpi(buildMode, 'passive')
    yIsLogAll = zeros(9, 1);
else
    yIsLogAll = [zeros(33, 1); ones(5, 1)];
end

% List indices
IDX_VREC = 1;
IDX_VSOMA = 2;
IDX_VDEND1 = 3;
IDX_VDEND2 = 4;
IDX_ITOTAL = 5;
IDX_ISTIM = 6;
IDX_GGABAB = 7;
IDX_IINT = 8;
IDX_IPAS = 9;
IDX_IT = 10;
IDX_IH = 11;
IDX_IA = 12;
IDX_IKIR = 13;
IDX_INAP = 14;
IDX_ITA = 15;
IDX_IT_SOMA = 16;
IDX_IT_DEND1 = 17;
IDX_IT_DEND2 = 18;
IDX_IA_SOMA = 19;
IDX_IA_DEND1 = 20;
IDX_IA_DEND2 = 21;
IDX_MT_SOMA = 22;
IDX_MINFT_SOMA = 23;
IDX_HT_SOMA = 24;
IDX_HINFT_SOMA = 25;
IDX_MT_DEND1 = 26;
IDX_MINFT_DEND1 = 27;
IDX_HT_DEND1 = 28;
IDX_HINFT_DEND1 = 29;
IDX_MT_DEND2 = 30;
IDX_MINFT_DEND2 = 31;
IDX_HT_DEND2 = 32;
IDX_HINFT_DEND2 = 33;
IDX_M2H_DEND2 = 34;
IDX_MINF2HINF_DEND2 = 35;
IDX_M2HABSDIFF_DEND2 = 36;
IDX_M2HDIFF_DEND2 = 37;
IDX_M2HRATIO_DEND2 = 38;

% Error check
if numel(labelsAll) ~= IDX_M2HRATIO_DEND2
    error('Index numbers needs to be updated!');
end

% Select data to plot
if strcmpi(buildMode, 'passive')
    switch plotType
        case 'overlapped'
            indToPlot = IDX_VSOMA:numel(vecsAll);
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'essential'
            indToPlot = [IDX_VSOMA, IDX_VDEND1, IDX_VDEND2, IDX_ISTIM, IDX_IINT];
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'somaVoltage'
            indToPlot = IDX_VSOMA;
        case 'allVoltages'
            indToPlot = IDX_VSOMA:IDX_IINT;
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'allTotalCurrents'
            indToPlot = [IDX_VSOMA, IDX_ITOTAL, IDX_ISTIM, IDX_IINT, IDX_IPAS];
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case {'allComponentCurrents', ...
                'allITproperties', 'dend2ITproperties'}
            fprintf(['No currents or channel properties are ', ...
                        'saved in passive sim mode!\n']);
            return
        otherwise
            error('plotType unrecognized!');
    end
else
    switch plotType
        case 'overlapped'
            indToPlot = IDX_VSOMA:numel(vecsAll);
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'essential'
            indToPlot = [IDX_VSOMA, IDX_GGABAB, IDX_ISTIM, IDX_IT, ...
                            IDX_M2HDIFF_DEND2];
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'somaVoltage'
            indToPlot = [IDX_VSOMA, IDX_M2HDIFF_DEND2];
        case 'allVoltages'
            indToPlot = IDX_VSOMA:IDX_IINT;
            if ~isempty(vVecsRec)
                indToPlot = [IDX_VREC, indToPlot];
            end
        case 'allTotalCurrents'
            indToPlot = [IDX_IINT, IDX_ITA, IDX_IPAS:IDX_INAP];
        case 'allComponentCurrents'
            indToPlot = [IDX_IT, IDX_IT_SOMA:IDX_IT_DEND2, ...
                            IDX_IA, IDX_IA_SOMA:IDX_IA_DEND2];
        case 'allITproperties'
            indToPlot = IDX_MT_SOMA:IDX_M2HRATIO_DEND2;
        case 'dend2ITproperties'
            indToPlot = [IDX_IT_DEND2, IDX_MT_DEND2:IDX_M2HRATIO_DEND2];
        otherwise
            error('plotType unrecognized!');
    end
end

% Extract data to plot
[dataForOverlapped, yLabelsOverlapped, yIsLog] = ...
    argfun(@(x) x(indToPlot), vecsAll, labelsAll, yIsLogAll);

% Construct matching time vectors
tVecsForOverlapped = repmat({tVecs}, size(dataForOverlapped));

% Decide on figure title and file name
figTitle = sprintf('Simulated traces for Experiment %s', expStrForTitle);

%% Plots
% Print to standard output
fprintf('Plotting figure of overlapped traces for %s ...\n', expStr);

% Plot overlapped traces
handles = plot_traces(tVecsForOverlapped, dataForOverlapped, ...
                    'Verbose', false, 'PlotMode', 'parallel', ...
                    'SubplotOrder', 'list', 'ColorMode', 'byTraceInPlot', ...
                    'LegendLocation', 'suppress', ...
                    'ColorMap', colorMap, 'XLimits', xLimits, ...
                    'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                    'YLabel', yLabelsOverlapped, ...
                    'FigTitle', figTitle, 'LineWidth', lineWidth, ...
                    otherArguments);

% Update y axis scale
subPlots = handles.subPlots;
for iAx = 1:numel(subPlots)
    if yIsLog(iAx)
        set(subPlots(iAx), 'YScale', 'log');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_m2h (simData, buildMode, ...
                                    xLimits, colorMap, lineWidth, ...
                                    expStr, expStrForTitle, otherArguments)

%% Hard-coded parameters
% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
% IT_M_SOMA = 13;
% IT_MINF_SOMA = 14;
% IT_H_SOMA = 15;
% IT_HINF_SOMA = 16;
IT_M_DEND2 = 47;
IT_MINF_DEND2 = 48;
IT_H_DEND2 = 49;
IT_HINF_DEND2 = 50;

% Only do this for active mode
if strcmpi(buildMode, 'passive')
    handles = struct;
    return
end

%% Process data
% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
[tVecs, itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim] = ...
    extract_columns(simData, [TIME_COL_SIM, IT_M_DEND2, IT_MINF_DEND2, ...
                    IT_H_DEND2, IT_HINF_DEND2]);
    % extract_columns(simData, [TIME_COL_SIM, IT_M_SOMA, IT_MINF_SOMA, ...
    %                 IT_H_SOMA, IT_HINF_SOMA]);

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
figTitle = sprintf('m2h in dend2 for Experiment %s', expStrForTitle);

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
    plot_traces(tVecs, itminf2hinfVecsSim, 'PlotOnly', true, ...
                'LineStyle', ':', 'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'ColorMap', colorMap, 'XLimits', xLimits, otherArguments);

% set(gca, 'YLim', [1e-6, 1]);

% Set the y axis to be log-scaled
set(gca, 'YScale', 'log');

handles.handlesInstantaneous = handlesInstantaneous;
handles.handlesSteadyState = handlesSteadyState;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = m3ha_plot_voltage_vs_opd (simData, buildMode, ...
                                timeLimits, xLimits, colorMap, lineWidth, ...
                                expStr, expStrForTitle, otherArguments)

%% Hard-coded parameters
% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
IDCLAMP_COL_SIM = 5;
IT_M_DEND2 = 47;
IT_MINF_DEND2 = 48;
IT_H_DEND2 = 49;
IT_HINF_DEND2 = 50;
INAP_DEND2 = 60;

itm2hDiffLowerLimit = 1e-8;
selectedMarkerSize = 6;
stimStartMs = 3000;
barRelValue = 0.95;

% Labels
timeLabel = 'Time (ms)';
voltageLabel = 'Voltage (mV)';
itm2hDiffLabel = 'm_{T}^2h_{T} - m_{\infty,T}^2h_{\infty,T}';
otherVecsLabel = 'I_{NaP} (nA)';
figTitle1 = sprintf('Voltage vs time');
figTitle2 = sprintf('Voltage vs m2hdiff');
figTitle3 = sprintf('m2hdiff vs time');
figTitle4 = sprintf('INaP in dendrite 2 vs time');

% Only do this for active mode
if strcmpi(buildMode, 'passive')
    handles = struct;
    return
end

%% Preparation
itm2hDiffLeftBound = xLimits(1);

% Count the number of traces
nTraces = numel(simData);

%% Process data
% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
[tVecs, vVecsSim, iVecsSim, ...
        itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
        otherVecs] = ...
    extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, IDCLAMP_COL_SIM, ...
                    IT_M_DEND2, IT_MINF_DEND2, IT_H_DEND2, IT_HINF_DEND2, ...
                    INAP_DEND2]);

% Find the indices of the time-axis limit endpoints
endPointsForPlots = find_window_endpoints(timeLimits, tVecs);

% Restrict to those endpoints
[tVecs, vVecsSim, iVecsSim, itmVecsSim, ...
        itminfVecsSim, ithVecsSim, ithinfVecsSim, otherVecs] = ...
    argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
            tVecs, vVecsSim, iVecsSim, itmVecsSim, ...
            itminfVecsSim, ithVecsSim, ithinfVecsSim, otherVecs);

% Compute m2hDiff
itm2h = (itmVecsSim .^ 2) .* ithVecsSim;
itminf2hinf = (itminfVecsSim .^ 2) .* ithinfVecsSim;
itm2hDiff = itm2h - itminf2hinf;
itm2hDiff(itm2hDiff < itm2hDiffLowerLimit) = itm2hDiffLowerLimit;
        
% Find endpoint for just the LTS region
methodNumber = 4;
switch methodNumber
case {1, 3, 4}
    switch methodNumber
    case 1
        % Find all LTS regions
        ltsParams = parse_lts(vVecsSim, 'StimStartMs', stimStartMs, ...
                                'tVec0s', tVecs, 'iVec0s', iVecsSim, ...
                                'Verbose', false);
        idxPeakStart = ltsParams.idxPeakStart;
        idxPeakEnd = ltsParams.idxPeakEnd;
    case 3
        % Parse maximum peak from itm2hDiff, 
        %   using itm2hDiffLeftBound as the peak lower bound
        peakParams = ...
            vecfun(@(x) parse_peaks(x, 'ParseMode', 'max', ...
                                    'PeakLowerBound', itm2hDiffLeftBound), ...
                    itm2hDiff, 'UniformOutput', true);

        % Extract index peak starts and ends
        [idxPeakStart, idxPeakEnd] = ...
            argfun(@(x) extract_fields(peakParams, x), ...
                    'idxPeakStart', 'idxPeakEnd');
    case 4
        % Parse maximum peak from itm2hDiff, 
        %   using itm2hDiffLeftBound as the peak lower bound
        peakParams = ...
            vecfun(@(x) parse_peaks(x, 'ParseMode', 'maxOfAll', ...
                                    'PeakLowerBound', itm2hDiffLeftBound), ...
                    itm2hDiff, 'UniformOutput', true);

        % Extract index peak starts and ends
        [idxPeakStart, idxPeakEnd] = ...
            argfun(@(x) extract_fields(peakParams, x), ...
                    'idxPeakStart', 'idxPeakEnd');
    end

    endPointsPeak = transpose([idxPeakStart, idxPeakEnd]);
case 2
    % Extract peak endpoints 
    endPointsPeak = vecfun(@find_lts_endpoints, vVecsSim, ...
                            'UniformOutput', false);
    idxPeakStart = transpose(endPointsPeak(1, :));
    idxPeakEnd = transpose(endPointsPeak(2, :));
end

% Restrict to just the LTS region
[tVecsLts, vVecsLts, itm2hDiffLts, otherVecsLts] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsPeak), ...
            tVecs, vVecsSim, itm2hDiff, otherVecs);

% Restrict to just the pre-LTS region
[tVecsPreLts, vVecsPreLts, itm2hDiffPreLts, otherVecsPreLts] = ...
    argfun(@(x) extract_subvectors(x, 'IndexEnd', idxPeakStart - 1), ...
            tVecs, vVecsSim, itm2hDiff, otherVecs);

% Restrict to just the post-LTS region
[tVecsPostLts, vVecsPostLts, itm2hDiffPostLts, otherVecsPostLts] = ...
    argfun(@(x) extract_subvectors(x, 'IndexStart', idxPeakEnd + 1), ...
            tVecs, vVecsSim, itm2hDiff, otherVecs);

% Extract the LTS start and end times
[timePeakStart, timePeakEnd] = ...
    argfun(@(x) extract_elements(tVecs, 'specific', 'Index', x), ...
            idxPeakStart, idxPeakEnd);

% Find the minimum and maximum times for the LTS regions
minTimePeakStart = min(timePeakStart);
maxTimePeakStart = max(timePeakStart);
minTimePeakEnd = min(timePeakEnd);
maxTimePeakEnd = max(timePeakEnd);

% Find the LTS windows
ltsWindowFlank1 = [minTimePeakStart, maxTimePeakStart];
ltsWindow = [maxTimePeakStart, minTimePeakEnd];
ltsWindowFlank2 = [minTimePeakEnd, maxTimePeakEnd];

% Find inflection points
indInflection = ...
    find_inflection_points_voltage_vs_opd(tVecsLts, vVecsLts, itm2hDiffLts, ...
                                            itm2hDiffLowerLimit);

%% Plots
% Print to standard output
fprintf('Plotting figure of voltage vs m2hdiff for %s ...\n', expStr);

% Create subplots
[fig, ax] = create_subplots(2, 2);

% Link voltage axes
linkaxes(ax([1, 2]), 'y');

% Link time axes
linkaxes(ax([1, 3, 4]), 'x');

% Create same color map but faded
colorMapFaded = decide_on_colormap(colorMap, nTraces, 'FadePercentage', 30);
barColorMap = decide_on_colormap('DarkGreen', 1);
flankColorMap = decide_on_colormap(barColorMap, 'OriginalNColors', true, ...
                                    'FadePercentage', 30);

% Plot voltage vs time
subplot(ax(1))
handles.tracesPre1 = ...
    plot_traces(tVecsPreLts, vVecsPreLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.tracesPost1 = ...
    plot_traces(tVecsPostLts, vVecsPostLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.traces1 = ...
    plot_traces(tVecsLts, vVecsLts, 'XLimits', timeLimits, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'LegendLocation', 'suppress', 'ColorMap', colorMap, ...
                'XLabel', timeLabel, 'YLabel', voltageLabel, ...
                'FigTitle', figTitle1);

% Plot horizontal bars
handles.bar1 = plot_window_boundaries(ltsWindowFlank1, ...
                                        'BoundaryType', 'horizontalBars', ...
                                        'ColorMap', flankColorMap, ...
                                        'BarRelValue', barRelValue);
handles.bar2 = plot_window_boundaries(ltsWindow, ...
                                        'BoundaryType', 'horizontalBars', ...
                                        'ColorMap', barColorMap, ...
                                        'BarRelValue', barRelValue);
handles.bar3 = plot_window_boundaries(ltsWindowFlank2, ...
                                        'BoundaryType', 'horizontalBars', ...
                                        'ColorMap', flankColorMap, ...
                                        'BarRelValue', barRelValue);

% Plot markers for inflection points
handles.selected1 = ...
    plot_selected(tVecsLts, vVecsLts, indInflection, ...
                'ColorMap', colorMap, 'LineWidth', lineWidth, ...
                'Marker', 'o', 'MarkerSize', selectedMarkerSize);

% Plot voltage vs m2hdiff
subplot(ax(2))
handles.traces2 = ...
    plot_traces(itm2hDiffLts, vVecsLts, 'XLimits', xLimits, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'LegendLocation', 'suppress', 'ColorMap', colorMap, ...
                'XLabel', itm2hDiffLabel, 'YLabel', voltageLabel, ...
                'FigTitle', figTitle2, otherArguments);

% Plot markers for inflection points
handles.selected2 = ...
    plot_selected(itm2hDiffLts, vVecsLts, indInflection, ...
                'ColorMap', colorMap, 'LineWidth', lineWidth, ...
                'Marker', 'o', 'MarkerSize', selectedMarkerSize);

% Set the x axis to be log-scaled
set(ax(2), 'XScale', 'log');

% Plot m2hdiff vs time
subplot(ax(3))
handles.tracesPre3 = ...
    plot_traces(tVecsPreLts, itm2hDiffPreLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.tracesPost3 = ...
    plot_traces(tVecsPostLts, itm2hDiffPostLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.traces3 = ...
    plot_traces(tVecsLts, itm2hDiffLts, 'XLimits', timeLimits, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'LegendLocation', 'suppress', 'ColorMap', colorMap, ...
                'XLabel', timeLabel, 'YLabel', itm2hDiffLabel, ...
                'FigTitle', figTitle3);

% Plot markers for inflection points
handles.selected3 = ...
    plot_selected(tVecsLts, itm2hDiffLts, indInflection, ...
                'ColorMap', colorMap, 'LineWidth', lineWidth, ...
                'Marker', 'o', 'MarkerSize', selectedMarkerSize);

% Set the y axis to be log-scaled
set(ax(3), 'YScale', 'log');
set(ax(3), 'YLim', get(ax(2), 'XLim'));

% Plot otherVecs vs time
subplot(ax(4))
handles.tracesPre4 = ...
    plot_traces(tVecsPreLts, otherVecsPreLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.tracesPost4 = ...
    plot_traces(tVecsPostLts, otherVecsPostLts, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'PlotOnly', true, 'ColorMap', colorMapFaded);
handles.traces4 = ...
    plot_traces(tVecsLts, otherVecsLts, 'XLimits', timeLimits, ...
                'Marker', '.', 'LineStyle', 'none', ...
                'LineWidth', lineWidth, ...
                'Verbose', false, 'PlotMode', 'overlapped', ...
                'LegendLocation', 'suppress', 'ColorMap', colorMap, ...
                'XLabel', timeLabel, 'YLabel', otherVecsLabel, ...
                'FigTitle', figTitle4);

% Plot markers for inflection points
handles.selected4 = ...
    plot_selected(tVecsLts, otherVecsLts, indInflection, ...
                'ColorMap', colorMap, 'LineWidth', lineWidth, ...
                'Marker', 'o', 'MarkerSize', selectedMarkerSize);

% Create overarching title
suptitle(expStrForTitle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indInflection = ...
                find_inflection_points_voltage_vs_opd (tVecs, vVecsSim, ...
                                                itm2hDiff, itm2hDiffLowerLimit)

% Hard-coded parameters
methodNumber = 11;
filtWidthMs = 10;            % in ms

% Compute sampling interval
siMs = compute_sampling_interval(tVecs);

% Compute dv/dt
[dvdtVecs, t1Vecs] = compute_derivative_trace(vVecsSim, tVecs);

% Smooth dv/dt over filtWidthMs
if ismember(methodNumber, [2, 3, 4, 5])
    dvdtVecs = movingaveragefilter(dvdtVecs, filtWidthMs, siMs);
end

% Compute d2v/dt2
if ismember(methodNumber, [3, 4, 5])
    d2vdt2Vecs = compute_derivative_trace(dvdtVecs, t1Vecs);
end

% Compute log10(itm2hDiff) or x
if ismember(methodNumber, [6, 7, 8, 9, 10, 11])
    % Force as column cell arrays
    itm2hDiffCell = force_column_cell(itm2hDiff);

    % Compute the logarithm
    logItm2hDiff = cellfun(@log10, itm2hDiffCell, 'UniformOutput', false);

    % Smooth logItm2hDiff over filtWidthMs
    % logItm2hDiffSmooth = movingaveragefilter(logItm2hDiff, filtWidthMs, siMs);
end

% Compute dv/dx in the phase plot
if ismember(methodNumber, [6, 7])
    dvdxVecs = compute_derivative_trace(vVecsSim, logItm2hDiff);
end

% Compute dx/dt
if ismember(methodNumber, [8, 9, 10, 11])
    [dxdtVecs, t1Vecs] = compute_derivative_trace(logItm2hDiff, tVecs);
end

% Compute d2x/dt2
if ismember(methodNumber, [9, 10, 11])
    % Smooth dx/dt over filtWidthMs
    dxdtVecs = movingaveragefilter(dxdtVecs, filtWidthMs, siMs);

    % Compute d2x/dt2
    d2xdt2Vecs = compute_derivative_trace(dxdtVecs, t1Vecs);
end

% Compute index of maximum concavity for logItm2hDiff 
%   before the maximum of logItm2hDiff
if ismember(methodNumber, [11, 12])
    % Extract the index of maximum value for logItm2hDiff 
    [~, ind2MaxValue] = extract_elements(logItm2hDiff, 'max');

    % Restrict to the part of logItm2hDiff before the maximum is reached
    d2xdt2VecsLeft = extract_subvectors(d2xdt2Vecs, 'IndexEnd', ind2MaxValue);

    % Extract the index of maximum concavity for logItm2hDiff 
    %   before the maximum of logItm2hDiff
    [~, ind2MaxConcavityBeforeMax] = extract_elements(d2xdt2VecsLeft, 'max');
end

switch methodNumber
case 1
    % Method 1: maximum of dv/dt
    % Extract the point of maximum slope for each voltage vector
    [~, ind1MaxSlope] = extract_elements(dvdtVecs, 'max');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind1MaxSlope + 1);
case 2
    % Method 2: local minima and maxima of dv/dt
    % Force as a cell array of column vectors
    dvdtVecsCell = force_column_cell(dvdtVecs);

    % Prominence must be at least 1/10th of the range
    minPeakProminence = compute_stats(dvdtVecsCell, 'range') / 10;

    % Find the peaks and troughs for each slope vector
    peakParams = cellfun(@(x, y) parse_peaks(x, 'ParseMode', 'all', ...
                                'MinPeakProminence', y), ...
                        dvdtVecsCell, num2cell(minPeakProminence), ...
                        'UniformOutput', false);

    % Extract the indices
    %   Note: must be consistent with parse_peaks.m
    [ind1Peaks, ind1PeakStarts, ind1PeakEnds] = ...
        argfun(@(x) extract_vars(peakParams, x), ...
                'idxPeak', 'idxPeakStart', 'idxPeakEnd');

    % Combine peaks and troughs
    ind1PeakTroughs = cellfun(@(a, b, c) unique_custom([a; b; c], 'sorted', ...
                                                        'IgnoreNan', true), ...
                                ind1Peaks, ind1PeakStarts, ind1PeakEnds, ...
                                'UniformOutput', false);

    % Find the inflection points in each voltage vector
    indInflection = cellfun(@(x) x + 1, ind1PeakTroughs, ...
                            'UniformOutput', false);
case 3
    % Method 3: zeros of d2v/dt2
    % Find the indices with d2v/dt2 closest to zero
    ind2Inflection = find_zeros(d2vdt2Vecs);

    % Force as column cell arrays
    ind2Inflection = force_column_cell(ind2Inflection);

    % Find the corresponding indices in the voltage and m2hDiff vectors
    indInflection = cellfun(@(x) x + 1, ind2Inflection, 'UniformOutput', false);
case 4
    % Method 4: maximum of d2v/dt2
    % Extract the point of maximum concavity for each voltage vector
    [~, ind2MaxConcavity] = extract_elements(d2vdt2Vecs, 'max');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind2MaxConcavity + 1);
case 5
    % Method 5: half-maximum of d2v/dt2
    % Extract the maximum concavity for each voltage vector
    valMaxSlope = extract_elements(d2vdt2Vecs, 'max');

    % Force as column cell arrays
    d2vdt2VecsCell = force_column_cell(d2vdt2Vecs);

    % Find the indices with d2v/dt2 closest to half-maximum
    ind2Inflection = cellfun(@(x, y) find_zeros(x - y/2), ...
                            d2vdt2VecsCell, num2cell(valMaxSlope), ...
                            'UniformOutput', false);

    % Force as column cell arrays
    ind2Inflection = force_column_cell(ind2Inflection);

    % Find the corresponding indices in the voltage and m2hDiff vectors
    indInflection = cellfun(@(x) x + 1, ind2Inflection, 'UniformOutput', false);
case 6
    % Method 6: maximum of dv/dx
    % Extract the point of maximum slope for the phase plot
    [~, ind1MaxSlope] = extract_elements(dvdxVecs, 'max');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind1MaxSlope + 1);
case 7
    % Method 7: first local maximum of dv/dx
    % Smooth dv/dx over 3 sample points
    dvdxVecs = movingaveragefilter(dvdxVecs, 3);

    % Force as a column cell array
    dvdxVecsCell = force_column_cell(dvdxVecs);

    % Prominence must be at least 1/10th of the range
    minPeakProminence = compute_stats(dvdxVecsCell, 'range') / 50;

    % Find the peaks and troughs for each slope vector
    peakParams = cellfun(@(x, y) parse_peaks(x, 'ParseMode', 'first', ...
                                'MinPeakProminence', y), ...
                        dvdxVecsCell, num2cell(minPeakProminence));

    % Extract the first local maximum of the phase plot
    ind1FirstMaxSlope = extract_fields(peakParams, 'idxPeak');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind1FirstMaxSlope + 1);
case 8
    % Method 8: maximum of dx/dt
    % Extract the point of maximum slope for logItm2hDiff
    [~, ind1MaxSlope] = extract_elements(dxdtVecs, 'max');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind1MaxSlope + 1);
case 9
    % Method 9: maximum of d2x/dt2
    % Extract the point of maximum concavity for logItm2hDiff
    [~, ind2MaxConcavity] = extract_elements(d2xdt2Vecs, 'max');

    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind2MaxConcavity + 1);
case 10
    % Method 10: last zero before the maximum of d2x/dt2
    % Extract the point of maximum concavity for logItm2hDiff
    [~, ind2MaxConcavity] = extract_elements(d2xdt2Vecs, 'max');

    % Find the indices with d2x/dt2 closest to zero
    ind2Zeros = find_zeros(d2xdt2Vecs);

    % Find the index with d2x/dt2 closest to zero before the last maximum
    ind2LastZeroBeforeMax = ...
        cellfun(@(x, y) x(find(x < y, 1, 'last')), ...
                ind2Zeros, num2cell(ind2MaxConcavity), 'UniformOutput', false);

    % Find the corresponding indices in the voltage and m2hDiff vectors
    indInflection = cellfun(@(x) x + 1, ind2LastZeroBeforeMax, ...
                            'UniformOutput', false);
case 11
    % Method 11: maximum of d2x/dt2 before maximum of x
    % Find the inflection point in each voltage vector
    indInflection = num2cell(ind2MaxConcavityBeforeMax + 1);
case 12
    % Method 11: last zero before the maximum of d2x/dt2 before maximum of x
end

% Force as column cell arrays
% itm2hDiff = force_column_cell(itm2hDiff);

% Remove values that are not in view
% indInflection = cellfun(@(a, b) setdiff(a, find(b == itm2hDiffLowerLimit)), ...
%                         indInflection, itm2hDiff, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecs = prepare_for_plotting(vecs, endPointsForPlots)
%% Prepare vectors for plotting

% Restrict vectors to xLimits to save time on plotting
vecs = extract_subvectors(vecs, 'Endpoints', endPointsForPlots);

% Combine vectors into matrices
vecs = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endPoints = find_lts_endpoints(vVec)
% TODO: Merge with parse_lts.m?

% Set a minimum peak prominence
minPeakProminence = range(vVec) / 20;

% Find all the troughs
[negAmpTroughs, indTroughs] = ...
    findpeaks(-vVec, 'MinPeakProminence', minPeakProminence);
ampTroughs = -negAmpTroughs;

% Find the index for peak start
if ~isempty(ampTroughs)
    % Find the trough with smallest amplitude
    [~, iTrough1] = min(ampTroughs);

    % Index peak start
    idxPeakStart = indTroughs(iTrough1);

    % Remove this trough
    ampTroughs(iTrough1) = [];
    indTroughs(iTrough1) = [];
else
    idxPeakStart = 1;
end

% Find the index for peak end
if ~isempty(ampTroughs)
    [~, iTrough2] = min(ampTroughs);

    % Index peak end
    idxPeakEnd = indTroughs(iTrough2);
else
    idxPeakEnd = numel(vVec);
end

% Return endpoints
endPoints = [idxPeakStart; idxPeakEnd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

dataForOverlapped = {vVecsSim; gCmdSim; iExtSim; ...
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

figHandle = set_figure_properties('Visible', visibleStatus, ...
                'AlwaysNew', true, 'FigExpansion', figExpansion, ...
                'Name', 'All traces');

% Extract the simulation numbers (still in string form)
simNumStrs = extractAfter(simStrs, 'sim');

% Convert numeric strings to numbers
simNums = str2double(simNumStrs);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
