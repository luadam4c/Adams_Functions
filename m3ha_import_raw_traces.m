function [data, sweepInfo, dataAll] = m3ha_import_raw_traces (fileNames, varargin)
%% Imports raw traces from .mat files in the m3ha format
% Usage: [data, sweepInfo, dataAll] = m3ha_import_raw_traces (fileNames, varargin)
% Explanation:
%       TODO
%
% Examples:
%       [data, sweepInfo, dataAll] = m3ha_import_raw_traces('D091710_0000_20');
%
% Outputs:
%       data        - imported raw traces; 
%                       Do this to extract:
%                       [tVecs, vVecs, iVecs, gVecs] = extract_columns(data, 1:4);
%       sweepInfo   - sweep info table, with fields:
%                       fileBases
%                       fileNames
%                       siMs
%                       currentPulseAmplitude
%                       vHold
%                       holdPotential
%                       holdCurrent
%                       baseNoise
%                       holdCurrentNoise
%                       sweepWeights
%                       gababAmp
%                       gababTrise
%                       gababTfallFast
%                       gababTfallSlow
%                       gababWeight
%       dataAll     - all imported raw traces
%
% Arguments:
%       fileNames   - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'CreateLog': whether to create a log file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ImportMode': import mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'full'      - full trace without padding
%                       'passive'   - current pulse response with padding
%                       'active'    - IPSC response with padding
%                   default == 'full'
%                   - 'ToParsePulse': whether to parse pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set based on importMode
%                   - 'ToMedianFilter': whether to median filter data
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set based on importMode
%                   - 'ToResample': whether to resample data
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set based on importMode
%                   - 'ToCorrectDcSteps': whether to correct unbalanced
%                                            bridges in the pulse responses 
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true (only has effect if toParsePulse is true)
%                   - 'ToAverageByVhold': whether to average pulse responses 
%                                           according to VHold
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true (only has effect if toParsePulse is true)
%                   - 'ToBootstrapByVhold': whether to bootstrap average 
%                                           pulse responses within each VHold
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true (only has effect if toParsePulse is true)
%                   - 'Directory': a directory for .mat files
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == fullfile(dataDir, matFilesDirName);
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'TimeToPad': time to pad in the beginning of the trace
%                                   (same units as the sampling interval 
%                                       in the time data)
%                   must be a nonnegative scalar
%                   default == set based on importMode
%                   - 'ResponseWindow': window(s) of response vector
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == set based on importMode
%                   - 'StimStartWindow': window(s) of stimulation start
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == set based on importMode
%                   - 'EpasEstimate': estimate for the reversal potential (mV)
%                   must be a numeric scalar
%                   default == []
%                   - 'RinEstimate': estimate for the input resistance (MOhm)
%                   must be a nonnegative scalar
%                   default == []
%                   - 'SweepInfoAll': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == m3ha_load_sweep_info
%                   - 'InitialSlopesPath': path to the initial slopes .mat file
%                                   that contains the following variables:
%                                       'filenamesSorted'
%                                       'iThreshold1Balanced'
%                                       'iThreshold2Balanced'
%                   must be a string scalar or a character vector
%                   default == ~/m3ha/data_dclamp/take4/initial_slopes_nSamplesForPlot_2_threeStdMainComponent.mat
%                   - 'PassiveParamsPath': path to passive parameters .xlsx file
%                                   that contains the following columns
%                                       'epas'
%                                       'Rin'
%                   must be a string scalar or a character vector
%                   default == ~/m3ha/data_dclamp/take4/dclampPassiveParams_byCells_tofit.xlsx
%
% Requires:
%       cd/apply_or_return.m
%       cd/argfun.m
%       cd/compute_combined_data.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sampling_interval.m
%       cd/construct_and_check_fullpath.m
%       cd/convert_to_samples.m
%       cd/correct_unbalanced_bridge.m
%       cd/create_time_vectors.m
%       cd/extract_columns.m
%       cd/extract_common_prefix.m
%       cd/extract_elements.m
%       cd/extract_fileparts.m
%       cd/extract_subvectors.m
%       cd/extract_vars.m
%       cd/find_in_strings.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/force_string_end.m
%       cd/print_cellstr.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_locate_homedir.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/set_default_flag.m
%
% Used by:
%       cd/m3ha_plot_example_jitter.m.m
%       cd/m3ha_simulate_population.m
%       cd/m3ha_xolotl_plot.m
%       cd/singleneuronfitting47.m and later versions

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-12-21 Changed tabs to spaces
% 2018-05-16 Now passes in vrow
% 2018-05-18 Changed all variables to camelback case
% 2018-05-18 Averaged the current pulse responses according to vHold
% 2018-05-21 Added outparams.baseNoiseIpscr and outparams.baseNoiseCpr
% 2018-06-20 Now uses initialSlopes to take 
%               out cpr traces with out-of-balance bridges
% 2018-06-21 Fixed bugs
% 2018-07-09 Added nSwpsCpr as output
% 2018-07-09 Now uses compute_rms_error() instead of rms_Gaussian()
% 2018-07-31 Use correct_unbalanced_bridge.m to correct for out-of-balance bridges
% 2018-08-09 Now computes sweep weights here
% 2018-08-10 baseNoiseIpscr and baseNoiseCpr are now column vectors
% 2018-09-12 Added outparams.ToAverageByVhold
% 2018-11-15 Moved to Adams_Functions
% 2018-11-28 Now pads vvecsIpscr with NaN instead of with holdPotentialIpscr
% 2019-01-01 Made all arguments optional except fileNames
%               and consolidated the different types of responses
% 2019-01-12 Now uses compute_combined_data.m
% 2019-10-15 Fixed conversion of current pulse amplitude from swpInfo
% 2019-11-13 Added siMs to swpInfo
% 2019-11-14 Added 'ImportMode' as an optional parameter
% 2019-12-04 Now imports the GABA-B IPSC parameters used in dynamic clamp
% TODO: Use m3ha_correct_unbalanced_bridge.m to not correct for all sweeps

%% Hard-coded constants
NS_PER_US = 1000;
PA_PER_NA = 1000;

%% Hard-coded parameters
validImportModes = {'full', 'active', 'passive'};

% Parameters to be consistent with find_passive_params.m
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations
dataDirName = fullfile('data_dclamp', 'take4');
matFilesDirName = 'matfiles';
initialSlopesFileName = 'initial_slopes_nSamplesForPlot_2_threeStdMainComponent.mat';
passiveFileName = 'dclampPassiveParams_byCells_tofit.xlsx';

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 2000;     % time to make sure initial value of simulation is stabilized

% The following must be consistent with dclampDataExtractor.m
ipscTimeOrig = 1000;                % time of IPSC application (ms), original
cpStartWindowOrig = [95, 105];      % window in which the current pulse start would lie (ms) 
                                    %   (Supposed to be 100 ms but there will be offset)

% The following must be consistent with both dclampDataExtractor.m & singleneuron4compgabab.hoc
cprWinOrig = [0, 360];              % window in which the current pulse response would lie (ms), original
ipscrWinOrig = [0, 8000];           % window in which the IPSC response would lie (ms), original

% Parameters used for data reorganization
%   Note: should be consistent with m3ha_resave_sweeps.m
rsims = 1;  % resampling interval in ms (1 kHz)
mfw1 = 2.5; % width in ms for the median filter for PClamp noise (conductance traces)
mfw2 = 10;  % width in ms for the median filter for corrupted data (current traces)
mfw3 = 30;  % width in ms for the median filter for spikes (voltage traces)

%% Default values for optional arguments
importModeDefault = 'full';     % import full trace by default
verboseDefault = true;          % print to standard output by default
createLogDefault = false;       % don't create log file by default
toParsePulseDefault = [];       % set later
toMedianFilterDefault = [];     % set later
toResampleDefault = [];         % set later
toCorrectDcStepsDefault = [];   % set later
toAverageByVholdDefault = [];   % set later
toBootstrapByVholdDefault = []; % set later
directoryDefault = '';
outFolderDefault = pwd;
timeToPadDefault = [];          % set later
responseWindowDefault = [];     % set later
stimStartWindowDefault = [];    % set later
epasEstimateDefault = [];
rinEstimateDefault = [];
swpInfoDefault = [];
initialSlopesPathDefault = '';
passiveParamsPathDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'fileNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['fileNames must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CreateLog', createLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ImportMode', importModeDefault, ...
    @(x) any(validatestring(x, validImportModes)));
addParameter(iP, 'ToParsePulse', toParsePulseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToMedianFilter', toMedianFilterDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToResample', toResampleDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToCorrectDcSteps', toCorrectDcStepsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToAverageByVhold', toAverageByVholdDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToBootstrapByVhold', toBootstrapByVholdDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeToPad', timeToPadDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'ResponseWindow', responseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['ResponseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'StimStartWindow', stimStartWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['StimStartWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EpasEstimate', epasEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RinEstimate', rinEstimateDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'SweepInfoAll', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'InitialSlopesPath', initialSlopesPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PassiveParamsPath', passiveParamsPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the input Parser
parse(iP, fileNames, varargin{:});
verbose = iP.Results.Verbose;
createLog = iP.Results.CreateLog;
importMode = validatestring(iP.Results.ImportMode, validImportModes);
toParsePulse = iP.Results.ToParsePulse;
toMedianFilter = iP.Results.ToMedianFilter;
toResample = iP.Results.ToResample;
toCorrectDcSteps = iP.Results.ToCorrectDcSteps;
toAverageByVhold = iP.Results.ToAverageByVhold;
toBootstrapByVhold = iP.Results.ToBootstrapByVhold;
matFilesDir = iP.Results.Directory;
outFolder = iP.Results.OutFolder;
timeToPad = iP.Results.TimeToPad;
responseWindowOrig = iP.Results.ResponseWindow;
stimStartWindowOrig = iP.Results.StimStartWindow;
epasEstimate = iP.Results.EpasEstimate;
RinEstimate = iP.Results.RinEstimate;
swpInfoAll = iP.Results.SweepInfoAll;
initialSlopesPath = iP.Results.InitialSlopesPath;
passiveParamsPath = iP.Results.PassiveParamsPath;

%% Prepare
% Print message
fprintf('Preparing for import ... \n');

% Decide whether to parse pulse
if isempty(toParsePulse)
    switch importMode
        case 'passive'
            toParsePulse = true;
        case {'full', 'active'}
            toParsePulse = false;
        otherwise
            error('importMode unrecognized!');
    end
end

% Decide whether to median filter traces
if isempty(toMedianFilter)
    switch importMode
        case 'active'
            toMedianFilter = true;
        case {'full', 'passive'}
            toMedianFilter = false;
        otherwise
            error('importMode unrecognized!');
    end
end

% Decide whether to resample traces
if isempty(toResample)
    switch importMode
        case 'active'
            toResample = true;
        case {'full', 'passive'}
            toResample = false;
        otherwise
            error('importMode unrecognized!');
    end
end

% Decide whether to modify current pulse responses
[toCorrectDcSteps, toAverageByVhold, toBootstrapByVhold] = ...
    argfun(@(x) set_default_flag(x, toParsePulse), ...
            toCorrectDcSteps, toAverageByVhold, toBootstrapByVhold);

% Decide on the amount of time to pad in ms
if isempty(timeToPad)
    switch importMode
        case 'full'
            timeToPad = 0;
        case {'active', 'passive'}
            timeToPad = timeToStabilize;
        otherwise
            error('importMode unrecognized!');
    end
end

% Decide on the response window in ms
if isempty(responseWindowOrig)
    switch importMode
        case 'full'
            responseWindowOrig = [0, Inf];
        case 'active'
            responseWindowOrig = ipscrWinOrig;
        case 'passive'
            responseWindowOrig = cprWinOrig;
        otherwise
            error('importMode unrecognized!');
    end
end

% Decide on the stimulation start window in ms
if isempty(stimStartWindowOrig)
    switch importMode
        case 'full'
            stimStartWindowOrig = ipscTimeOrig;
        case 'active'
            stimStartWindowOrig = ipscTimeOrig;
        case 'passive'
            stimStartWindowOrig = cpStartWindowOrig;
        otherwise
            error('importMode unrecognized!');
    end
end

% Locate the m3ha home directory if needed
if isempty(matFilesDir) || isempty(initialSlopesPath) || isempty(swpInfoAll)
    homeDirectory = m3ha_locate_homedir;
end

% Locate the data directory if needed
if isempty(matFilesDir) || isempty(initialSlopesPath)
    dataDir = fullfile(homeDirectory, dataDirName);
end

% Locate the matfiles directory if not provided
if isempty(matFilesDir)
    matFilesDir = fullfile(dataDir, matFilesDirName);
end

% Locate the initial slopes file if not provided
if isempty(initialSlopesPath)
    initialSlopesPath = fullfile(dataDir, initialSlopesFileName);
end

% Load sweep information if not provided
%   Note: the file names are read in as row names
if isempty(swpInfoAll)
    swpInfoAll = m3ha_load_sweep_info('Directory', dataDir);
end

% Locate the passive parameters file
if isempty(passiveParamsPath)
    passiveParamsPath = fullfile(dataDir, passiveFileName);
end

% Make sure fileNames is a column cell array
fileNames = force_column_cell(fileNames);

% Make sure fileNames all end in .mat
fileNames = force_string_end(fileNames, '.mat');

% Extract file bases
fileBases = extract_fileparts(fileNames, 'base');

% Count the total number of sweeps to import
nSwps = numel(fileBases);

% Construct full paths to matfiles and check if they exist
[filePaths, pathExists] = ...
    construct_and_check_fullpath(fileNames, 'Directory', matFilesDir, ...
                                'Verbose', verbose, 'Extension', 'mat');

% Get the cell name
cellName = m3ha_extract_cell_name(fileBases);

% Extract estimated epas if available
if isempty(epasEstimate) && isfile(passiveParamsPath) && ~isempty(cellName)
    epasEstimate = extract_vars(passiveParamsPath, 'epasEstimate', ...
                                'RowConditions', {'cellName', cellName});
end

% Extract estimated Rin if available
if isempty(RinEstimate) && isfile(passiveParamsPath) && ~isempty(cellName)
    RinEstimate = extract_vars(passiveParamsPath, 'Rinput', ...
                                'RowConditions', {'cellName', cellName});

end

% If a matfile does not exist, return
if ~all(pathExists)
    fprintf('One or more matfiles does not exist!\n');
    data = [];
    sweepInfo = [];
    dataAll = [];
    return
end

% Find the expected stimulation start time
stimStartExpectedOrig = mean(stimStartWindowOrig);

% Find the length of the stimulation start time window
if numel(stimStartWindowOrig) == 1
    stimStartWindowOrigLength = 0;
else
    stimStartWindowOrigLength = diff(stimStartWindowOrig);
end
    
% Baseline window in simulation time
baseWindow = timeToPad + [0, stimStartExpectedOrig];

% Get the name of the output folder
[~, outFolderName] = fileparts(outFolder);

% Create log file if requested
if createLog
    % Create log file name
    logFileName = sprintf('%s_%s.log', outFolderName, mfilename);

    % Create log file path
    logPath = fullfile(outFolder, logFileName);

    % Open the log file
    fid = fopen(logPath, 'w');
end

%% Extract information from dynamic clamp experiments if available
% Extract GABA-B IPSC waveform used in dynamic clamp
%   Note: GABA-B conductance amplitudes are converted to from nS to uS
gababAmp = swpInfoAll{fileBases, 'gabab_amp'} / NS_PER_US;
gababTrise = swpInfoAll{fileBases, 'gabab_Trise'};
gababTfallFast = swpInfoAll{fileBases, 'gabab_TfallFast'};
gababTfallSlow = swpInfoAll{fileBases, 'gabab_TfallSlow'};
gababWeight = swpInfoAll{fileBases, 'gabab_w'};

% Extract the holding voltage conditions for each file
vHold = swpInfoAll{fileBases, 'vrow'};

%% Import raw traces
% Print message
fprintf('Importing raw traces for this cell ... \n');

% Print usage message
if createLog
    for iSwp = 1:nSwps
        fprintf(fid, 'Using trace %s ... \n', fileBases{iSwp});
    end
end

% Open the matfiles
matFiles = cellfun(@matfile, filePaths, 'UniformOutput', false);

if toMedianFilter && toResample
    % Extract median-filtered & resampled data
    dataOrig = cellfun(@(x) x.d_mfrs, matFiles, 'UniformOutput', false);
else
    % Extract original data
    dataOrig = cellfun(@(x) x.d_orig, matFiles, 'UniformOutput', false);
end

% Extract data vectors
[tVecsOrig, gVecsOrig, iVecsOrig, vVecsOrig] = extract_columns(dataOrig, 1:4);

% Convert conductance vectors from nS to uS
% TODO: Make function convert_units(data, oldUnits, newUnits)
%       convert_units(gVecs, 'nS', 'uS')
gVecsOrig = cellfun(@(x) x / NS_PER_US, gVecsOrig, 'UniformOutput', false);

% Convert current vectors from pA to nA
%       TODO: convert_units(iVecs, 'pA', 'nA')
iVecsOrig = cellfun(@(x) x / PA_PER_NA, iVecsOrig, 'UniformOutput', false);

% Compute the sampling intervals in ms
siMs = compute_sampling_interval(tVecsOrig);

% Update the response window upper limit from the time vectors
if isinf(responseWindowOrig(2))
    % Get the last time point for each time vector
    timeLasts = extract_elements(tVecsOrig, 'last');

    % Set the upper limit of the response window to be 
    %   the minimum of all last time points
    responseWindowOrig(2) = min(timeLasts);
end

%% Restrict traces to the approximate response window
% Expand by stimStartWindowOrigLength to get the approximate response window
approxWindowOrig = [responseWindowOrig(1); ...
                    responseWindowOrig(2) + stimStartWindowOrigLength];

% Find the indices for the approximate current pulse response
endPointsApprox = find_window_endpoints(approxWindowOrig, tVecsOrig);

% Extract the approximate current pulse response regions
[vVecsApprox, iVecsApprox, gVecsApprox] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsApprox), ...
            vVecsOrig, iVecsOrig, gVecsOrig);

%% Process current pulses from the data vectors
if toParsePulse
    % Print message
    fprintf('Processing current pulses from the original data vectors ... \n');

    % Parse the pulse vectors
    pulseParams = parse_pulse(iVecsApprox);

    % Extract the index right after the start of the pulse for each vector
    idxStimStart = pulseParams.idxBeforeStart;

    % Extract the current pulse amplitudes
    currentPulseAmplitude = pulseParams.pulseAmplitude;
    % pulseWidth = pulseParams.pulseWidthSamples .* siMs;

    %{
    % Parse the pulse response, using the indices from the pulse
    responseParams = ...
        parse_pulse_response(vVecsApprox, siMs, ...
                                'PulseVector', ivecsApproxCpr, ...
                                'SameAsPulse', true, ...
                                'MeanValueWindowMs', meanVoltageWindow);

    % Extract the holding potentials (mV)
    holdPotential = responseParams.baseValue;

    % Extract the maximum voltage values (mV)
    maxVoltage = responseParams.maxValue;

    % Extract the voltage changes (mV)
    voltageChange = responseParams.steadyAmplitude;

    % Decide whether each trace will be used
    toUse = pulseWidth >= 0 & ...
            sign(pulseAmplitude) == sign(voltageChange) & ...
            maxVoltage <= spikeThresholdInit;

    % Count the number of traces to be used
    nToUse = sum(toUse);

    %}
else
    idxStimStart = convert_to_samples(stimStartExpectedOrig, siMs);
    currentPulseAmplitude = swpInfoAll{fileBases, 'currpulse'} / PA_PER_NA;
end

%% Fix current pulse response traces that may have out-of-balance bridges
if toParsePulse && toCorrectDcSteps
% TODO: Fix initial_slopes to output a lookup table and use it here
    % Print message
    fprintf('Fixing current pulse response traces that may have out-of-balance bridges ... \n');

    % Load find_initial_slopes.m results
    load(initialSlopesPath, 'filenamesSorted', ...
        'iThreshold1Balanced', 'iThreshold2Balanced'); 

    % Find the index of file in all files sorted by initial slope
    %   in descending order
    ftemp = @(x) find_in_strings(x, filenamesSorted);

    % Determine whether the initial slopes exceed threshold
    %   Note: These may have out-of-balance bridges
    isOutOfBalance = ...
        cellfun(@(x) ftemp(x) < iThreshold2Balanced || ...
                    ftemp(x) > iThreshold1Balanced, fileNames);

    % Print out an appropriate message
    if createLog
        if any(isOutOfBalance)
            fprintf(fid, ['The following current pulse responses will be ', ...
                            'corrected due to out-of-balance bridges:\n']);
            print_cellstr(fileNames(isOutOfBalance), ...
                          'FileID', fid, 'OmitBraces', true, ...
                          'OmitQuotes', true, 'Delimiter', '\n');
        else
            fprintf(fid, ['There are no current pulse response traces ', ...
                            'with out-of-balance bridges!\n']);
        end
    end

    % Correct for traces that may have out-of-balance bridges
    vVecsApprox = cellfun(@(x, y, z) ...
                    apply_or_return(x, @correct_unbalanced_bridge, y, z), ...
                    num2cell(isOutOfBalance), vVecsApprox, iVecsApprox, ...
                    'UniformOutput', false);
end

%% Reshape data vectors to be compared with simulations
% Print message
fprintf('Reshaping data vectors to be compared with simulations ... \n');

% Compute the expected baseline window length in ms
baseLengthMs = stimStartExpectedOrig - responseWindowOrig(1);

% Compute the expected response window length in ms
responseLengthMs = responseWindowOrig(2) - stimStartExpectedOrig;

% Convert times to samples
[baseLengthSamples, responseLengthSamples, ...
    nSamplesToPadForStabilization] = ...
    argfun(@(x) convert_to_samples(x, siMs), ...
            baseLengthMs, responseLengthMs, timeToPad);

% Compute the index to start each response vector
%   Note: this might be less than 1
idxStart = idxStimStart - baseLengthSamples + 1;

% Compute the index to end each response vector
idxEnd = idxStimStart + responseLengthSamples;

% If idxStart is less than 1, compute the number of indices to pad for baseline
nSamplesToPadBase = 1 - idxStart;
nSamplesToPadBase(nSamplesToPadBase < 0) = 0;

% Compute the total number of samples to pad before data
nSamplesToPadTotal = nSamplesToPadForStabilization + nSamplesToPadBase;

% Generate vectors to pad
vvecsToPad = arrayfun(@(x) NaN * ones(x, 1), nSamplesToPadTotal, ...
                        'UniformOutput', false);
ivecsToPad = arrayfun(@(x) zeros(x, 1), nSamplesToPadTotal, ...
                        'UniformOutput', false);
gvecsToPad = arrayfun(@(x) zeros(x, 1), nSamplesToPadTotal, ...
                        'UniformOutput', false);

% Compute the actual index to start the current pulse response for each vector
idxStartToExtract = idxStart;
idxStartToExtract(idxStart < 1) = 1;

% Put the endpoints to extract in cell array form
endPointsToExtract = ...
    force_column_vector(transpose([idxStartToExtract, idxEnd]), ...
                       'IgnoreNonVectors', false);

% Extract the response regions from the approximate response regions
[vVecsExtracted, iVecsExtracted, gVecsExtracted] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsToExtract), ...
            vVecsApprox, iVecsApprox, gVecsApprox);

% Concatenate the vectors to pad with the extracted data
vVecs = cellfun(@(x, y) vertcat(x, y), vvecsToPad, vVecsExtracted, ...
                    'UniformOutput', false);
iVecs = cellfun(@(x, y) vertcat(x, y), ivecsToPad, iVecsExtracted, ...
                    'UniformOutput', false);
gVecs = cellfun(@(x, y) vertcat(x, y), gvecsToPad, gVecsExtracted, ...
                    'UniformOutput', false);

% Compute the total number of samples for each vector
%   Note: Corresponding vectors in vVecs, iVecs, gVecs should have equal length!
nSamples = cellfun(@length, vVecs);

% Create time vectors for the final responses
tVecs = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                            'TimeUnits', 'ms', 'BoundaryMode', 'leftadjust', ...
                            'TimeStart', idxStartToExtract .* siMs, ...
                            'ForceCellOutput', true);

% Combine the vectors for output
data = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                tVecs, vVecs, iVecs, gVecs, 'UniformOutput', false);

% Save all the data
dataAll = data;

%% Average the current pulse responses according to vHold
if toParsePulse && (toAverageByVhold || toBootstrapByVhold)
    % Print message
    fprintf('Averaging the current pulse responses according to vHold ... \n');

    % Average the data by holding voltage conditions
    if toAverageByVhold
        [data, vUnique] = ...
            compute_combined_data(data, 'mean', 'Grouping', vHold, ...
                                'ColNumToCombine', 2);

        % Create a new file prefix
        filePrefix = strcat(cellName, '_vhold');

        % Define file names by the unique vHold level
        fileBases = create_labels_from_numbers(vUnique, 'Prefix', filePrefix);
        fileNames = strcat(fileBases, '.mat');
    elseif toBootstrapByVhold
        [data, vUnique] = ...
            compute_combined_data(data, 'bootmean', 'Grouping', vHold, ...
                                'ColNumToCombine', 2);

        % Create a new file prefix
        filePrefix = strcat(cellName, '_resampled_vhold_');

        % Rename files by the vHold level for each file 
        fileBases = create_labels_from_numbers(vHold, 'Prefix', filePrefix);
        fileNames = strcat(fileBases, '.mat');
    else
        error('Code logic error!');
    end

    % Update the number of sweeps
    nSwps = numel(data);

    % Extract data vectors
    [tVecs, vVecs, iVecs, gVecs] = extract_columns(data, 1:4);

    % Update the sampling interval (esp. the number of rows)
    siMs = compute_sampling_interval(tVecs);

    % Use an averaged current pulse amplitude for 
    %   the new set of current pulse responses
    currentPulseAmplitude = mean(currentPulseAmplitude) * ones(nSwps, 1);

    % Update holding potential conditions
    vHold = vUnique;

    % Set GABAB IPSC parameters to all zeros
    gababAmp = zeros(nSwps, 1);
    gababTrise = zeros(nSwps, 1);
    gababTfallFast = zeros(nSwps, 1);
    gababTfallSlow = zeros(nSwps, 1);
    gababWeight = zeros(nSwps, 1);
end

%% Compute the actual holding potentials 
% Print message
fprintf('Computing the actual holding potentials ... \n');

if toParsePulse
    % Parse the (may be averaged) current pulse responses
    responseParams = parse_pulse_response(vVecs, siMs, ...
                                        'PulseVector', iVecs, ...
                                        'SameAsPulse', true, ...
                                        'MeanValueWindowMs', meanVoltageWindow);

    % Extract the holding potentials (mV) for the current pulse responses
    holdPotential = responseParams.baseValue;
else
    % Extract actual holding potentials for each file from swpInfoAll
    %   Note: This was averaged over 20 ms before IPSC start 
    %           for the median filtered trace
    %           See m3ha_find_lts.m for specifics.
    holdPotential = swpInfoAll{fileBases, 'actVhold'};
end

%% Compute the baseline noise and sweep weights
% Print message
fprintf('Computing the baseline noise and sweep weights ... \n');

% Compute baseline noise and sweep weights
[~, ~, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tVecs, vVecs, 'BaseWindow', baseWindow);

%% Estimate the holding current and holding current noise
% Print message
fprintf('Estimate the holding current and holding current noise ... \n');

if ~isempty(epasEstimate) && ~isempty(RinEstimate)
    % Estimate the holding currents (nA) based on 
    %   estimated input resistance (MOhm) and resting membrane potential (mV)
    %   I = (V - epas) / R
    holdCurrent = (holdPotential - epasEstimate) ./ RinEstimate;

    % Estimate the corresponding variations in holding current (nA)
    holdCurrentNoise = (baseNoise ./ RinEstimate) * 5;
else
    holdCurrent = NaN(nSwps, 1);
    holdCurrentNoise = NaN(nSwps, 1);
end

%% Output results
% Print message
fprintf('Putting results into a table ... \n');

% Output in sweepInfo tables
sweepInfo = table(fileBases, fileNames, siMs, currentPulseAmplitude, ...
                    vHold, holdPotential, holdCurrent, baseNoise, ...
                    holdCurrentNoise, sweepWeights, ...
                    gababAmp, gababTrise, gababTfallFast, ...
                    gababTfallSlow, gababWeight, 'RowNames', fileBases);

% Close log file
if createLog
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cellName = m3ha_extract_cell_name (fileNames)
%% Extract the cell name from a list of file names

% Extract the common prefix
commonPrefix = extract_common_prefix(fileNames);

% 
if numel(commonPrefix) < 7
    error('Cannot find common cell name!');
end

% The cell name is the first 7 characters
cellName = commonPrefix(1:7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Unpack individual data vectors
[tVecs, vVecs, iVecs, gVecs] = extract_columns(data, 1:4);

% Find unique vHold values
vUnique = unique(vHold, 'sorted');

% Count the number of unique holding voltage conditions
nVHoldCond = length(vUnique);

% Group the voltage traces by unique vHold values, 
%   then average the grouped traces
% TODO: Use compute_average_trace with modifications to pass group
vVecsAveraged = cell(nVHoldCond, 1);
for iVhold = 1:nVHoldCond
    % Get the current vHold value
    vnow = vUnique(iVhold);

    % Collect all cpr traces with this vHold value
    vVecsGroupedThisVhold = vVecs(vHold == vnow);

    % Average the voltage traces from this vHold group
    vVecAveragedThis = compute_average_trace(vVecsGroupedThisVhold);

    % Save in arrays
    vVecsAveraged{iVhold} = vVecAveragedThis;
end
vVecs = vVecsAveraged;

% For time, current and conductance, 
%   take the first copy and repeat nVHoldCond times
%   Note: Since pulse widths are not necessarily the same, 
%           current vectors should not be averaged
[tVecs, iVecs, gVecs] = ...
    argfun(@(x) repmat(x(1), nVHoldCond, 1), tVecs, iVecs, gVecs);

% Re-combine the vectors for output
data = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                    tVecs, vVecs, iVecs, gVecs, ...
                    'UniformOutput', false);

fileNames = arrayfun(@(x) strcat(fileNames{1}(1:7), ...
                        '_vhold', num2str(x)), vUnique, ...
                        'UniformOutput', false);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
