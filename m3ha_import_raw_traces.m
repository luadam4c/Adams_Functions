function [dataCpr, dataIpscr, sweepInfoCpr, sweepInfoIpscr, dataCprAll] = ...
                m3ha_import_raw_traces (fileNames, swpInfo, initialSlopes, ...
                    cpStartWindowOrig, cprWinOrig, ...
                    ipscTimeOrig, ipscDur, ...
                    epasEstimate, RinEstimate, varargin)
%% Imports raw traces from .mat files in the m3ha format
% Usage:
% Examples:
%       TODO
% Outputs:
%       TODO
% Arguments:
%       fileNames   - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': a full directory path, 
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'OutFolder': directory to place outputs
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'CorrectDcStepsFlag': whether to correct 
%                                           unbalanced bridges in the traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AverageByVholdFlag': whether to average responses 
%                                           according to VHold
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TimeToPad': time to pad in the beginning of the trace
%                                   (same units as the sampling interval 
%                                       in the time data)
%                   must be a nonnegative scalar
%                   default == 0
%                   - 'Windows': windows to extract 
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == []
%
% Requires:
%       cd/apply_or_return.m
%       cd/argfun.m
%       cd/compute_average_trace.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sampling_interval.m
%       cd/construct_and_check_fullpath.m
%       cd/convert_to_samples.m
%       cd/correct_unbalanced_bridge.m
%       cd/create_time_vectors.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/find_ind_str_in_cell.m
%       cd/find_window_endpoints.m
%       cd/force_column_cell.m
%       cd/force_column_numeric.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/singleneuronfitting47.m and later versions
%
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
% 2018-09-12 Added outparams.AverageByVholdFlag
% 2018-11-15 Moved to Adams_Functions
% 2018-11-28 Now pads vvecsIpscr with NaN instead of with holdPotentialIpscr

%% Hard-coded constants
NS_PER_US = 1000;
PA_PER_NA = 1000;

%% Hard-coded parameters
% Parameters to be consistent with find_passive_params.m
meanVoltageWindow = 0.5;    % width in ms for calculating mean voltage 
                            %   for input resistance calculations

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
directoryDefault = '';
outFolderDefault = pwd;
correctDcStepsFlagDefault = true;
averageByVholdFlagDefault = true;
timeToPadDefault = 0;
windowsDefault = [];            % extract entire trace(s) by default

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

% Add required inputs to an input Parser
addRequired(iP, 'fileNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['fileNames must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'CorrectDcStepsFlag', correctDcStepsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AverageByVholdFlag', averageByVholdFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TimeToPad', timeToPadDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));
addParameter(iP, 'Windows', windowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Read from the input Parser
parse(iP, fileNames, varargin{:});
verbose = iP.Results.Verbose;
matFilesDir = iP.Results.Directory;
outFolder = iP.Results.OutFolder;
correctDcStepsFlag = iP.Results.CorrectDcStepsFlag;
averageByVholdFlag = iP.Results.AverageByVholdFlag;
timeToPad = iP.Results.TimeToPad;
windows = iP.Results.Windows;

%% Prepare
% Print message
fprintf('Preparing for import ... \n');

% Make sure fileNames is a column cell array
fileNames = force_column_cell(fileNames);

% Count the total number of sweeps to import
nSwpsIpscr = numel(fileNames);

% Construct full paths to matfiles and check if they exist
[filePaths, pathExists] = ...
    construct_and_check_fullpath(fileNames, 'Directory', matFilesDir, ...
                                    'Verbose', verbose);

% If a matfile does not exist, return
if ~all(pathExists)
    fprintf('One or more matfiles does not exist!\n');
    dataCpr = [];
    dataIpscr = [];
    sweepInfoCpr = [];
    sweepInfoIpscr = [];
    dataCprAll = [];
    return
end

% Find the expected current pulse start time
cpStartExpectedOrig = mean(cpStartWindowOrig);

% Compute the expected current pulse baseline window length in ms
cprBaseLengthMs = cpStartExpectedOrig - cprWinOrig(1);

% Compute the expected current pulse response window length in ms
cprResponseLengthMs = cprWinOrig(2) - cpStartExpectedOrig;

% Find the length of the current pulse start time window
cpStartWindowOrigLength = cpStartWindowOrig(2) - cpStartWindowOrig(1);

% Expand by cpStartWindowOrigLength to get 
%   the approximate current pulse response window
approxCprWindow = [cprWinOrig(1); cprWinOrig(2) + cpStartWindowOrigLength];

% Baseline window for the IPSC response in original time
baseWindowIpscrOrig = [0, ipscTimeOrig];

% Baseline window for the current pulse response in simulation time
baseWindowCpr = timeToPad + [0, cpStartExpectedOrig];

% Get the name of the output folder
[~, outFolderName] = fileparts(outFolder);

% Create log file name
logFileName = sprintf('%s_%s.log', outFolderName, mfilename);

% Create log file path
logPath = fullfile(outFolder, logFileName);

% Open the log file
fid = fopen(logPath, 'w');

%% Import raw traces
% Print message
fprintf('Importing raw traces for this cell ... \n');

% Print usage message
for iSwp = 1:nSwpsIpscr
    fprintf(fid, 'Using trace %s ... \n', fileNames{iSwp});
end

% Open the matfiles
matFiles = cellfun(@matfile, filePaths, 'UniformOutput', false);

% Extract original data
dataOrig = cellfun(@(x) x.d_orig, matFiles, 'UniformOutput', false);

% Extract median-filtered & resampled data
dataMfrs = cellfun(@(x) x.d_mfrs, matFiles, 'UniformOutput', false);

% Extract original data vectors
[tvecsOrig, gvecsOrig, ivecsOrig, vvecsOrig] = extract_columns(dataOrig, 1:4);

% Extract median-filtered & resampled data vectors
[tvecsMfrs, gvecsMfrs, ivecsMfrs, vvecsMfrs] = extract_columns(dataMfrs, 1:4);

% Convert conductance vectors from nS to uS
% TODO: Make function convert_units(data, oldUnits, newUnits)
%       convert_units(gvecsOrig, 'nS', 'uS')
gvecsOrig = cellfun(@(x) x / NS_PER_US, gvecsOrig, 'UniformOutput', false);
gvecsMfrs = cellfun(@(x) x / NS_PER_US, gvecsMfrs, 'UniformOutput', false);

% Convert current vectors from pA to nA
ivecsOrig = cellfun(@(x) x / PA_PER_NA, ivecsOrig, 'UniformOutput', false);
ivecsMfrs = cellfun(@(x) x / PA_PER_NA, ivecsMfrs, 'UniformOutput', false);

% Compute the sampling intervals in ms
siMsOrig = compute_sampling_interval(tvecsOrig);
siMsMfrs = compute_sampling_interval(tvecsMfrs);

%% Process current pulses from the original data vectors
% Print message
fprintf(['Processing current pulses from the original data vectors ... \n']);

% Find the indices for the approximate current pulse response
endPointsApproxCpr = find_window_endpoints(approxCprWindow, tvecsOrig);

% Extract the approximate current pulse response regions
[vvecsApproxCpr, ivecsApproxCpr, gvecsApproxCpr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsApproxCpr), ...
            vvecsOrig, ivecsOrig, gvecsOrig);

% Parse the pulse vectors
[pulseParams, ~] = parse_pulse(ivecsApproxCpr);

% Extract the index right after the start of the pulse for each vector
idxCpStart = pulseParams.idxAfterStart;

% Extract the current pulse widths and amplitudes
pulseWidth = pulseParams.pulseWidthSamples .* siMsOrig;
pulseAmplitude = pulseParams.pulseAmplitude;

% Store these current pulse widths and amplitudes for the IPSC response
currentPulseAmplitudeIpscr = pulseAmplitude;

%{
% Parse the pulse response, using the indices from the pulse
[responseParams, responseData] = ...
    parse_pulse_response(vvecsApproxCpr, siMsOrig, ...
                            'PulseVector', ivecsApproxCpr, ...
                            'SameAsPulse', true, ...
                            'MeanValueWindowMs', meanVoltageWindow);

% Extract the holding potentials (mV)
holdPotentialCpr = responseParams.baseValue;

% Extract the maximum voltage values (mV)
maxVoltageCpr = responseParams.maxValue;

% Extract the voltage changes (mV)
voltageChangeCpr = responseParams.steadyAmplitude;

% Decide whether each trace will be used
toUse = pulseWidth >= 0 & ...
        sign(pulseAmplitude) == sign(voltageChange) & ...
        maxVoltage <= spikeThresholdInit;

% Count the number of traces to be used
nToUse = sum(toUse);

%}

%% Fix current pulse response traces that may have out-of-balance bridges
if correctDcStepsFlag
% TODO: Fix initial_slopes to output a lookup table and use it here
    % Print message
    fprintf('Fixing current pulse response traces that may have out-of-balance bridges ... \n');

    % Read from initialSlopes
    initialSlopeFilenames = initialSlopes.filenamesSorted;
    initialSlopeThreshold1IndexBalanced = initialSlopes.iThreshold1Balanced;
    initialSlopeThreshold2IndexBalanced = initialSlopes.iThreshold2Balanced;

    % Find the index of file in all files sorted by initial slope
    %   in descending order
    ftemp = @(x) find_ind_str_in_cell(x, initialSlopeFilenames);

    % Determine whether the initial slopes exceed threshold
    %   Note: These may have out-of-balance bridges
    isOutOfBalance = ...
        cellfun(@(x) ftemp(x) < initialSlopeThreshold2IndexBalanced || ...
                    ftemp(x) > initialSlopeThreshold1IndexBalanced, fileNames);

    % Print out an appropriate message
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

    % Correct for traces that may have out-of-balance bridges
    vvecsApproxCpr = ...
        cellfun(@(x, y, z) ...
                apply_or_return(x, @correct_unbalanced_bridge, y, z), ...
                num2cell(isOutOfBalance), vvecsApproxCpr, ivecsApproxCpr, ...
                'UniformOutput', false);
end

%% Construct current pulse response vectors to be compared with simulations
% Print message
fprintf('Constructing current pulse response vectors to be compared with simulations ... \n');

% Convert times to samples
[cprBaseLengthSamples, cprResponseLengthSamples, ...
    nSamplesToPadForStabilization] = ...
    argfun(@(x) convert_to_samples(x, siMsOrig), ...
            cprBaseLengthMs, cprResponseLengthMs, timeToPad);

% Compute the index to start the current pulse response for each vector
%   Note: this might be less than 1
idxCprStart = idxCpStart - cprBaseLengthSamples;

% Compute the index to end the current pulse response for each vector
idxCprEnd = idxCpStart + cprResponseLengthSamples - 1;

% If idxCprStart is less than 1, 
%   compute the number of indices to pad before current pulse response data
nSamplesToPadCprBase = 1 - idxCprStart;
nSamplesToPadCprBase(nSamplesToPadCprBase < 0) = 0;

% Compute the total number of samples to pad before data
nSamplesToPadCpr = nSamplesToPadForStabilization + nSamplesToPadCprBase;

% Generate vectors to pad
vvecsToPadCpr = arrayfun(@(x) NaN * ones(x, 1), nSamplesToPadCpr, ...
                        'UniformOutput', false);
ivecsToPadCpr = arrayfun(@(x) zeros(x, 1), nSamplesToPadCpr, ...
                        'UniformOutput', false);
gvecsToPadCpr = arrayfun(@(x) zeros(x, 1), nSamplesToPadCpr, ...
                        'UniformOutput', false);

% Compute the actual index to start the current pulse response for each vector
idxCprStartToExtract = idxCprStart;
idxCprStartToExtract(idxCprStart < 1) = 1;

% Put the endpoints to extract in cell array form
endPointsToExtractCpr = ...
    force_column_numeric(transpose([idxCprStartToExtract, idxCprEnd]), ...
                       'IgnoreNonVectors', false);

% Extract the current pulse response regions from 
%   approximate current pulse response regions
[vvecsExtractedCpr, ivecsExtractedCpr, gvecsExtractedCpr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsToExtractCpr), ...
            vvecsApproxCpr, ivecsApproxCpr, gvecsApproxCpr);

% Concatenate the vectors to pad with the extracted data
vvecsCpr = cellfun(@(x, y) vertcat(x, y), vvecsToPadCpr, vvecsExtractedCpr, ...
                    'UniformOutput', false);
ivecsCpr = cellfun(@(x, y) vertcat(x, y), ivecsToPadCpr, ivecsExtractedCpr, ...
                    'UniformOutput', false);
gvecsCpr = cellfun(@(x, y) vertcat(x, y), gvecsToPadCpr, gvecsExtractedCpr, ...
                    'UniformOutput', false);

% Compute the total number of samples for each vector
%   Note: Corresponding vectors in vvecsCpr, ivecsCpr, gvecsCpr 
%           should have the same length!
nSamplesCpr = cellfun(@length, vvecsCpr);

% Create time vectors for the final current pulse responses
tvecsCpr = create_time_vectors(nSamplesCpr, 'TimeUnits', 'ms', ...
                                'SamplingIntervalMs', siMsOrig);

% Combine the vectors for output
dataCpr = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                    tvecsCpr, vvecsCpr, ivecsCpr, gvecsCpr, ...
                    'UniformOutput', false);

% Save all the data
dataCprAll = dataCpr;

%% Average the current pulse responses according to vHold
% TODO: Pull out to its own function and use in m3ha_run_neuron_once.m
if AverageByVholdFlag
    % Print message
    fprintf('Averaging the current pulse responses according to vHold ... \n');

    % Extract holding voltage conditions for each file from swpInfo
    vHoldCond = swpInfo{fileNames, 'vrow'};

    % Unpack individual data vectors
    [tvecsCpr, vvecsCpr, ivecsCpr, gvecsCpr] = extract_columns(dataCpr, 1:4);

    % Find unique vHold values
    vUnique = unique(vHoldCond, 'sorted');

    % Count the number of unique holding voltage conditions
    nVHoldCond = length(vUnique);

    % Group the voltage traces by unique vHold values, 
    %   then average the grouped traces
    % TODO: Use compute_average_trace with modifications to pass group
    vvecsCprAveraged = cell(nVHoldCond, 1);
    for iVhold = 1:nVHoldCond
        % Get the current vHoldCond value
        vnow = vUnique(iVhold);

        % Collect all cpr traces with this vHoldCond value
        vvecsGroupedThisVhold = vvecsCpr(vHoldCond == vnow);

        % Average the voltage traces from this vHoldCond group
        vvecCprAveragedThis = compute_average_trace(vvecsGroupedThisVhold);

        % Save in arrays
        vvecsCprAveraged{iVhold} = vvecCprAveragedThis;
    end
    vvecsCpr = vvecsCprAveraged;

    % For time, current and conductance, 
    %   take the first copy and repeat nVHoldCond times
    %   Note: Since pulse widths are not necessarily the same, 
    %           current vectors should not be averaged
    tvecsCpr = repmat({tvecsCpr{1}}, nVHoldCond, 1);
    ivecsCpr = repmat({ivecsCpr{1}}, nVHoldCond, 1);
    gvecsCpr = repmat({gvecsCpr{1}}, nVHoldCond, 1);

    % Re-combine the vectors for output
    dataCpr = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                        tvecsCpr, vvecsCpr, ivecsCpr, gvecsCpr, ...
                        'UniformOutput', false);

    % Define file names by the vhold level
    fileNamesCpr = arrayfun(@(x) strcat(fileNames{1}(1:7), ...
                            '_vhold', num2str(x)), vUnique, ...
                            'UniformOutput', false);

    % Count the number of sweeps for the current pulse response
    nSwpsCpr = numel(dataCpr);

    % Use an averaged current pulse amplitude for 
    %   the averaged current pulse response
    currentPulseAmplitudeCpr = mean(pulseAmplitude) * ones(nSwpsCpr, 1);
else
    dataCpr = dataCpr;
    fileNamesCpr = fileNames;
    currentPulseAmplitudeCpr = pulseAmplitude;
end

%% Construct IPSC response vectors to be compared with simulations
% Print message
fprintf('Constructing IPSC response vectors to be compared with simulations ... \n');

% Convert times to samples
[ipscTimeOrigSamples, ipscDurSamples, nSamplesToPadIpscr] = ...
    argfun(@(x) convert_to_samples(x, siMsMfrs), ...
            ipscTimeOrig, ipscDur, timeToPad);

% Generate vectors to pad
vvecsToPadIpscr = arrayfun(@(x) NaN * ones(x, 1), nSamplesToPadIpscr, ...
                            'UniformOutput', false);
ivecsToPadIpscr = arrayfun(@(x) zeros(x, 1), nSamplesToPadIpscr, ...
                            'UniformOutput', false);
gvecsToPadIpscr = arrayfun(@(x) zeros(x, 1), nSamplesToPadIpscr, ...
                            'UniformOutput', false);
                        
% Decide on the endpoints to extract for the IPSC response
endPointsToExtractIpscr = ...
    arrayfun(@(x, y) [1; x + y], ipscTimeOrigSamples, ipscDurSamples, ...
                    'UniformOutput', false);

% Extract the IPSC response regions from the original data vectors
[vvecsExtractedIpscr, ivecsExtractedIpscr, gvecsExtractedIpscr] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPointsToExtractIpscr), ...
            vvecsMfrs, ivecsMfrs, gvecsMfrs);

% Concatenate the vectors to pad with the extracted data
vvecsIpscr = cellfun(@(x, y) vertcat(x, y), vvecsToPadIpscr, vvecsExtractedIpscr, ...
                    'UniformOutput', false);
ivecsIpscr = cellfun(@(x, y) vertcat(x, y), ivecsToPadIpscr, ivecsExtractedIpscr, ...
                    'UniformOutput', false);
gvecsIpscr = cellfun(@(x, y) vertcat(x, y), gvecsToPadIpscr, gvecsExtractedIpscr, ...
                    'UniformOutput', false);

% Compute the total number of samples for each vector
%   Note: Corresponding vectors in vvecsIpscr, ivecsIpscr, gvecsIpscr 
%           should have the same length!
nSamplesIpscr = cellfun(@length, vvecsIpscr);

% Create time vectors for the final current pulse responses
tvecsIpscr = create_time_vectors(nSamplesIpscr, 'TimeUnits', 'ms', ...
                                'SamplingIntervalMs', siMsMfrs);

% Combine the vectors for output
dataIpscr = cellfun(@(x, y, z, w) horzcat(x, y, z, w), ...
                    tvecsIpscr, vvecsIpscr, ivecsIpscr, gvecsIpscr, ...
                    'UniformOutput', false);

% Store the file names
fileNamesIpscr = fileNames;

%% Compute the actual holding potentials 
% Print message
fprintf('Computing the actual holding potentials ... \n');

% Temporary measure before parse_pulse_response is fixed
% TODO:
siMsOrigTemp = mean(siMsOrig);

% Parse the (may be averaged) current pulse responses
responseParams = parse_pulse_response(vvecsCpr, siMsOrigTemp, ...
                                    'PulseVector', ivecsCpr, ...
                                    'SameAsPulse', true, ...
                                    'MeanValueWindowMs', meanVoltageWindow);

% Extract the holding potentials (mV) for the current pulse responses
holdPotentialCpr = responseParams.baseValue;

% Extract actual holding potentials for each file from swpInfo
%   Note: This was averaged over 20 ms before IPSC start 
%           for the median filtered trace
%           See find_LTS.m for specifics.
holdPotentialIpscr = swpInfo{fileNames, 'actVhold'};

%% Compute the baseline noise and sweep weights
% Print message
fprintf('Computing the baseline noise and sweep weights ... \n');

% Compute baseline noise and sweep weights for the current pulse response
%   Note: Can't use original data because they might be averaged
[~, ~, baseNoiseCpr, sweepWeightsCpr] = ...
    compute_default_sweep_info(tvecsCpr, vvecsCpr, 'BaseWindow', baseWindowCpr);

% Compute baseline noise and sweep weights for the IPSC response
[~, ~, baseNoiseIpscr, sweepWeightsIpscr] = ...
    compute_default_sweep_info(tvecsMfrs, vvecsMfrs, ...
            'BaseWindow', baseWindowIpscrOrig);

%% Estimate the holding current and holding current noise
% Print message
fprintf('Estimate the holding current and holding current noise ... \n');

% Estimate the holding currents (nA) based on 
%   estimated input resistance (MOhm) and resting membrane potential (mV)
%   I = (V - epas) / R
holdCurrentIpscr = (holdPotentialIpscr - epasEstimate) / RinEstimate;
holdCurrentCpr = (holdPotentialCpr - epasEstimate) / RinEstimate;

% Estimate the corresponding variations in holding current (nA)
holdCurrentNoiseIpscr = (baseNoiseIpscr / RinEstimate) * 5;
holdCurrentNoiseCpr = (baseNoiseCpr / RinEstimate) * 5;

%% Save results
% Print message
fprintf('Saving results ... \n');

% Save in sweepInfo tables
sweepInfoCpr = table(fileNamesCpr, currentPulseAmplitudeCpr, ...
                    holdPotentialCpr, holdCurrentCpr, baseNoiseCpr, ...
                    holdCurrentNoiseCpr, sweepWeightsCpr);
sweepInfoIpscr = table(fileNamesIpscr, currentPulseAmplitudeIpscr, ...
                    holdPotentialIpscr, holdCurrentIpscr, baseNoiseIpscr, ...
                    holdCurrentNoiseIpscr, sweepWeightsIpscr);

% Close log file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find approximate current pulse response window
% Must be sure to include current pulse start and at least 500 ms after current pulse start
acprwinBegin = find(m.d_orig(:, 1) >= cpStartWindowOrig(1), 1);
cpStartWindowOrigLength = cpStartWindowOrig(2) - cpStartWindowOrig(1);
acprwinEnd = find(m.d_orig(:, 1) <= cprWinOrig(2) + cpStartWindowOrigLength/2, 1, 'last');
acprwinInd = acprwinBegin:acprwinEnd;       % indices corresponding to approximate cpr window

% Find indices for current pulse response relative to tvecCpr
idxCprStart = idxCpStart;
cprWinLengthSamples = cprWinOrig(2) - cprWinOrig(1);
idxCprEnd = idxCprStart + round(cprWinLengthSamples/siMsOrig);
cprInd = (idxCprStart + 1):idxCprEnd;       % since we are padding later, we need the +1

% Find holding potential right before current pulse application
mvind = 1:round(meanVoltageWindow/siMsOrig);               % indices for taking the mean of voltages
baseInd = idxCpStart - round(0.5/siMsOrig) - fliplr(mvind);      % base indices
cprBaselineVoltage = mean(vvecCpr(baseInd));
outparams.holdPotentialCpr(ct, 1) = cprBaselineVoltage;

% Pad raw traces so that current pulses lie at 2100-2110 ms and response goes on until 2500 ms
about2100 = round(outparams.cprWindow(1)/siMsOrig) * siMsOrig;       % about 2100 ms
                % this is needed because siMsOrig is often 0.099 ms
tvecToPad = (siMsOrig:siMsOrig:about2100)';
vvecToPad = outparams.holdPotentialCpr(ct, 1) * ones(1, length(tvecToPad))';
ivecToPad = zeros(1, length(tvecToPad))';
gvecToPad = zeros(1, length(tvecToPad))';
tvecShifted = tvecCpr(cprInd) - tvecCpr(idxCprStart - 1) + about2100;
vvecShifted = vvecCpr(cprInd);
ivecShifted = ivecCpr(cprInd);
gvecShifted = gvecCpr(cprInd);
dataCpr{ct}(:, 1) = vertcat(tvecToPad, tvecShifted);    % time vector of current pulse response
dataCpr{ct}(:, 2) = vertcat(vvecToPad, vvecShifted);    % voltage vector of current pulse response
dataCpr{ct}(:, 3) = vertcat(ivecToPad, ivecShifted);    % current vector of current pulse response
dataCpr{ct}(:, 4) = vertcat(gvecToPad, gvecShifted);    % conductance vector of current pulse response

holdPotentialCpr = zeros(nSweeps, 1);          % stores the holding potentials right before cp application (mV)

%       /home/Matlab/Kojis_Functions/rms_Gaussian.m
    addpath(fullfile(functionsDirectory, 'Kojis_Functions')); 
                                            % for rms_Gaussian.m
    baseNoiseIpscr(ct) = rms_Gaussian(baseline);
    cprBaselineNoise = rms_Gaussian(baselineCprThis);

% Take out traces that may have out-of-balance bridges if requested
if outparams.correctDcStepsFlag
    % Determine whether the initial slopes exceed threshold
    %   Note: These may have out-of-balance bridges
    isOutOfBalance = ...
        cellfun(@(x) find_ind_str_in_cell(x, initialSlopeFilenames) < ...
                                        initialSlopeThresholdIndexToInclude, ...
                    fileNamesThisVhold);

    % Print out appropriate message
    if any(isOutOfBalance)
        fprintf(fid, ['For Vhold == %g, ', ...
                 'the following current pulse responses will be taken out ', ...
                 'due to out-of-balance bridges:\n'], vnow);
        print_cellstr(fileNamesThisVhold(isOutOfBalance), ...
                      'FileID', fid, 'OmitBraces', true, ...
                      'OmitQuotes', true, 'Delimiter', '\n');
    else
        fprintf(fid, ['For Vhold == %g, ', ...
                 'there are no out-of-balance bridges!\n'], vnow);
    end

    % Take out traces that may have out-of-balance bridges
    dCprGroupedThisVhold = dCprGroupedThisVhold(~isOutOfBalance);
end

initialSlopeThresholdIndexToInclude = initialSlopes.iThreshold1Balanced;

holdCurrentIpscr = zeros(nSweeps, 1);                % stores the holding currents recorded in Christine's Excel sheet (nA)

% Average the holding currents
holdCurrentThisVhold = mean(holdCurrentIpscr(vHoldCond == vnow));

baseNoiseCpr = zeros(1, nVHoldCond);

baseNoiseIpscr = zeros(1, nSweeps);                  % stores baseline noise in voltage recordings

isOutOfBalance = ...
    cellfun(@(x) ftemp(x) < initialSlopeThreshold2IndexBalanced || ...
                ftemp(x) > initialSlopeThreshold1IndexBalanced, ...
                fileNamesThisVhold);

% Correct for traces that may have out-of-balance bridges
parfor iTrace = 1:nTracesThisVhold
    if isOutOfBalance(iTrace)
        % Get the old trace
        vvecOld = dCprGroupedThisVhold{iTrace}(:, 2);
        ivecOld = dCprGroupedThisVhold{iTrace}(:, 3);

        % Correct any unbalanced bridge in the trace
        vvecNew = correct_unbalanced_bridge(vvecOld, ivecOld);

        % Store the new trace
        dCprGroupedThisVhold{iTrace}(:, 2) = vvecNew;
    end
end

print_cellstr(fileNamesThisVhold(isOutOfBalance), ...
              'FileID', fid, 'OmitBraces', true, ...
              'OmitQuotes', true, 'Delimiter', '\n');

% Get the corresponding fileNames
fileNamesThisVhold = fileNames(vHoldCond == vnow);


% Number of current pulse response traces to fit
nSwpsCpr = nVHoldCond;

% Find the number of traces
nTracesThisVhold = numel(dCprGroupedThisVhold);

%% Get other information for each sweep
% Find the holding current used by Christine for each sweep, 
%   converting from pA to nA
holdCurrentIpscr = actIhold(swpIdxIpscr)' / 1000;

% Compute the regression coefficient
% RinRegression = (holdCurrentIpscr - mean(holdCurrentIpscr)) \ ...
%                     (holdPotentialIpscr - mean(holdPotentialIpscr));

%% Estimate epas and Rin with linear least squares
% Define the matrices
% Ohm's Law: IR + epas = V
%     units: [nA] * [MOhm] + [mV] = [mV]
% X * w = V
X = ones(nSweeps, 2);
X(:, 1) = holdCurrentIpscr;
V = holdPotentialIpscr;

% Compute the estimates
w = pinv(X) * V;
RinEstimate = w(1);                 % input resistance (MOhm)
epasEstimate = w(2);                % resting membrane potential (mV)

% Check if the values make sense
if RinEstimate <= 0
    colorRin = 'r';
else
    colorRin = 'iCol';
end

% Construct a vector of holding currents
holdCurrentToPlot = linspace(min(holdCurrentIpscr), max(holdCurrentIpscr), 1000);

% Compute predicted values
holdPotentialPredicted = holdCurrentToPlot * RinEstimate + epasEstimate;

% Plot holding potential versus holding current
h = figure('Visible', 'off');
clf(h);
hold on;
plot(holdCurrentIpscr, holdPotentialIpscr, 'o', 'LineWidth', 2);
plot(holdCurrentToPlot, holdPotentialPredicted, 'r')
text(0.1, 0.9, ['Rin = ', num2str(RinEstimate), ' MOhm'], ...
    'Units', 'normalized', 'Color', colorRin);
text(0.1, 0.85, ['epas = ', num2str(epasEstimate), ' mV'], ...
    'Units', 'normalized', 'Color', 'iCol');
% text(0.1, 0.8, ['Slope = ', num2str(RinRegression), ' MOhm'], ...
%     'Units', 'normalized', 'Color', 'r');
title(['Voltage-Current relationship for ', outparams.cellName]);
ylabel('Holding potential (mV)');
xlabel('Holding current (nA)');
figName = fullfile(outparams.outFolderName, ...
                    [outparams.prefix, ...
                    '_voltage-current-relationship.png']);
saveas(h, figName);
close(h)

global outparams

outparams.fileNames = fileNames;
outparams.swpIdxIpscr = swpIdxIpscr;
outparams.currentPulseAmplitude = currentPulseAmplitude;

outparams.holdPotentialCpr = holdPotentialCpr;
outparams.holdCurrentCpr = holdCurrentCpr;
outparams.baseNoiseCpr = baseNoiseCpr;
outparams.holdCurrentNoiseCpr = holdCurrentNoiseCpr;
outparams.sweepWeightsCpr = sweepWeightsCpr;

outparams.holdPotentialIpscr = holdPotentialIpscr;
outparams.holdCurrentIpscr = holdCurrentIpscr;
outparams.baseNoiseIpscr = baseNoiseIpscr;
outparams.holdCurrentNoiseIpscr = holdCurrentNoiseIpscr;
outparams.sweepWeightsIpscr = sweepWeightsIpscr;

if ~isempty(swpIndG200P{1})
    swpIdx = swpIndG200P{iRow}(iCol);
elseif ~isempty(swpIndPCond{1})
    swpIdx = swpIndPCond{iRow}(iCol);
elseif ~isempty(swpIndRow{1})
    swpIdx = swpIndRow{iRow}(iCol);
end

% Get the current file name
fileName = fnrow{swpIdx};

fileNames = cell(nSweeps, 1);                   % stores file names of traces used

% Store in arrays
swpIdxIpscr(ct) = swpIdx;
fileNames{ct} = fileName;

swpIdxIpscr = zeros(nSweeps, 1);                % stores the sweep index for each trace used

for iRow = 1:nRows
    for iCol = 1:nColumns
        % Increment the sweep count
        ct = ct + 1;

        % Get the current file name
        fileName = fileNamesRowCol{iRow, iCol};

    end
end

if ct ~= nSweeps
    error('Number of traces imported incorrect!!');
else
    fprintf('\n');
end


% Get the swpIdx from the file name
swpIdx = find_ind_str_in_cell(fileName, fnrow);

fnrow = swpInfo.fnrow;

nSweeps = nRows * nColumns;

filePath = filePaths{iSwp};
% Open the matfile
m = matfile(filePath);

% Use original data
tvecOrig = m.d_orig(:, 1);

indCpr = transpose(endPoints(1):endPoints(2));

currpulse = swpInfo.currpulse;

currentPulseAmplitudeIpscr = currentPulseAmplitude;
currentPulseAmplitudeCpr = currentPulseAmplitude;

% Compute desired window in which the current pulse response would lie (ms)
cprWindow = cprWinOrig + timeToStabilizeMs;

% Compute desired time of IPSC application (ms)
timeIPSC = ipscTimeOrig + timeToStabilizeMs;

nSamplesToPadForStabilization = round(cprWindow(1) ./ siMsOrig);

vrow = swpInfo.vrow;
actVhold = swpInfo.actVhold;

% Compute the expected current pulse baseline window length in samples
cprBaseLengthSamples = round(cprBaseLengthMs ./ siMsOrig);

% Compute the expected current pulse response window length in samples
cprResponseLengthSamples = round(cprResponseLengthMs ./ siMsOrig);

% Compute the number of samples to pad
nSamplesToPadForStabilization = round(timeToStabilizeMs ./ siMsOrig);

% Initialize output variables
dataCpr = cell(nSweeps, 1);                     % stores current pulse response traces for fitting
dataIpscr = cell(nSweeps, 1);                   % stores IPSC response traces for fitting
currentPulseAmplitude = zeros(nSweeps, 1);      % stores the current pulse amplitude (nA)
% cpstart = zeros(nSweeps, 1);                  % stores the time of current pulse application (ms)
baseNoiseIpscr = zeros(nSweeps, 1);             % stores baseline noise in voltage recordings
holdCurrentIpscr = zeros(nSweeps, 1);           % stores estimated holding currents
holdPotentialIpscr = zeros(nSweeps, 1);         % stores the holding potentials right before IPSC application (mV)
ct = 0;                                         % counts number of raw traces imported

for iSwp = 1:nSweeps
    % Print message
    fprintf('Using trace %s ... \n', fileNames{iSwp});

    m = matFiles{iSwp};

    % Find current pulse response window
    %   Must be sure to include current pulse start and cpStartWindowOrigLength adjustment
    acprwinBegin = find(tvecOrig >= cprWinOrig(1), 1);
    acprwinEnd = find(tvecOrig <= cprWinOrig(2) + cpStartWindowOrigLength, 1, 'last');
    acprwinInd = acprwinBegin:acprwinEnd;       % indices corresponding to approximate cpr window

    % Use original data for current pulse response
    siMsOrig = tvecOrig(2) - tvecOrig(1);        % sampling interval in ms
    tvecCpr = tvecOrig(acprwinInd);             % time vector of original data in ms
    gvecCpr = m.d_orig(acprwinInd, 2);          % conductance vector of original data in nS
    gvecCpr = gvecCpr / 1000;                   % conductance vector in uS
    ivecCpr = m.d_orig(acprwinInd, 3);          % current vector of original data in pA
    ivecCpr = ivecCpr / 1000;                   % current vector in nA
    vvecCpr = m.d_orig(acprwinInd, 4);          % voltage vector of original data in mV

    % Find current pulse amplitude (convert to nA) and start of current pulse application
    currentPulseAmplitude(ct, 1) = currpulse(swpIdx) / PA_PER_NA;
    cpStartWindowOrigEnd = find(tvecCpr <= cpStartWindowOrig(2), 1, 'last');
    idxCpStart = find(ivecCpr(1:cpStartWindowOrigEnd) > currentPulseAmplitude(ct, 1) * 0.25, 1, 'last');
%         cpstart(ct, 1) = tvecCpr(idxCpStart);                 % current pulse start in ms

    % Find the expected baseline window length in samples
    cprBaseLengthSamples = round((cpStartExpectedOrig - cprWinOrig(1))/siMsOrig);

    % Find the number of indices to pad before current pulse response data
    idxCprStart = idxCpStart - cprBaseLengthSamples;
    if idxCprStart < 1
        nSamplesToPadCprBase = 1 - idxCprStart;
    else
        nSamplesToPadCprBase = 0;
    end

    % Find indices for current pulse response relative to tvecCpr
    cprWinLengthSamples = round((cprWinOrig(2) - cprWinOrig(1))/siMsOrig);
    idxCprEnd = idxCprStart + cprWinLengthSamples - 1;
    if idxCprStart < 1
        cprInd = 1:idxCprEnd;
    else
        cprInd = idxCprStart:idxCprEnd;
    end

    % Get full time vector of current pulse response
    tvecCprLength = round(cprWindow(2)/siMsOrig);
    tvecCprFull = siMsOrig * (1:tvecCprLength)';
    dataCpr{ct}(:, 1) = tvecCprFull;                           % time vector of current pulse response

    % Pad raw traces so that current pulses lie at 2100-2110 ms
    nSamplesToPadForStabilization = round(cprWindow(1)/siMsOrig);
    nSamplesToPadCpr = nSamplesToPadForStabilization + nSamplesToPadCprBase;
    vvecToPadCpr = NaN * ones(1, nSamplesToPadCpr)';
    ivecToPadCpr = zeros(1, nSamplesToPadCpr)';
    gvecToPadCpr = zeros(1, nSamplesToPadCpr)';
    dataCpr{ct}(:, 2) = vertcat(vvecToPadCpr, vvecCpr(cprInd));   % voltage vector of current pulse response
    dataCpr{ct}(:, 3) = vertcat(ivecToPadCpr, ivecCpr(cprInd));   % current vector of current pulse response
    dataCpr{ct}(:, 4) = vertcat(gvecToPadCpr, gvecCpr(cprInd));   % conductance vector of current pulse response

    % Use median-filtered & resampled data for IPSC response
    tvec = m.d_mfrs(:, 1);  % time vector of median-filtered then resampled data in ms
    gvec = m.d_mfrs(:, 2);  % conductance vector of median-filtered then resampled data in nS
    gvec = gvec / 1000;     % conductance vector in uS
    ivec = m.d_mfrs(:, 3);  % current vector of median-filtered then resampled data in pA
    ivec = ivec / 1000;     % current vector in nA
    vvec = m.d_mfrs(:, 4);  % voltage vector of median-filtered then resampled data in mV

    % Holding potential was already extracted during data analysis
    holdPotentialIpscr(ct, 1) = actVhold(swpIdx);    

    % Estimate the holding currents (nA) based on estimated
    %   input resistance (MOhm) and resting membrane potential (mV)
    %   I = (V - epas) / R
    holdCurrentIpscr = (holdPotentialIpscr - epasEstimate) / RinEstimate;

    % Pad raw traces so that IPSCs begin at 3000 ms (timeIPSC) 
    %            and goes on for 7000 ms (ipscDur)
    sims = tvec(2) - tvec(1);       % Should be 1 ms
    about3000 = round(timeIPSC/sims)*sims;
    tvecToPadIpscr = (sims:sims:about3000)';
    vvecToPadIpscr = holdPotentialIpscr(ct, 1) * ones(1, length(tvecToPadIpscr))';
    ivecToPadIpscr = zeros(1, length(tvecToPadIpscr))';
    gvecToPadIpscr = zeros(1, length(tvecToPadIpscr))';
    indofipsc = round(ipscTimeOrig/sims);
    indofend = round((ipscTimeOrig + ipscDur)/sims);
    tvecShifted = tvec(indofipsc:indofend) - tvec(indofipsc - 1) + about3000;
    vvecShifted = vvec(indofipsc:indofend);
    ivecShifted = ivec(indofipsc:indofend);
    gvecShifted = gvec(indofipsc:indofend);
    dataIpscr{ct}(:, 1) = vertcat(tvecToPadIpscr, tvecShifted);        % time vector of IPSC response
    dataIpscr{ct}(:, 2) = vertcat(vvecToPadIpscr, vvecShifted);        % voltage vector of IPSC response
    dataIpscr{ct}(:, 3) = vertcat(ivecToPadIpscr, ivecShifted);        % current vector of IPSC response
    dataIpscr{ct}(:, 4) = vertcat(gvecToPadIpscr, gvecShifted);        % conductance vector of IPSC response

    % Find the baseline noise
    indBaseline = 1:indofipsc-1;
    baseline = vvec(indBaseline);
    baseNoiseIpscr(ct) = compute_rms_error(baseline);
end

% Find the holding voltage assigned for each sweep
vHoldCond = vrow(swpIdxIpscr);


endPointsToExtractIpscr = [1, ipscTimeOrigSamples + ipscDurSamples];

endPointsToExtractIpscr = [ones(nSweeps, 1), ...
                            ipscTimeOrigSamples + ipscDurSamples];
endPointsToExtractIpscr = ...
    force_column_numeric(transpose(endPointsToExtractIpscr), ...
                       'IgnoreNonVectors', false);

% vvecsToPadIpscr = arrayfun(@(x) x * ones(nSamplesToPadIpscr, 1), ...
%                             holdPotentialIpscr, 'UniformOutput', false);

holdCurrentCpr = zeros(nSwpsCpr, 1);
% Estimate the holding currents (nA) based on estimated
%   input resistance and resting membrane potential
%   I = (V - epas) / R
cprHoldCurrent = (cprBaselineVoltage - epasEstimate) / RinEstimate;
holdCurrentCpr(iSwp) = cprHoldCurrent;

% Compute sweep weights based on baseline noise
sweepWeightsCpr = 1 ./ baseNoiseCpr;
sweepWeightsIpscr = 1 ./ baseNoiseIpscr;

%% Add directories to search path for required functions across servers
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for compute_rms_error.m, correct_unbalanced_bridge, 
    %   find_ind_str_in_cell.m, print_cellstr.m
    addpath(fullfile(functionsDirectory, 'Adams_Functions')); 
end

% Group the traces by unique vHold values, then average the grouped traces
dataCprGrouped = cell(nVHoldCond, 1);
dataCprAveraged = cell(nVHoldCond, 1);
for iVhold = 1:nVHoldCond
    % Get the current vHold value
    vnow = vUnique(iVhold);

    % Collect all cpr traces with this vHold value
    dCprGroupedThisVhold = dataCpr(vHoldCond == vnow);

    % Preallocate
    dCprAveragedThisVhold = zeros(ndpCpr, 4);

    % Take the time, current and conductance traces from the first trace
    dCprAveragedThisVhold(:, 1) = dCprGroupedThisVhold{1}(:, 1);
    dCprAveragedThisVhold(:, 3) = dCprGroupedThisVhold{1}(:, 3);
    dCprAveragedThisVhold(:, 4) = dCprGroupedThisVhold{1}(:, 4);

    % Average the voltage traces
    temp1 = cellfun(@(x) x(:, 2), dCprGroupedThisVhold, ...
                    'UniformOutput', false);
    vCprGroupedThisVhold = cell2mat(temp1');
    vCprAveragedThis = nanmean(vCprGroupedThisVhold, 2);
    dCprAveragedThisVhold(:, 2) = vCprAveragedThis;

    % Save in arrays
    dataCprGrouped{iVhold} = dCprGroupedThisVhold;
    dataCprAveraged{iVhold} = dCprAveragedThisVhold;
end

% Find unique vHold values
vUnique = force_column_numeric(sort(unique(vHoldCond)));

% Rename as new dataCpr
dataCpr = dataCprAveraged;

% Compute the baseline noise for the IPSC response
baseNoiseIpscr = compute_baseline_noise(vvecsMfrs, tvecsMfrs, ...
                                        baseWindowIpscrOrig);

% Find the baseline noise
baseNoiseCpr = zeros(nSwpsCpr, 1);
for iSwp = 1:nSwpsCpr
    % Find the baseline noise
    indBaselineCpr = nSamplesToPadForStabilization + (1:cprBaseLengthSamples);
    baselineCprThis = vCprAveragedThis(indBaselineCpr);
    cprBaselineNoise = compute_rms_error(baselineCprThis);

    baseNoiseCpr(iSwp) = cprBaselineNoise;
end

% Count the number of sweeps for the current pulse response
nSwpsCpr = numel(dataCpr);

% Find the holding potential
holdPotentialCpr = zeros(nSwpsCpr, 1);
for iSwp = 1:nSwpsCpr
    % Get the current voltage trace
    vCprAveragedThis = dataCpr{iSwp}(:, 2);

    % Find holding potential right before current pulse application
    mvind = 1:round(meanVoltageWindow/siMsOrig);   % indices for taking the mean of voltages
    baseInd = nSamplesToPadForStabilization + idxCpStart - ...
                round(0.5/siMsOrig) - fliplr(mvind);      % base indices
    cprBaselineVoltage = nanmean(vCprAveragedThis(baseInd));

    % Save in arrays
    holdPotentialCpr(iSwp) = cprBaselineVoltage;
end

%% Process IPSC responses from the median-filtered & resampled data vectors
% Print message
fprintf(['Processing IPSC responses from ', ...
            'the median-filtered & resampled data vectors ... \n']);

% Get the number of data points for the current pulse response
ndpCpr = size(dataCpr{1}, 1);

% Correct for traces that may have out-of-balance bridges
parfor iSwp = 1:nSwpsIpscr
    if isOutOfBalance(iSwp)
        % Get the old trace
        vvecOld = dataCpr{iSwp}(:, 2);
        ivecOld = dataCpr{iSwp}(:, 3);

        % Correct any unbalanced bridge in the trace
        vvecNew = correct_unbalanced_bridge(vvecOld, ivecOld);

        % Store the new trace
        dataCpr{iSwp}(:, 2) = vvecNew;
    end
end

function [dataCpr, dataIpscr, sweepInfoCpr, sweepInfoIpscr, dataCprAll] = ...
                m3ha_import_raw_traces (fileNames, swpInfo, initialSlopes, ...
                    cpStartWindowOrig, cprWinOrig, timeToStabilizeMs, ...
                    ipscTimeOrig, ipscDur, ...
                    epasEstimate, RinEstimate, ...
                    correctDcStepsFlag, oldAverageCprFlag, generateDataFlag, ...
                    varargin)

% Estimate the corresponding variations in holding current (nA)
if generateDataFlag
    holdCurrentNoiseIpscr = (baseNoiseIpscr / RinEstimate) * 5;
    holdCurrentNoiseCpr = (baseNoiseCpr / RinEstimate) * 5;
else
    holdCurrentNoiseIpscr = (baseNoiseIpscr / RinEstimate) * 0;
    holdCurrentNoiseCpr = (baseNoiseCpr / RinEstimate) * 0;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%