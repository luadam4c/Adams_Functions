function errors = compute_single_neuron_errors (vSim, vRec, varargin)
%% Computes the average total error for a single neuron
% Usage: errors = compute_single_neuron_errors (vSim, vRec, varargin)
% Explanation:
%       Let a = (1 / (1 + lts2SweepErrorRatio))
%           b = (1 / (1 + match2FeatureErrorRatio))
%           c = featureWeights(1)/sum(featureWeights)
%           d = featureWeights(2)/sum(featureWeights)
%           e = featureWeights(3)/sum(featureWeights)
%
%       Then:
%       totalError = ...
%           avgSwpError * a + ...
%           ltsMatchError * (1 - a) * (1 - b) + ...
%           avgLtsAmpError * (1 - a) * b * c + ...
%           avgLtsDelayError * (1 - a) * b * d + ...
%           avgLtsSlopeError * (1 - a) * b * e
%
%       Note: If any of the 5 types of errors is NaN, 
%               the weights are renormalized across the remaining values
%           Therefore, in some cases the following will not be true
%       totalError = avgSwpError * a + avgLtsError * (1 - a)
%
%
% Example(s):
%       TODO
%
% Outputs:
%       errors      - a structure of all the errors computed, with fields:
%                       totalError
%                       errorsToAverage
%                       errorWeights
%                       fields returned by compute_sweep_errors.m
%                       fields returned by compute_lts_errors.m
%                       baseWindow
%                       fitWindow
%                       baseNoise
%                       sweepWeights
%                   specified as a scalar structure
%
% Arguments:    
%       vSim        - simulated voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       vRec        - recorded voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'ErrorMode': error mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'SweepOnly' - compute sweep errors only
%                       'Sweep&LTS' - compute sweep & LTS errors only
%                   default == 'passive'
%                   - 'TimeVecs': common time vectors
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == create_time_vectors(nSamples)
%                   - 'IvecsSim': simulated current traces
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == []
%                   - 'IvecsRec': recorded current traces
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == []
%                   - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set in compute_default_sweep_info.m
%                   - 'FitWindow': time window to fit for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == set in compute_default_sweep_info.m
%                   - 'BaseNoise': baseline noise value(s)
%                   must be a numeric vector
%                   default == set in compute_default_sweep_info.m
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_default_sweep_info.m
%                   - 'LtsFeatureWeights': LTS feature weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == [1, 1, 1]
%                   - 'MissedLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the existence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   - 'FalseLtsError': a dimensionless error that penalizes 
%                                       a misprediction of the absence of LTS
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_lts_errors.m
%                   - 'Match2FeatureErrorRatio': ratio of LTS match error to 
%                                                   LTS feature error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 1
%                   - 'Lts2SweepErrorRatio': ratio of LTS error to sweep error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 3
%                   - 'NormalizeError': whether to normalize errors 
%                                       by an initial error
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'InitSwpError': initial average sweep error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == []
%                   - 'InitLtsError': initial average LTS error
%                   must be empty or a numeric vector with length == nSweeps
%                   default == []
%                   - 'IpscTime': current stimulation start time (ms), 
%                           Note: this is the time relative to which 
%                                       the peak delay is computed
%                   must be a positive scalar
%                   default == set in parse_ipsc.m
%                   - 'IpscPeakWindow': window (ms) to look for IPSC peak
%                   must be a numeric vector
%                   default == set in parse_ipsc.m
%                   - 'FileBase': base of filename (without extension) 
%                                   corresponding to each vector
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == set in decide_on_filebases.m
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix to prepend to output file names
%                   must be a character array
%                   default == set in parse_ipsc.m and parse_lts.m
%                   - 'SaveLtsInfoFlag': whether to save LTS info
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveLtsStatsFlag': whether to save LTS statistics
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotIpeakFlag': whether to current peak analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotLtsFlag': whether to plot vtrace/LTS/burst analyses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PlotStatisticsFlag': whether to plot LTS statistics
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/argfun.m
%       cd/compute_default_sweep_info.m
%       cd/compute_lts_errors.m
%       cd/compute_sweep_errors.m
%       cd/compute_weighted_average.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_time_vectors.m
%       cd/decide_on_filebases.m
%       cd/extract_common_prefix.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_row_count.m
%       cd/merge_structs.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/set_default_flag.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 2018-10-28 Now uses extract_subvectors.m
% 2019-11-15 Added 'IpscTime' as an optional argument
% 2019-11-15 Added 'IpscPeakWindow' as an optional argument
% 2019-11-15 Added 'FileBase' as an optional argument
% 2019-11-16 Added 'LtsFeatureWeights' as an optional argument
% 2019-11-16 Added 'LtsExistError' as an optional argument
% 2019-11-16 Added 'Lts2SweepErrorRatio' as an optional argument
% 2019-11-17 Added 'OutFolder' as an optional parameter
% 2019-11-18 Added 'Prefix' as an optional parameter
% 2019-12-18 Added 'Match2FeatureErrorRatio' as an optional parameter
% 2019-12-18 Added errorWeights in output
% 2019-12-19 Fixed errorWeights when LTS is not present or all mismatched
% 2019-12-22 Now sets all default errorWeights
% 2020-01-03 Now uses fitWindow as LTS search window
% TODO:
%   Implement saveLtsStatsFlag
%   Implement plotStatisticsFlag
%   Fix plotIpeakFlag

%% Hard-coded parameters
validErrorModes = {'SweepOnly', 'Sweep&LTS'};

% Consistent with singleneuronfitting71.m
defaultLtsFeatureWeights = [1; 1; 1];   % default weights for optimizing 
                                        %   LTS statistics
defaultMatch2FeatureErrorRatio = 1;     % default error ratio between
                                        %   match error and avg feature error
defaultLts2SweepErrorRatio = 3;         % default error ratio of LTS error 
                                        %   to sweep error



%% Default values for optional arguments
errorModeDefault = 'SweepOnly'; %'Sweep&LTS'; % compute sweep & LTS errors by default
timeVecsDefault = [];           % set later
ivecsSimDefault = [];           % not provided by default
ivecsRecDefault = [];           % not provided by default
baseWindowDefault = [];         % set later
fitWindowDefault = [];          % set later
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set later
ltsFeatureWeightsDefault = [];  % set later
missedLtsErrorDefault = [];     % set in compute_lts_errors.m
falseLtsErrorDefault = [];      % set in compute_lts_errors.m
match2FeatureErrorRatioDefault = [];    % set later
lts2SweepErrorRatioDefault = [];        % set later
normalizeErrorDefault = false;  % don't normalize errors by default
initSwpErrorDefault = [];       % no initial sweep error values by default
initLtsErrorDefault = [];       % no initial LTS error values by default
ipscTimeDefault = [];           % set later
ipscPeakWindowDefault = [];     % set later
fileBaseDefault = {};           % set later
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
prefixDefault = '';             % set later
saveLtsInfoFlagDefault = false;  % don't save LTS info by default
saveLtsStatsFlagDefault = false; % don't save LTS statistics by default
plotIpeakFlagDefault = false;
plotLtsFlagDefault = false;
plotStatisticsFlagDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vSim', ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['vSim must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'vRec', ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['vRec must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ErrorMode', errorModeDefault, ...
    @(x) any(validatestring(x, validErrorModes)));
addParameter(iP, 'TimeVecs', timeVecsDefault, ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['TimeVecs must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'IvecsSim', ivecsSimDefault, ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['IvecsSim must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'IvecsRec', ivecsRecDefault, ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['IvecsRec must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
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
addParameter(iP, 'LtsFeatureWeights', ltsFeatureWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'LtsFeatureWeights must be a numeric vector!'));
addParameter(iP, 'MissedLtsError', missedLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'MissedLtsError must be a numeric vector!'));
addParameter(iP, 'FalseLtsError', falseLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'FalseLtsError must be a numeric vector!'));
addParameter(iP, 'Match2FeatureErrorRatio', match2FeatureErrorRatioDefault, ...
    @(x) assert(isnumericvector(x), 'Match2FeatureErrorRatio must be a numeric vector!'));
addParameter(iP, 'Lts2SweepErrorRatio', lts2SweepErrorRatioDefault, ...
    @(x) assert(isnumericvector(x), 'InitLtsError must be a numeric vector!'));
addParameter(iP, 'NormalizeError', normalizeErrorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'InitSwpError', initSwpErrorDefault, ...
    @(x) assert(isnumericvector(x), 'InitSwpError must be a numeric vector!'));
addParameter(iP, 'InitLtsError', initLtsErrorDefault, ...
    @(x) assert(isnumericvector(x), 'InitLtsError must be a numeric vector!'));
addParameter(iP, 'IpscTime', ipscTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'IpscPeakWindow', ipscPeakWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileBase must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveLtsInfoFlag', saveLtsInfoFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveLtsStatsFlag', saveLtsStatsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotIpeakFlag', plotIpeakFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotLtsFlag', plotLtsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotStatisticsFlag', plotStatisticsFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vSim, vRec, varargin{:});
errorMode = validatestring(iP.Results.ErrorMode, validErrorModes);
tBoth = iP.Results.TimeVecs;
iSim = iP.Results.IvecsSim;
iRec = iP.Results.IvecsRec;
baseWindow = iP.Results.BaseWindow;
fitWindow = iP.Results.FitWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;
ltsFeatureWeights = iP.Results.LtsFeatureWeights;
missedLtsError = iP.Results.MissedLtsError;
falseLtsError = iP.Results.FalseLtsError;
match2FeatureErrorRatio = iP.Results.Match2FeatureErrorRatio;
lts2SweepErrorRatio = iP.Results.Lts2SweepErrorRatio;
normalizeError = iP.Results.NormalizeError;
initSwpError = iP.Results.InitSwpError;
initLtsError = iP.Results.InitLtsError;
ipscTime = iP.Results.IpscTime;
ipscPeakWindow = iP.Results.IpscPeakWindow;
fileBase = iP.Results.FileBase;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;
saveLtsInfoFlag = iP.Results.SaveLtsInfoFlag;
saveLtsStatsFlag = iP.Results.SaveLtsStatsFlag;
plotIpeakFlag = iP.Results.PlotIpeakFlag;
plotLtsFlag = iP.Results.PlotLtsFlag;
plotStatisticsFlag = iP.Results.PlotStatisticsFlag;

%% Preparation
% Decide on LTS feature weights
if isempty(ltsFeatureWeights)
    ltsFeatureWeights = defaultLtsFeatureWeights;
end

% Decide on LTS match to feature error ratio
if isempty(match2FeatureErrorRatio)
    match2FeatureErrorRatio = defaultMatch2FeatureErrorRatio;
end

% Decide on LTS to sweep error ratio
if isempty(lts2SweepErrorRatio)
    lts2SweepErrorRatio = defaultLts2SweepErrorRatio;
end

% Determine whether LTS errors are computed
computeLtsError = set_default_flag([], strcmp(errorMode, 'Sweep&LTS') && ...
                                    lts2SweepErrorRatio ~= 0);

% Determine whether voltage traces needs to be parsed
parseVoltage = set_default_flag([], strcmp(errorMode, 'Sweep&LTS') && ...
                        (computeLtsError || saveLtsInfoFlag || ...
                        saveLtsStatsFlag || plotLtsFlag || plotStatisticsFlag));

% Determine whether current traces needs to be parsed
parseCurrent = set_default_flag([], strcmp(errorMode, 'Sweep&LTS') && ...
                                (parseVoltage || plotIpeakFlag));

% Count the number of samples
nSamples = count_samples(vSim);

% Set default time vector(s)
if isempty(tBoth)
    tBoth = create_time_vectors(nSamples, 'TimeUnits', 'ms');
end

% Compute default windows, noise and weights
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tBoth, vRec, ...
            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
            'BaseNoise', baseNoise, 'SweepWeights', sweepWeights);

% Force data vectors to become column numeric vectors
[tBoth, vSim, vRec, iSim, iRec, fitWindow] = ...
    argfun(@(x) force_column_vector(x, 'ForceCellOutput', true), ...
            tBoth, vSim, vRec, iSim, iRec, fitWindow);

% Count the number of sweeps
nSweeps = count_vectors(vSim);

% Create file bases if not provided
fileBase = decide_on_filebases(fileBase, nSweeps);

% Decide on prefix if not provided
if isempty(prefix)
    prefix = extract_common_prefix(fileBase);
end

% Match row counts for sweep-dependent variables with the number of sweeps
[fitWindow, vRec, tBoth, iSim, iRec] = ...
    argfun(@(x) match_row_count(x, nSweeps), ...
            fitWindow, vRec, tBoth, iSim, iRec);

% Create stuff for simulated and recorded separately
if parseCurrent || parseVoltage
    % Create differenct file bases
    [fileBaseRec, fileBaseSim] = ...
        argfun(@(x) strcat(fileBase, x), '_rec', '_sim');

    % Create differenct prefixes
    [prefixRec, prefixSim] = ...
        argfun(@(x) strcat(prefix, x), '_rec', '_sim');
end

%% Parse current traces
if parseCurrent
    % Analyze IPSC traces
    ipscTableSim = parse_ipsc(iSim, 'tVecs', tBoth, ...
                'StimStartMs', ipscTime, 'PeakWindowMs', ipscPeakWindow, ...
                'Verbose', false, 'PlotFlag', plotIpeakFlag, ...
                'FileBase', fileBaseSim, 'OutFolder', outFolder, ...
                'Prefix', prefixSim);
    ipscTableRec = parse_ipsc(iRec, 'tVecs', tBoth, ...
                'StimStartMs', ipscTime, 'PeakWindowMs', ipscPeakWindow, ...
                'Verbose', false, 'PlotFlag', plotIpeakFlag, ...
                'FileBase', fileBaseRec, 'OutFolder', outFolder, ...
                'Prefix', prefixRec);
end

%% Parse voltage traces
if parseVoltage
    % Find the peak delay of the IPSCs
    [ipscDelaySim, ipscDelayRec] = ...
        argfun(@(x) x.peakDelayMs, ipscTableSim, ipscTableRec);

    % Find and compute low-threshold spike features
    ltsFeaturesSim = parse_lts(vSim, 'MinPeakDelayMs', ipscDelaySim, ...
                        'StimStartMs', ipscTime, 'tVec0s', tBoth, ...
                        'NoiseWindowMsOrMaxNoise', baseWindow, ...
                        'Verbose', false, 'PlotFlag', plotLtsFlag, ...
                        'SaveSheetFlag', saveLtsInfoFlag, ...
                        'FileBase', fileBaseSim, 'OutFolder', outFolder, ...
                        'Prefix', prefixSim, 'SearchWindowMs', fitWindow);
    ltsFeaturesRec = parse_lts(vRec, 'MinPeakDelayMs', ipscDelayRec, ...
                        'StimStartMs', ipscTime, 'tVec0s', tBoth, ...
                        'NoiseWindowMsOrMaxNoise', baseWindow, ...
                        'Verbose', false, 'PlotFlag', plotLtsFlag, ...
                        'SaveSheetFlag', saveLtsInfoFlag, ...
                        'FileBase', fileBaseRec, 'OutFolder', outFolder, ...
                        'Prefix', prefixRec, 'SearchWindowMs', fitWindow);
end

%% Compute errors
% Compute sweep errors
swpErrors = compute_sweep_errors(vSim, vRec, 'TimeVecs', tBoth, ...
                                'FitWindow', fitWindow, ...
                                'SweepWeights', sweepWeights, ...
                                'NormalizeError', normalizeError, ...
                                'InitSwpError', initSwpError);

% Compute pulse response feature errors
% TODO

% Compute IPSC response feature errors
if computeLtsError
    % Compute low-threshold spike errors
    ltsErrors = compute_lts_errors(ltsFeaturesSim, ltsFeaturesRec, ...
                                    'BaseNoise', baseNoise, ...
                                    'SweepWeights', sweepWeights, ...
                                    'FeatureWeights', ltsFeatureWeights, ...
                                    'MissedLtsError', missedLtsError, ...
                                    'FalseLtsError', falseLtsError, ...
                                    'NormalizeError', normalizeError, ...
                                    'InitLtsError', initLtsError);
else
    % Set as NaN for other errors
    ltsErrors.ltsAmpErrors = NaN;
    ltsErrors.ltsDelayErrors = NaN;
    ltsErrors.ltsSlopeErrors = NaN;
    ltsErrors.ltsSweepWeights = NaN;
    ltsErrors.avgLtsAmpError = NaN;
    ltsErrors.avgLtsDelayError = NaN;
    ltsErrors.avgLtsSlopeError = NaN;
    ltsErrors.avgLtsFeatureError = NaN;
    ltsErrors.ltsMatchError = NaN;
    ltsErrors.avgLtsError = NaN;
end

%% Compute normalized error weights
% Combine the errors
errorsToAverage = [swpErrors.avgSwpError; ltsErrors.ltsMatchError; ...
                    ltsErrors.avgLtsAmpError; ltsErrors.avgLtsDelayError; ...
                    ltsErrors.avgLtsSlopeError];

% Compute coefficients
%   Note:
%       Let a = (1 / (1 + lts2SweepErrorRatio))
%           b = (1 / (1 + match2FeatureErrorRatio))
%           c = featureWeights(1)/sum(featureWeights)
%           d = featureWeights(2)/sum(featureWeights)
%           e = featureWeights(3)/sum(featureWeights)
a = (1 / (1 + lts2SweepErrorRatio));
b = (1 / (1 + match2FeatureErrorRatio));
c = ltsFeatureWeights(1)/sum(ltsFeatureWeights);
d = ltsFeatureWeights(2)/sum(ltsFeatureWeights);
e = ltsFeatureWeights(3)/sum(ltsFeatureWeights);

% Compute weights
%   Note: Assuming all errors are not NaN, we have
%       totalError = ...
%           avgSwpError * a + ...
%           ltsMatchError * (1 - a) * (1 - b) + ...
%           avgLtsAmpError * (1 - a) * b * c + ...
%           avgLtsDelayError * (1 - a) * b * d + ...
%           avgLtsSlopeError * (1 - a) * b * e
errorWeights = [a; (1 - a) * (1 - b); ...
                (1 - a) * b * c; ...
                (1 - a) * b * d; ...
                (1 - a) * b * e];

% Renormalize weights if some are NaN
if any(isnan(errorsToAverage))
    errorWeights(isnan(errorsToAverage)) = 0;
    errorWeights = errorWeights ./ sum(errorWeights);
end

% Total error (dimensionless) is the weighted average of sweep error, 
%   LTS match error, LTS amp error, LTS time error and LTS slope error
totalError = compute_weighted_average(errorsToAverage, 'IgnoreNan', true, ...
                        'Weights', errorWeights, 'AverageMethod', 'linear');

%% Store in output errors structure
errors.totalError = totalError;
errors.errorsToAverage = errorsToAverage;
errors.errorWeights = errorWeights;
errors.lts2SweepErrorRatio = lts2SweepErrorRatio;
errors = merge_structs(errors, ltsErrors);
errors = merge_structs(errors, swpErrors);
errors.baseWindow = baseWindow;
errors.fitWindow = fitWindow;
errors.baseNoise = baseNoise;
errors.sweepWeights = sweepWeights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract the start and end indices of the time vector for fitting
endPoints = find_window_endpoints(fitWindow, tBoth);

% Extract the regions to fit
[tBoth, vSim, vRec, iSim, iRec] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPoints), ...
            tBoth, vSim, vRec, iSim, iRec);

% Combine the errors and weights
errorsToAverage = [swpErrors.avgSwpError; ltsErrors.avgLtsError];
weightsForErrors = [1; lts2SweepErrorRatio];

% Total error (dimensionless) is the weighted average of sweep error 
%   and LTS error, weighted by lts2SweepErrorRatio
totalError = compute_weighted_average(errorsToAverage, 'IgnoreNan', true, ...
                        'Weights', weightsForErrors, 'AverageMethod', 'linear');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
