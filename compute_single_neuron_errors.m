function errors = compute_single_neuron_errors (vSim, vReal, varargin)
%% Computes all errors for single neuron data
% Usage: errors = compute_single_neuron_errors (vSim, vReal, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       errors      - a structure of all the errors computed, with fields:
%                       totalError
%                       fields returned by compute_sweep_errors.m
%                   specified as a scalar structure
% Arguments:    
%       vSim        - simulated voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       vReal       - recorded voltage traces
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
%                   - 'IvecsReal': recorded current traces
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == []
%                   - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == first half of the trace
%                   - 'FitWindow': time window to fit for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == second half of the trace
%                   - 'BaseNoise': baseline noise value(s)
%                   must be a numeric vector
%                   default == apply compute_baseline_noise.m
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 1 ./ baseNoise
%                   - 'NormalizeError': whether to normalize errors 
%                                       by an initial error
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'InitSwpError': initial sweep errors
%                   must be empty or a numeric vector with length == nSweeps
%                   default == []
%                   - 'IpscTime': current stimulation start time (ms), 
%                           Note: this is the time relative to which 
%                                       the peak delay is computed
%                   must be a positive scalar
%                   default == detected
%                   - 'FileBase': base of filename (without extension) 
%                                   corresponding to each vector
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == 'unnamed_1', 'unnamed_2', ...
%
% Requires:
%       cd/argfun.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_time_vectors.m
%       cd/decide_on_filebases.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_row_count.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 2018-10-28 Now uses extract_subvectors.m
% 2019-11-15 TODO: Add 'IpscTime' as an optional argument
% 

%% Hard-coded parameters
validErrorModes = {'SweepOnly', 'Sweep&LTS'};

%% Default values for optional arguments
errorModeDefault = 'SweepOnly'; %'Sweep&LTS'; % compute sweep & LTS errors by default
timeVecsDefault = [];           % set later
ivecsSimDefault = [];           % not provided by default
ivecsRealDefault = [];          % not provided by default
baseWindowDefault = [];         % set later
fitWindowDefault = [];          % set later
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set later
normalizeErrorDefault = false;  % don't normalize errors by default
initSwpErrorDefault = [];       % no initial error values by default
ipscTimeDefault = [];           % set later
ipscPeakWindowDefault = [];     % set later
fileBaseDefault = {};           % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vSim', ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['vSim must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'vReal', ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['vReal must be either a numeric vector ', ...
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
addParameter(iP, 'IvecsReal', ivecsRealDefault, ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['IvecsReal must be either a numeric vector ', ...
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
addParameter(iP, 'NormalizeError', normalizeErrorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'InitSwpError', initSwpErrorDefault, ...
    @(x) assert(isnumericvector(x), 'InitSwpError must be a numeric vector!'));
addParameter(iP, 'IpscTime', ipscTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'IpscPeakWindow', ipscPeakWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'FileBase', fileBaseDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['FileBase must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Read from the Input Parser
parse(iP, vSim, vReal, varargin{:});
errorMode = validatestring(iP.Results.ErrorMode, validErrorModes);
tBoth = iP.Results.TimeVecs;
iSim = iP.Results.IvecsSim;
iReal = iP.Results.IvecsReal;
baseWindow = iP.Results.BaseWindow;
fitWindow = iP.Results.FitWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;
normalizeError = iP.Results.NormalizeError;
initSwpError = iP.Results.InitSwpError;
ipscTime = iP.Results.IpscTime;
ipscPeakWindow = iP.Results.IpscPeakWindow;
fileBase = iP.Results.FileBase;

%% Preparation
% Count the number of samples
nSamples = count_samples(vSim);

% Set default time vector(s)
if isempty(tBoth)
    tBoth = create_time_vectors(nSamples, 'TimeUnits', 'ms');
end


% Compute default windows, noise and weights
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tBoth, vReal, ...
            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
            'BaseNoise', baseNoise, 'SweepWeights', sweepWeights);

% Force data vectors to become column numeric vectors
[tBoth, vSim, vReal, iSim, iReal, fitWindow] = ...
    argfun(@(x) force_column_vector(x, 'ForceCellOutput', true), ...
            tBoth, vSim, vReal, iSim, iReal, fitWindow);

% Count the number of sweeps
nSweeps = count_vectors(vSim);

% Create file bases if not provided
fileBase = decide_on_filebases(fileBase, nSweeps);

% Match row counts for sweep-dependent variables with the number of sweeps
[fitWindow, vReal, tBoth, iSim, iReal] = ...
    argfun(@(x) match_row_count(x, nSweeps), ...
            fitWindow, vReal, tBoth, iSim, iReal);

% Extract the start and end indices of the time vector for fitting
endPoints = find_window_endpoints(fitWindow, tBoth);

% Extract the regions to fit
[tBoth, vSim, vReal, iSim, iReal] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPoints), ...
            tBoth, vSim, vReal, iSim, iReal);

%% Compute errors
% Compute sweep errors
swpErrors = compute_sweep_errors(vSim, vReal, 'TimeVecs', tBoth, ...
                                'SweepWeights', sweepWeights, ...
                                'NormalizeError', normalizeError, ...
                                'InitSwpError', initSwpError);

% Compute pulse response feature errors
% TODO

% Compute IPSC response feature errors
% TODO
switch errorMode
    case 'SweepOnly'
        % Set as NaN for other errors
        ltsErrors.ltsAmpErrors = NaN;
        ltsErrors.ltsDelayErrors = NaN;
        ltsErrors.ltsSlopeErrors = NaN;
        ltsErrors.avgLtsAmpError = NaN;
        ltsErrors.avgLtsDelayError = NaN;
        ltsErrors.avgLtsSlopeError = NaN;
        ltsErrors.avgLtsError = NaN;
    case 'Sweep&LTS'
        % Analyze IPSC traces
        [ipscTableSim, ipscTableReal] = ...
            argfun(@(x) parse_ipsc(x, 'tVecs', tBoth, 'FileBase', fileBase, ...
                                    'StimStartMs', ipscTime, ...
                                    'PeakWindowMs', ipscPeakWindow), ...
                    iSim, iReal);

        % Find the peak delay of the IPSCs
        [ipscDelaySim, ipscDelayReal] = ...
            argfun(@(x) x.peakDelayMs, ipscTableSim, ipscTableReal);

        % Find and compute low-threshold spike features
        ltsFeaturesSim = parse_lts(vSim, 'MinPeakDelayMs', ipscDelaySim, ...
                            'StimStartMs', ipscTime, 'tVec0s', tBoth, ...
                            'FileBase', fileBase);
        ltsFeaturesReal = parse_lts(vReal, 'MinPeakDelayMs', ipscDelayReal, ...
                            'StimStartMs', ipscTime, 'tVec0s', tBoth, ...
                            'FileBase', fileBase);

        % TODO: compute_lts_errors.m
        % ltsErrors = compute_lts_errors(ltsFeaturesSim, ltsFeaturesReal, ...
        %                                 'SweepWeights', sweepWeights, ...
        %                                 'NormalizeError', normalizeError, ...
        %                                 'InitLtsError', initLtsError);
        % Set as NaN for other errors
        ltsErrors.ltsAmpErrors = NaN;
        ltsErrors.ltsDelayErrors = NaN;
        ltsErrors.ltsSlopeErrors = NaN;
        ltsErrors.avgLtsAmpError = NaN;
        ltsErrors.avgLtsDelayError = NaN;
        ltsErrors.avgLtsSlopeError = NaN;
        ltsErrors.avgLtsError = NaN;

        error('Not implemented yet!');

    otherwise
        error('code logic error!');
end

% Combine errors
switch errorMode
    case 'SweepOnly'
        % Use the average sweep error as the total error
        totalError = swpErrors.avgSwpError;
    case 'Sweep&LTS'
        %% TODO: Combine sweep and LTS errors
        totalError = swpErrors.avgSwpError;
        % totalError = combine_sweep_lts_errors();
    otherwise
        error('code logic error!');
end

%% Store in output errors structure
errors.totalError = totalError;
errors = merge_structs(errors, swpErrors);
errors = merge_structs(errors, ltsErrors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Make sure all windows, if a vector and not a matrix, are row vectors
if isvector(fitWindow)
    fitWindow = force_row_vector(fitWindow);
end
endPoints = find_window_endpoints(transpose(fitWindow), tBoth);

%       cd/force_row_vector.m

%                   default == ones(nSweeps, 1)
% Set default sweep weights for averaging
if isempty(sweepWeights)
    sweepWeights = ones(nSweeps, 1);
end

% Force data arrays to become column cell arrays of column numeric vectors
[tBoth, vSim, vReal, iSim, iReal, fitWindow] = ...
    argfun(@force_column_cell, tBoth, vSim, vReal, iSim, iReal, fitWindow);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
