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
%
% Requires:
%       cd/argfun.m
%       cd/compute_default_sweep_info.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_time_vectors.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_numeric.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_row_count.m
%
% Used by:
%       cd/m3ha_run_neuron_once.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 2018-10-28 Now uses extract_subvectors.m
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
initSwpErrorDefault = [];   % no initial error values by default

% errorModeDefault = ;          % default TODO: Description of param1
% baseWindowDefault = ;       % default TODO: Description of param1
% baseNoiseDefault = ;        % default TODO: Description of param1

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

%% Preparation
% Count the number of samples
nSamples = count_samples(vSim);

% Set default time vector(s)
if isempty(tBoth)
    tBoth = create_time_vectors(nSamples);
end

% Compute default windows, noise and weights
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    compute_default_sweep_info(tBoth, vReal, ...
            'BaseWindow', baseWindow, 'FitWindow', fitWindow, ...
            'BaseNoise', baseNoise, 'SweepWeights', sweepWeights);

% Force data vectors to become column numeric vectors
[tBoth, vSim, vReal, iSim, iReal, fitWindow] = ...
    argfun(@force_column_numeric, tBoth, vSim, vReal, iSim, iReal, fitWindow);

% Force data arrays to become column cell arrays of column numeric vectors
[tBoth, vSim, vReal, iSim, iReal, fitWindow] = ...
    argfun(@force_column_cell, tBoth, vSim, vReal, iSim, iReal, fitWindow);

% Count the number of sweeps
nSweeps = count_vectors(vSim);

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

% Compute LTS errors
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
%     errors = merge_structs(swpErrors, )
    otherwise
        error('code logic error!');
end

% Combine errors
switch errorMode
    case 'SweepOnly'
        % Use the average sweep error as the total error
        totalError = swpErrors.avgSwpError;
    case 'Sweep&LTS'
%     errors = merge_structs(swpErrors, )
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
    fitWindow = force_row_numeric(fitWindow);
end
endPoints = find_window_endpoints(transpose(fitWindow), tBoth);

%       cd/force_row_numeric.m

%                   default == ones(nSweeps, 1)
% Set default sweep weights for averaging
if isempty(sweepWeights)
    sweepWeights = ones(nSweeps, 1);
end


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%