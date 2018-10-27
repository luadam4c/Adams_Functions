function errors = compute_single_neuron_errors (vSim, vReal, varargin)
%% Computes all errors for single neuron data
% Usage: errors = compute_single_neuron_errors (vSim, vReal, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       errors      - a structure of all the errors computed, with fields:
%                       swpError        - all sweep errors
%                       avgSwpError     - average sweep error
%                       normAvgSwpError - normalized average sweep error
%                       initSwpError    - initial sweep error
%                   specified as a scalar structure
% Arguments:    
%       vSim        - simulated voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       vReal       - recorded voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'TimeVecs': TODO: Description of TimeVecs
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsSim': TODO: Description of IvecsSim
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsReal': TODO: Description of IvecsReal
%                   must be a TODO
%                   default == TODO
%                   - 'SimMode': TODO: Description of SimMode
%                   must be a TODO
%                   default == TODO
%                   - 'FitWindow': TODO: Description of FitWindow
%                   must be a TODO
%                   default == TODO
%                   - 'SweepWeights': TODO: Description of SweepWeights
%                   must be a TODO
%                   default == TODO
%                   - 'NormalizeError': normalize errors by initial error
%                   must be a TODO
%                   default == TODO
%                   - 'InitSwpError': TODO: Description of BaseNoise
%                   must be a TODO
%                   default == TODO
% 
%                   - 'BaseWindow': TODO: Description of BaseWindow
%                   must be a TODO
%                   default == TODO
%                   - 'BaseNoise': TODO: Description of BaseNoise
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/argfun.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/compute_rms_error.m
%       cd/create_time_vectors.m
%       cd/find_window_endpoints.m
%       cd/force_row_numeric.m
%       cd/iscellnumericvector.m
%       cd/match_row_count.m
%       cd/normalize_by_initial_value.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeVecsDefault = [];       % set later
ivecsSimDefault = ;         % default TODO: Description of param1
ivecsRealDefault = ;        % default TODO: Description of param1
simModeDefault = ;          % default TODO: Description of param1
fitWindowDefault = ;        % default TODO: Description of param1
sweepWeightsDefault = ;     % default TODO: Description of param1
normalizeErrorDefault = false;  % don't normalize errors by default
initSwpErrorDefault = ;     % default TODO: Description of param1

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
    @(x) assert(isnumeric(x) && isvector(x) || iscellnumericvector(x), ...
                ['vSim must be either a numeric vector', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'vReal', ...
    @(x) assert(isnumeric(x) && isvector(x) || iscellnumericvector(x), ...
                ['vReal must be either a numeric vector', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVecs', timeVecsDefault, ...
    % TODO: validation function %);
addParameter(iP, 'IvecsSim', ivecsSimDefault, ...
    % TODO: validation function %);
addParameter(iP, 'IvecsReal', ivecsRealDefault, ...
    % TODO: validation function %);
addParameter(iP, 'SimMode', simModeDefault, ...
    % TODO: validation function %);
addParameter(iP, 'FitWindow', fitWindowDefault, ...
    % TODO: validation function %);
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    % TODO: validation function %);
addParameter(iP, 'NormalizeError', normalizeErrorDefault, ...
    % TODO: validation function %);
addParameter(iP, 'InitSwpError', initSwpErrorDefault, ...
    % TODO: validation function %);

% addParameter(iP, 'BaseWindow', baseWindowDefault, ...
%     % TODO: validation function %);
% addParameter(iP, 'BaseNoise', baseNoiseDefault, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vSim, vReal, varargin{:});
tBoth = iP.Results.TimeVecs;
iSim = iP.Results.IvecsSim;
iReal = iP.Results.IvecsReal;
simMode = iP.Results.SimMode;
fitWindow = iP.Results.FitWindow;
sweepWeights = iP.Results.SweepWeights;
normalizeError = iP.Results.NormalizeError;
initSwpError = iP.Results.InitSwpError;

% baseWindow = iP.Results.BaseWindow;
% baseNoise = iP.Results.BaseNoise;

% Make sure all windows, if a vector and not a matrix, are row vectors
if isvector(fitWindow)
    fitWindow = force_row_numeric(fitWindow);
end

%% Preparation
% Count the number of samples
nSamples = count_samples(vSim);

% Count the number of sweeps
nSweeps = count_vectors(vSim);

% Set default time vector(s)
tBoth = create_time_vectors(nSamples);

% Force data vectors to become column cell arrays of column numeric vectors
[tBoth, vSim, vReal, iSim, iReal] = ...
    argfun(@force_column_cell, tBoth, vSim, vReal, iSim, iReal);

% Match row counts for sweep-dependent variables with the number of sweeps
[fitWindow, tBoth, iSim, iReal] = ...
    argfun(@(x) match_row_count(x, nSweeps), fitWindow, tBoth, iSim, iReal);

% Simulation mode specific operations
switch simMode
    case 'passive'
    case 'active'
    otherwise
        error('simMode unrecognized!');
end

% Extract the lower and upper bounds of the fit windows
fitWindowLowerBounds = fitWindow(:, 1);
fitWindowUpperBounds = fitWindow(:, 2);

% Extract the start and end indices of the time vector for fitting
[idxStarts, idxEnds] = find_window_endpoints(tBoth);

% Extract the regions to fit
[tBoth, vSim, vReal, iSim, iReal] = ...
    argfun(@(w) cellfun(@(x, y, z) x(y:z), w, idxStarts, idxEnds), ...
            tBoth, vSim, vReal, iSim, iReal);

% Compute root-mean-square sweep errors over all sample points
swpError = cellfun(@(x, y) compute_rms_error(x, y), vSim, vReal);

% Compute a weighted root-mean-square error over all sweeps
avgSwpError = compute_weighted_average(swpError, 'Weights', sweepWeights, ...
                                        'AverageMethod', 'root-mean-square')

% If requested, make errors dimensionless by 
%   storing or dividing by an initial error value
if normalizeError
    [normAvgSwpError, initSwpError] = ...
        normalize_by_initial_value(avgSwpError, initSwpError);
end

%% Store in output errors structure
errors.swpError = swpError;
errors.avgSwpError = avgSwpError;
if normalizeError
    errors.normAvgSwpError = normAvgSwpError;
    errors.initSwpError = initSwpError;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fitWindowLowerBounds = match_vector_counts(fitWindowLowerBounds, [nSweeps, 1]);
fitWindowUpperBounds = match_vector_counts(fitWindowUpperBounds, [nSweeps, 1]);

if isempty(tBoth)
    tBoth = arrayfun(@(x) transpose(1:x), nSamples);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%