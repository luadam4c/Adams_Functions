function errorStruct = compute_sweep_errors (vSim, vReal, varargin)
%% Computes all errors for single neuron data
% Usage: errorStruct = compute_sweep_errors (vSim, vReal, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       errorStruct - a structure of all the errors computed, with fields:
%                       normalizeError
%                       normAvgSwpError - normalized average sweep error
%                       initSwpError    - initial sweep error
%                       avgSwpError     - average sweep error
%                       swpErrors       - all sweep errors
%                       fitWindow
%                       sweepWeights
%                   specified as a scalar structure
% Arguments:    
%       vSim        - simulated voltage traces
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric vector or a cell array of numeric vectors
%       vReal       - recorded voltage traces
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'TimeVecs': common time vectors
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric vector or a cell array of numeric vectors
%                   default == create_time_vectors(nSamples)
%                   - 'FitWindow': time window to fit for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == entire trace
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == set in compute_weighted_average.m
%                   - 'NormalizeError': whether to normalize errors 
%                                       by an initial error
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'InitSwpError': initial sweep errors
%                   must be empty or a numeric vector with length == nSweeps
%                   default == []
%   TODO:
%                   - 'ReturnSweepErrors': whether to return sweep errors only
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/argfun.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/compute_rms_error.m
%       cd/compute_weighted_average.m
%       cd/create_time_vectors.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_row_count.m
%       cd/normalize_by_initial_value.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/plot_fitted_traces.m
%       cd/m3ha_xolotl_plot.m

% File History:
% 2018-10-28 Adapted from compute_single_neuron_errors.m
% 2018-10-28 Now uses extract_subvectors.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeVecsDefault = [];       % set later
fitWindowDefault = [];      % use entire trace(s) by default
sweepWeightsDefault = [];   % set in compute_weighted_average.m
normalizeErrorDefault = false;  % don't normalize errors by default
initSwpErrorDefault = [];   % no initial error values by default

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
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['vSim must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'vReal', ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['vReal must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVecs', timeVecsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['TimeVecs must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'FitWindow', fitWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['FitWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeights must be a numeric vector!'));
addParameter(iP, 'NormalizeError', normalizeErrorDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'InitSwpError', initSwpErrorDefault, ...
    @(x) assert(isnumericvector(x), 'InitSwpError must be a numeric vector!'));

% Read from the Input Parser
parse(iP, vSim, vReal, varargin{:});
tBoth = iP.Results.TimeVecs;
fitWindow = iP.Results.FitWindow;
sweepWeights = iP.Results.SweepWeights;
normalizeError = iP.Results.NormalizeError;
initSwpError = iP.Results.InitSwpError;

%% Preparation
% Count the number of samples
nSamples = count_samples(vSim);

% Set default time vector(s)
if isempty(tBoth)
    tBoth = create_time_vectors(nSamples, 'TimeUnits', 'ms');
end

% Force data vectors to become column cell arrays of column numeric vectors
[tBoth, vSim, vReal, fitWindow] = ...
    argfun(@(x) force_column_vector(x, 'ForceCellOutput', true), ...
            tBoth, vSim, vReal, fitWindow);

% Count the number of sweeps
nSweeps = count_vectors(vSim);

% Match row counts for sweep-dependent variables with the number of sweeps
[fitWindow, vReal, tBoth] = ...
    argfun(@(x) match_row_count(x, nSweeps), fitWindow, vReal, tBoth);

% Extract the start and end indices of the time vector for fitting
endPoints = find_window_endpoints(fitWindow, tBoth);

% Extract the regions to fit
[vSim, vReal] = ...
    argfun(@(x) extract_subvectors(x, 'Endpoints', endPoints), vSim, vReal);

%% Compute errors
% Compute root-mean-square sweep errors over all sample points
swpErrors = compute_rms_error(vSim, vReal);

% Compute a weighted root-mean-square error over all sweeps
avgSwpError = compute_weighted_average(swpErrors, 'Weights', sweepWeights, ...
                        'IgnoreNaN', true, 'AverageMethod', 'root-mean-square');

% If requested, make errors dimensionless by 
%   storing or dividing by an initial error value
if normalizeError
    [normAvgSwpError, initSwpError] = ...
        normalize_by_initial_value(avgSwpError, initSwpError);
end

%% Store in output errors structure
if normalizeError
    errorStruct.normalizeError = normalizeError;
    errorStruct.normAvgSwpError = normAvgSwpError;
    errorStruct.initSwpError = initSwpError;
end
errorStruct.avgSwpError = avgSwpError;
errorStruct.swpErrors = swpErrors;
errorStruct.fitWindow = fitWindow;
errorStruct.sweepWeights = sweepWeights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fitWindowLowerBounds = match_vector_counts(fitWindowLowerBounds, [nSweeps, 1]);
fitWindowUpperBounds = match_vector_counts(fitWindowUpperBounds, [nSweeps, 1]);

if isempty(tBoth)
    tBoth = arrayfun(@(x) transpose(1:x), nSamples);
end

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

[idxStarts, idxEnds] = find_window_endpoints(transpose(fitWindow), tBoth);

[tBoth, vSim, vReal, iSim, iReal] = ...
    argfun(@(w) cellfun(@(x, y, z) x(y:z), w, ...
                        num2cell(idxStarts), num2cell(idxEnds)), ...
            tBoth, vSim, vReal, iSim, iReal);

[tBoth, vSim, vReal, iSim, iReal] = ...
    argfun(@(w) if isempty(w); w; else; ...
            cellfun(@(x, y) x(y(1):y(2), w, endPoints)); end, ...
            tBoth, vSim, vReal, iSim, iReal);

swpErrors = cellfun(@(x, y) compute_rms_error(x, y), vSim, vReal);

%       cd/force_row_vector.m
% Make sure all windows, if a vector and not a matrix, are row vectors
if isvector(fitWindow)
    fitWindow = force_row_vector(fitWindow);
end
endPoints = find_window_endpoints(transpose(fitWindow), tBoth);

% Force data arrays to become column cell arrays of column numeric vectors
[tBoth, vSim, vReal, fitWindow] = ...
    argfun(@force_column_cell, tBoth, vSim, vReal, fitWindow);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
