function rmsErrors = compute_rms_error(vec1s, varargin)
%% Computes the root mean squared error(s) given one or two sets of vectors
% Usage: rmsErrors = compute_rms_error(vec1s, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       rmsErrors1 = compute_rms_error(1:5)
%       rmsErrors2 = compute_rms_error(1:5, 2:6)
%       rmsErrors3 = compute_rms_error({1:5, 2:6}, 'Endpoints', [1, 3])
%       rmsErrors4 = compute_rms_error({1:5, 2:6}, 'Endpoints', {[1, 3], [2, 4]})
%       rmsErrors2 = compute_rms_error({1:5, 2:4}, {2:6, 3:5})
%       rmsErrors5 = compute_rms_error(1:5, 2:6, 'Endpoints', {[1, 3], [2, 4]})
%       rmsErrors6 = compute_rms_error(magic(3), magic(3) + 1)
%       rmsErrors7 = compute_rms_error(magic(3))
%
% Outputs:
%       rmsErrors   - root mean squared error(s) for each vector
%                       or between each pair of vectors
%                   specified as a numeric vector
%
% Arguments:
%       vec1s       - the first set of vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       vec2s       - (opt) the second set of vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == nanmean(vec1s)
%       varargin    - 'Endpoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], vec1s)
%                   - 'IgnoreNan': whether to ignore NaN entries
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Weights': weights for averaging
%                   must be a numeric array
%                   default == set in compute_weighted_average.m
%                   - Any other parameter-value pair for 
%                           compute_weighted_average()
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/compute_weighted_average.m
%       cd/iscellnumeric.m
%       cd/iscellnumericvector.m
%       cd/isemptycell.m
%       cd/extract_subvectors.m
%       cd/force_matrix.m
%       cd/force_column_vector.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_baseline_noise.m
%       cd/compute_lts_errors.m
%       cd/compute_sweep_errors.m
%       cd/m3ha_simulate_population.m

% File History:
% 2018-07-09 Modified from the built-in rms.m
% 2018-10-28 Now vectors do not need to have equal lengths
% 2018-10-28 Now takes multiple vectors as arguments
% 2018-10-28 Added 'Endpoints' as an optional argument
% 2019-11-17 Added 'IgnoreNan' and 'Weights' as optional arguments
% 2020-01-01 Now forces vectors as a matrix instead

%% Default values for optional arguments
vecs2Default = [];              % set later
endPointsDefault = [];          % set in extract_subvectors.m
ignoreNanDefault = true;        % ignores NaN by default
weightsDefault = [];            % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec1s', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'vec2s', vecs2Default, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec2s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumericvector(x), ...
                ['EndPoints must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));
addParameter(iP, 'IgnoreNan', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Weights', weightsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'3d'}));

% Read from the Input Parser
parse(iP, vec1s, varargin{:});
vec2s = iP.Results.vec2s;
endPoints = iP.Results.EndPoints;
ignoreNan = iP.Results.IgnoreNan;
valueWeights = iP.Results.Weights;

% Keep unmatched arguments for the compute_weighted_average() function
otherArguments = iP.Unmatched;

%% Preparation
% Force row vectors as column vectors
[vec1s, vec2s] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), vec1s, vec2s);

% Force vec1s and vec2s to be numeric matrices,
%   unless the number of samples in each vector are different
if iscell(vec1s) || iscell(vec2s)
    [vec1s, vec2s] = ...
        argfun(@(x) force_matrix(x, 'AlignMethod', 'none', ...
                                'Verbose', false), ...
                vec1s, vec2s);
end

% If not provided, set default vec2s to be the means of each vector
%   Otherwise, restrict to the same endpoints
vec2s = set_vec2_if_empty(vec1s, vec2s);

% Restrict to the given end points if needed
if ~isempty(endPoints)
    [vec1s, vec2s] = ...
        argfun(@(x) extract_subvectors(x, 'Endpoints', endPoints), ...
                vec1s, vec2s);
end

%% Do the job
if iscell(vec1s) || iscell(vec2s)
    % Make sure vec1s and vec2s are both column cell arrays
    %   with the same number of vectors
    [vec1s, vec2s] = match_format_vector_sets(vec1s, vec2s, ...
                    'ForceCellOutput', true);

    % Compute the root-mean_square error for each pair of vectors
    rmsErrors = array_fun(@(x, y) compute_rms_error_helper(x, y, ...
                                ignoreNan, valueWeights, otherArguments), ...
                            vec1s, vec2s);
else
    rmsErrors = compute_rms_error_helper(vec1s, vec2s, ignoreNan, ...
                                        valueWeights, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vec2 = set_vec2_if_empty (vec1, vec2)
%% Set the default second vector to be the mean values of the first vector

if isempty(vec2) || iscell(vec2) && all(all(isemptycell(vec2)))
    if iscell(vec1)
        vec2 = array_fun(@(x, y) set_vec2_if_empty(x, y), ...
                        vec1, vec2, 'UniformOutput', false);
    else
        % Count the number of rows
        nRows = size(vec1, 1);

        % Replace each column with the nanmean of that column
        vec2 = repmat(nanmean(vec1, 1), nRows, 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmsError = compute_rms_error_helper (vec1, vec2, ignoreNan, ...
                                                valueWeights, otherArguments)
%% Compute the root-mean-square error between two vectors

% Compute errors at every sample point
errors = vec1 - vec2;

% Compute the root-mean-squared error
rmsError = compute_weighted_average(errors, 'IgnoreNan', ignoreNan, ...
                                    'Weights', valueWeights, ...
                                    'AverageMethod', 'root-mean-square', ...
                                    otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Check relationships between arguments
if length(vec1s) ~= length(vec2s)
    error('The two vectors must have the same length!');
end

% Compare each sample point to the mean of the entire vector
vec2s = nanmean(vec1s) * ones(size(vec1s));

% Compare each sample point to the mean of the entire vector
if iscell(vec1s)
    vec2s = cellfun(@nanmean, vec1s, 'UniformOutput', false);
else
    vec2s = nanmean(vec1s);
end

% Make sure vec2s is either a column vector or a cell array of column vectors
vec2s = force_column_vector(vec2s);

if iscell(vec1s)
else
    rmsErrors = compute_rms_error_helper(vec1s, vec2s);
end

%   TODO: implement dim as in rms.m

if isempty(vec2s)
    vec2s = cellfun(@nanmean, vec1s, 'UniformOutput', false);
else
    vec2s = extract_subvectors(vec2s, 'Endpoints', endPoints);
end

[vec1s, vec2s] = match_format_vector_sets(vec1s, vec2s, 'ForceCellOutputs', true);

% Compute the squared error
squaredError = errors .* conj(errors);

% Compute the mean-squared error
meanSquaredError = nanmean(squaredError);

% Compute the root-mean-squared error
rmsError = sqrt(meanSquaredError);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
