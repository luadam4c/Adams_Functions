function [idxClosest, valClosest] = find_closest (vecs, target, varargin)
%% Finds the element(s) in numeric vector(s) closest to target(s)
% Usage: [idxClosest, valClosest] = find_closest (vecs, target, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [i, v] = find_closest(2:9, 5.6)
%       [i, v] = find_closest([1, 2, 3; 4, 5, 6; 7, 8, 9], 5)
%       [i, v] = find_closest([1, 2, 3; 4, 5, 6; 7, 8, 9], 5, 'Direction', 'down')
%       [i, v] = find_closest([1, 2, 3; 4, 5, 6; 7, 8, 9], 5, 'Direction', 'up')
%       [i, v] = find_closest([1, 2, 3; 4, 5, 6; 7, 8, 9], 5, 'Direction', 'none')
%       [i, v] = find_closest(9:-2:1, 4)
%       [i, v] = find_closest({5:-1:1, 8:-2:0}, 4)
%       [i, v] = find_closest(5:-1:1, [3.2, 3.8])
%       [i, v] = find_closest(5:-1:1, [3.2, 3.8], 'Direction', 'none')
%       [i, v] = find_closest(1:100000, 5489)
%       [i, v] = find_closest(9:-2:1, 11)
%
% Outputs:
%       idxClosest  - index(ices) of the closest value(s)
%                   specified as a positive integer vector
%       valClosest  - the closest value(s)
%                   specified as a numeric vector
%
% Arguments:
%       vecs        - vector(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       target      - target value(s)
%                   must be a numeric vector
%       varargin    - 'Direction': rounding direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'nearest'   - round to 'nearest'
%                       'down'      - always round down
%                       'up'        - always round up
%                       'none'      - no rounding, use interpolation
%                   default == 'nearest'
%
% Requires:
%       cd/argfun.m
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/isemptycell.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_gabab_conductance.m
%       cd/detect_spikes_multiunit.m
%       cd/parse_atf_swd.m
%       cd/parse_ipsc.m
%       cd/parse_phase_info.m
%       cd/parse_stim.m
%       vIRt/virt_moore.m

% File History:
% 2019-11-14 Created by Adam Lu
% 2019-11-25 Added 'none' as a direction
% 2020-01-03 Improved performance
% 2020-06-29 Now does not require imput to be monotonic
% 2021-05-07 Fixed bug

%% Hard-coded parameters
validDirections = {'nearest', 'down', 'up', 'none'};

%% Default values for optional arguments
directionDefault = 'nearest';

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
addRequired(iP, 'vecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array ', ...
                    'or a cell array of numeric vectors!']));
addRequired(iP, 'target', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));

% Read from the Input Parser
parse(iP, vecs, target, varargin{:});
direction = validatestring(iP.Results.Direction, validDirections);

%% Preparation
% Make sure target is a row vector
target = force_row_vector(target);

% Make sure vecs is not a row vector
vecs = force_column_vector(vecs, 'IgnoreNonVectors', true);

% Match the number of vecs with the number of targets
[vecs, target] = ...
    match_format_vector_sets(vecs, target, 'TreatRowVecAsOne', false);

% Sort the values in descending order
%   Note: Must be descending for 'nearest' to work
[vecs, origInd] = apply_to_all_cells(@(x) sort(x, 'descend'), vecs);

%% Do the job
% Create "time windows"
if iscell(target)
    windows = cellfun(@(x) [x; x], target, 'UniformOutput', false);
else
    windows = repmat(target, [2, 1]);
end

% Find endpoints that "include" the target value
indClosest = find_window_endpoints(windows, vecs, 'BoundaryMode', 'inclusive');

% If not found, the candidates are the endpoints of the vecs
if isemptycell(indClosest)
    indClosest = find_window_endpoints([], vecs);
end

% Find corresponding values
valsClosest = extract_subvectors(vecs, 'Indices', indClosest);

% Decide on the closest
if iscell(target)
    [idxClosestSorted, valClosest] = ...
        array_fun(@(x, y, z) find_closest_helper(x, y, z, direction), ...
                    target, valsClosest, indClosest);
else
    [idxClosestSorted, valClosest] = ...
        find_closest_helper(target, valsClosest, indClosest, direction);
end

% Force outputs as column vectors
[idxClosestSorted, valClosest] = ...
    argfun(@force_column_vector, idxClosestSorted, valClosest);

% Use original indices
if iscell(origInd)
    idxClosest = cellfun(@(x, y) x(y), origInd, num2cell(idxClosestSorted));
else
    idxClosest = origInd(idxClosestSorted);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxClosest, valClosest] = ...
                find_closest_helper (target, valsClosest, indClosest, direction)
% Sort the values in descending order
%   Note: Must be descending for 'nearest' to work
[valsClosest, origInd] = sort(valsClosest, 'descend');

% Reorder the indices in the same order
indClosest = extract_subvectors(indClosest, 'Indices', origInd);

% Compute the absolute differences
absDiffValues = abs(target - valsClosest);

% Choose the endpoint that is "closest"
switch direction
    case 'nearest'
        % For each column, find the row with minimum absolute difference
        [~, iChoice] = min(absDiffValues, [], 1);

        % Choose the row for each column
        [idxClosest, valClosest] = ...
            argfun(@(x) extract_elements(x, 'specific', 'Index', iChoice), ...
                    indClosest, valsClosest);
    case 'down'
        idxClosest = indClosest(end, :);
        valClosest = valsClosest(end, :);
    case 'up'
        idxClosest = indClosest(1, :);
        valClosest = valsClosest(1, :);
    case 'none'
        valClosest = target;
        idxClosest = interp1_custom(valsClosest, indClosest, valClosest);
    otherwise
        error('direction unrecognized!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vq = interp1_custom (x, v, xq)
%% TODO: Pull out as its own function

if isvector(x)
    vq = interp1_custom_helper(x, v, xq);
else
    vq = array_fun(@(a) interp1_custom_helper(x(:, a), v(:, a), xq(a)), ...
                    1:size(x, 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vq = interp1_custom_helper (x, v, xq)
%% Uses interp1 but allows v to have only one sample

uniqueV = unique(v);

if numel(uniqueV) == 1
    vq = uniqueV;
else
    vq = interp1(x, v, xq);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%