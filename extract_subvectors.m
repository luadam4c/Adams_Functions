function subVecs = extract_subvectors (vecs, varargin)
%% Extracts subvectors from vectors, given either endpoints, value windows or a certain align mode ('leftadjust', 'rightadjust')
% Usage: subVecs = extract_subvectors (vecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       subVecs1 = extract_subvectors({1:5, 2:6}, 'EndPoints', [1, 3])
%       subVecs2 = extract_subvectors({1:5, 2:6}, 'EndPoints', {[1, 3], [2, 4]})
%       subVecs3 = extract_subvectors(1:5, 'EndPoints', {[1, 3], [2, 4]})
%       subVecs4 = extract_subvectors({1:5, 2:6}, 'Windows', [2.5, 6.5])
% Outputs:
%       subVecs     - subvectors extracted
%                   specified as a numeric array 
%                       or a cell array of numeric vectors
% Arguments:
%       vecs        - vectors to extract
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'Indices': indices for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == create_indices(endPoints)
%                   - 'EndPoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], vecs)
%                   - 'Windows': value windows to extract 
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == []
%                   - 'AlignMethod': method for alignment/truncation
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'leftAdjust'  - equal number of elements from the left
%                       'rightAdjust' - equal number of elements from the right
%                       'none'        - no alignment/truncation
%                   default == 'leftAdjust'
%
% Requires:
%       cd/argfun.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/find_window_endpoints.m
%       cd/force_column_numeric.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_average_trace.m
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_common_directory.m
%       cd/find_passive_params.m
%       cd/filter_and_extract_pulse_response.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_run_neuron_once.m
%       cd/parse_pulse_response.m
%       cd/plot_traces.m

% File History:
% 2018-10-28 Created by Adam Lu
% 2018-10-29 Now returns empty if input is empty
% 2018-12-07 Now allows vecs to be a numeric array
% 2018-12-07 Now allows endPoints to be empty
% 2018-12-17 Now allows subvectors to be extracted from arbitrary indices
% 2018-12-17 Now allows subvectors to be extracted from an align mode
% TODO: check if all endpoints have 2 elements
% 

%% Hard-coded parameters
validAlignMethods = {'leftadjust', 'rightadjust', 'none'};

%% Default values for optional arguments
indicesDefault = [];            % set later
endPointsDefault = [];          % set later
windowsDefault = [];            % extract entire trace(s) by default
alignMethodDefault  = 'none';   % no alignment/truncation by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'Windows', windowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));

% Read from the Input Parser
parse(iP, vecs, varargin{:});
indices = iP.Results.Indices;
endPoints = iP.Results.EndPoints;
windows = iP.Results.Windows;
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);

% If indices is provided and endPoints or windows is also provided, 
%   display warning
if ~isempty(indices) && (~isempty(endPoints) || ~isempty(windows))
    fprintf('Windows will be ignored because end points are provided!\n');
end

% If endPoints is provided and windows is also provided, display warning
if ~isempty(endPoints) && ~isempty(windows)
    fprintf('Windows will be ignored because end points are provided!\n');
end

% TODO: check if all endpoints have 2 elements

%% Preparation
% If vecs is empty, return empty
if isempty(vecs)
    subVecs = vecs;
    return
end

% Force as column vectors
[endPoints, windows] = ...
    argfun(@(x) force_column_numeric(x, 'IgnoreNonVectors', true), ...
            endPoints, windows);

% Find end points if not provided
if isempty(endPoints)
    % TODO: Check if vecs are nondecreasing. 
    %   If not, ask for endPoints or indices to be passed

    % Extract the start and end indices of the vectors for fitting
    %   Note: If first argument is empty, 
    %           the first and last indices will be returned
    endPoints = find_window_endpoints(windows, vecs);
end

% Create indices if not provided
if isempty(indices)
    indices = create_indices(endPoints);
end

% If there is an align/truncation method used, apply it to indices
indices = align_subvectors(indices, alignMethod);

% If one of indices and vecs is a cell function, match the formats of 
%   indices and vecs so that cellfun can be used
if iscell(indices) || iscell(vecs)
    [indices, vecs] = ...
        match_format_vector_sets(indices, vecs, 'ForceCellOutputs', true);
end

%% Do the job
if iscell(vecs)
    subVecs = cellfun(@(x, y) extract_subvectors_helper(x, y), ...
                        vecs, indices, 'UniformOutput', false);
else
    subVecs = extract_subvectors_helper(vecs, indices);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subVec = extract_subvectors_helper (vec, indices)
%% Extract a subvector from vector(s) if not empty

% If the time window is out of range, or if the vector is empty, 
%   return an empty vector
if isempty(indices) || isempty(vec)
    subVec = [];
else
    subVec = vec(indices, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = align_subvectors(indices, alignMethod)

switch alignMethod
case {'leftadjust', 'rightadjust'}
    % Count the number of elements in each vector
    nSamples = count_samples(indices, 'ForceColumnOutput', true);

    % Get unique nSamples
    uniqueNSamples = unique(nSamples);

    % Get the number of unique nSamples
    nUniqueNSamples = numel(uniqueNSamples);

    % If nSamples are all equal, do nothing
    if nUniqueNSamples == 1
        return
    end

    % Get the minimum number of samples
    minNSamples = min(uniqueNSamples);

    % Create new endpoints
    switch alignMethod
    case 'leftadjust'
        % Create new endpoints by repetition
        newEndPoints = match_format_vector_sets([1, minNSamples], indices);
    case 'rightadjust'
        % Get the new starting positions 
        newStarts = nSamples - minNSamples + 1;

        % Create new endpoints
        newEndPoints = transpose([newStarts, nSamples]);
    end

    % Extract indices from new endpoints
    indices = extract_subvectors(indices, 'EndPoints', newEndPoints, ...
                                'AlignMethod', 'none');
case 'none'
    % Do nothing
otherwise
    error_unrecognized(get_var_name(alignMethod), alignMethod, mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumericvector(x) || iscellnumericvector(x), ...
                ['EndPoints must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));

[endPoints, vecs] = ...
    match_format_vector_sets(endPoints, vecs, 'ForceCellOutputs', false);

% If one of endPoints and vecs is a cell function, match the formats of 
%   endPoints and vecs so that cellfun can be used
if iscell(endPoints) || iscell(vecs)
    [endPoints, vecs] = ...
        match_format_vector_sets(endPoints, vecs, 'ForceCellOutputs', true);
end

if iscell(vecs)
    subVecs = cellfun(@(x, y) extract_subvectors_helper(x, y), ...
                        vecs, endPoints, 'UniformOutput', false);
else
    subVecs = extract_subvectors_helper(vecs, endPoints);
end

if isempty(endPoints) || isempty(vec)
    subVec = [];
else
    subVec = vec(endPoints(1):endPoints(2), :);
end

newEndPoints = repmat({[1, minNSamples]}, size(indices));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
