function subVecs = extract_subvectors (vecs, varargin)
%% Extracts subvectors from vectors, given either endpoints, value windows or a certain align mode ('leftAdjust', 'rightAdjust')
% Usage: subVecs = extract_subvectors (vecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       subVecs1 = extract_subvectors({1:5, 2:6}, 'EndPoints', [1, 3])
%       subVecs2 = extract_subvectors({1:5, 2:6}, 'EndPoints', {[1, 3], [2, 4]})
%       subVecs3 = extract_subvectors(1:5, 'EndPoints', {[1, 3], [2, 4]})
%       subVecs4 = extract_subvectors({1:5, 2:6}, 'Windows', [2.5, 6.5])
%       subVecs5 = extract_subvectors(magic(3), 'EndPoints', {[1,Inf], [2,Inf], [3,Inf]})
%       subVecs6 = extract_subvectors({{1:5, 2:6}, {1:2}}, 'TreatCellNumAsArray', true)
%       extract_subvectors(3:3:18, 'Indices', [2.5, 5.5])
%       extract_subvectors(3:3:18, 'Indices', {3.5; [2.5, 5.5]})
%
% Outputs:
%       subVecs     - subvectors extracted
%                   specified as a numeric array 
%                       or a cell array of numeric vectors
% Arguments:
%       vecs        - vectors to extract
%                   must be an array
%       varargin    - 'Indices': indices for the subvectors to extract 
%                       Note: if provided, would override 'EndPoints'
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == create_indices(endPoints)
%                   - 'EndPoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], vecs)
%                   - 'IndexStart': first index to extract
%                   must be empty or a numeric vector
%                   default == ones(nVectors, 1)
%                   - 'IndexEnd': last index to extract
%                   must be empty or a numeric vector
%                   default == numel(vector) * ones(nVectors, 1)
%                   - 'MaxNum': maximum number of indices to extract
%                   must be a positive integer scalar or Inf
%                   default == Inf
%                   - 'Windows': value windows to extract 
%                       Note: this assumes that the values are nondecreasing
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == []
%                   - 'AlignMethod': method for truncation or padding
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'leftAdjust'  - align to the left and truncate
%                       'rightAdjust' - align to the right and truncate
%                       'leftAdjustPad'  - align to the left and pad
%                       'rightAdjustPad' - align to the right and pad
%                       'none'        - no alignment/truncation
%                   default == 'none'
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'ForceCellOutput': whether to force output as 
%                                           a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/apply_iteratively.m
%       cd/argfun.m
%       cd/create_default_endpoints.m
%       cd/count_samples.m
%       cd/create_empty_match.m
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%       cd/isnumericvector.m
%       cd/ispositiveintegervector.m
%       cd/match_format_vector_sets.m
%       cd/unique_custom.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_relative_event_times.m
%       cd/compute_rms_error.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sweep_errors.m
%       cd/extract_columns.m
%       cd/extract_common_directory.m
%       cd/find_passive_params.m
%       cd/filter_and_extract_pulse_response.m
%       cd/force_matrix.m
%       cd/iscellvector.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_individual_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/parse_multiunit.m
%       cd/parse_pulse_response.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_histogram.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/select_similar_values.m

% File History:
% 2018-10-28 Created by Adam Lu
% 2018-10-29 Now returns empty if input is empty
% 2018-12-07 Now allows vecs to be a numeric array
% 2018-12-07 Now allows endPoints to be empty
% 2018-12-17 Now allows subvectors to be extracted from arbitrary indices
% 2018-12-17 Now allows subvectors to be extracted from an align mode
% 2018-12-24 Added 'IndexStart', 'IndexEnd' as arguments
% 2018-12-27 Added the align methods 'leftAdjustPad', 'rightAdjustPad'
% 2019-01-03 Added usage of force_column_vector.m
% 2019-01-03 Moved code to create_empty_match.m
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 2019-01-04 Fixed bugs for cellstrs
% 2019-04-26 Added padarray_custom
% 2019-04-26 Fixed padding when there are empty vectors
% 2019-08-13 Added 'MaxNum' as an optional parameter
% 2019-08-21 Now allows indices to be non-integers (uses interpolation)
% 2019-09-07 Now matches vectors with the number of indices vectors
% 2019-09-07 Added 'ForceCellOutput' as an optional argument
% 2019-10-03 Added 'TreatCellNumAsArray' as an optional argument
% TODO: check if all endpoints have 2 elements
% 

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad', 'none'};

%% Default values for optional arguments
indicesDefault = [];            % set later
endPointsDefault = [];          % set later
indexStartDefault = [];         % set later
indexEndDefault = [];           % set later
maxNumDefault = Inf;            % no limit by default
windowsDefault = [];            % extract entire trace(s) by default
alignMethodDefault  = 'none';   % no alignment/truncation by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
forceCellOutputDefault = false;     % don't force as cell array by default

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
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IndexStart', indexStartDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexStart must be either empty or a numeric vector!'));
addParameter(iP, 'IndexEnd', indexEndDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexEnd must be either empty or a numeric vector!'));
addParameter(iP, 'MaxNum', maxNumDefault, ...
    @(x) assert(isinf(x) || isaninteger(x), ...
                'MaxNum must be either Inf or a positive integer scalar!'));
addParameter(iP, 'Windows', windowsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Windows must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs, varargin{:});
indices = iP.Results.Indices;
endPoints = iP.Results.EndPoints;
indexStart = iP.Results.IndexStart;
indexEnd = iP.Results.IndexEnd;
maxNum = iP.Results.MaxNum;
windows = iP.Results.Windows;
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
forceCellOutput = iP.Results.ForceCellOutput;

% If indices is provided and endPoints or windows is also provided, 
%   display warning
if ~isempty(indices) && (~isempty(endPoints) || ~isempty(windows))
    fprintf('Endpoints will be ignored because indices are provided!\n');
end

% If endPoints is provided and windows is also provided, display warning
if ~isempty(endPoints) && ~isempty(windows)
    fprintf('Windows will be ignored because end points are provided!\n');
end

% TODO: check if all endpoints have 2 elements

%% Preparation
% If vecs is empty, or if all options are empty, return vecs
if isempty(vecs) || ...
        isempty(indices) && isempty(endPoints) && isempty(windows) && ...
            isempty(indexStart) && isempty(indexEnd) && ...
            strcmpi('AlignMethod', 'none')
    subVecs = vecs;
    return
end

% Find end points if not provided
if isempty(endPoints)
    if ~isempty(windows)
        % Check if the input has the right type
        if ~isnumeric(vecs) && ~iscellnumeric(vecs)
            error('Input expected to be numeric when windows are provided!\n');
        end

        % TODO: Also check if vecs are nondecreasing. 
        %   If not, ask for endPoints or indices to be passed

        % Extract the start and end indices of the vectors for fitting
        endPoints = find_window_endpoints(windows, vecs);
    else
        % Count the number of elements in each vector
        nSamples = count_samples(vecs, 'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray);

        % Construct end points
        endPoints = create_default_endpoints(nSamples);
    end
end

% Create indices if not provided
if isempty(indices)
    indices = create_indices(endPoints, 'Vectors', vecs, ...
                            'IndexStart', indexStart, 'IndexEnd', indexEnd, ...
                            'MaxNum', maxNum, ...
                            'TreatCellAsArray', treatCellAsArray, ...
                            'TreatCellNumAsArray', treatCellNumAsArray, ...
                            'TreatCellStrAsArray', treatCellStrAsArray);
end

%% Do the job
if isempty(indices)
    % If no indices are found, return empty match
    if iscell(vecs)
        subVecs = create_empty_match(vecs);
    else
        subVecs = [];
    end
else
    % If there is a alignment method used, apply it to indices
    indices = align_subvectors(indices, alignMethod);

    % If one of indices and vecs is a cell array, 
    %   match the formats of indices and vecs so that cellfun can be used
    if iscell(indices) || iscell(vecs)
        [indices, vecs] = ...
            match_format_vector_sets(indices, vecs, ...
                                'ForceCellOutputs', false, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray);
    end

    % Extract subvectors
    if iscellnumeric(vecs) && ~treatCellNumAsArray
        subVecs = cellfun(@(x, y) extract_subvectors_helper(x, y), ...
                            vecs, indices, 'UniformOutput', false);
    elseif iscell(vecs) && ~treatCellAsArray && ...
                ~(iscellnumeric(vecs) && treatCellNumAsArray) && ...
                ~(iscellstr(vecs) && treatCellStrAsArray)
        subVecs = cellfun(@(x, y) extract_subvectors(x, 'Indices', y, ...
                                'TreatCellAsArray', treatCellAsArray, ...
                                'TreatCellNumAsArray', treatCellNumAsArray, ...
                                'TreatCellStrAsArray', treatCellStrAsArray), ...
                            vecs, indices, 'UniformOutput', false);
    else
        subVecs = extract_subvectors_helper(vecs, indices);
    end
end

% Force as cell array of column vectors if requested
if forceCellOutput
    subVecs = force_column_cell(subVecs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subVec = extract_subvectors_helper (vec, indices)
%% Extract a subvector from vector(s) if not empty

% If indices all NaNs, it is for padding, so return it as the subvector
%   If the time window is out of range, return an empty vector
if all(all(isnan(indices)))
    subVec = indices;
    return
elseif isempty(indices) || isempty(vec)
    subVec = [];
    return
end

% Make sure vectors and indices are in columns
[vec, indices] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true, ...
                                        'TreatCharAsScalar', true, ...
                                        'TreatCellAsArray', true, ...
                                        'TreatCellNumAsArray', true, ...
                                        'TreatCellStrAsArray', true), ...
            vec, indices);

% Count the number of different indices vectors
nIndices = size(indices, 2);

% Count the desired number of rows
nRows = size(indices, 1);

% Count the desired number of columns
nColumns = size(vec, 2);

% If there is only one vector, repmat it to match the number of indices
%   vectors
if nColumns == 1 && nIndices > 1 
    vec = repmat(vec, 1, nIndices);
    nColumns = nIndices;
end

% Make sure nIndices is either 1 or nColumns
if (nIndices ~= 1 && nIndices ~= nColumns)
    error('nIndices not correct!');
end

% Initialize subVec with NaNs
subVec = create_empty_match(vec, 'NRows', nRows, 'NColumns', nColumns);

% Find the parts of indices without NaNs
% TODO: Make this a function find_nonempty_indices.m
%       and return without NaTs if necessary
if nIndices == 1
    withoutNaNs = find(~isnan(indices));
else
    withoutNaNs = arrayfun(@(x) find(~isnan(indices(:, x))), ...
                            transpose(1:nIndices), 'UniformOutput', false);
end

% Extract the subvectors
if nIndices == 0
    % do nothing
elseif nIndices == 1
    % Replace the parts that are not NaNs
    subVec(withoutNaNs, :) = extract_one_subvector(vec, indices(withoutNaNs));
else
    for iCol = 1:nColumns
        % Get the indices for this column
        indicesThis = indices(:, iCol);

        % Get the parts without NaNs for the indices of this column
        withoutNaNsThis = withoutNaNs{iCol};

        % Extract this column
        vecThis = vec(:, iCol);

        % Replace the parts of this column that are not NaNs
        subVec(withoutNaNsThis, iCol) = ...
            extract_one_subvector(vecThis, indicesThis(withoutNaNsThis));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subVec = extract_one_subvector (vec, indices)
% Extract just one subvector

if ispositiveintegervector(indices)
    % Just get the subvector
    subVec = vec(indices, :);
else
    % Set up indices
    allIndices = 1:size(vec, 1);

    % Interpolate
    subVec = interp1(allIndices, vec, indices);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = align_subvectors(indices, alignMethod)
%% TODO: Description
% TODO: Pull out to its own function?

switch alignMethod
case {'leftAdjust', 'rightAdjust', 'leftAdjustPad', 'rightAdjustPad'}
    % Force as column vectors
    indices = force_column_vector(indices, 'IgnoreNonVectors', true);
    
    % Count the number of elements in each vector
    nSamples = count_samples(indices, 'ForceColumnOutput', true, ...
                                'CountMethod', 'nrows');

    % Get unique nSamples
    uniqueNSamples = apply_iteratively(@unique_custom, nSamples, ...
                                        'IgnoreNaN', true);

    % Get the number of unique nSamples
    nUniqueNSamples = numel(uniqueNSamples);

    % If nSamples are all equal, do nothing
    if nUniqueNSamples == 1
        return
    end
case 'none'
    % Do nothing
otherwise
    error_unrecognized(get_var_name(alignMethod), alignMethod, mfilename);
end

switch alignMethod
case {'leftAdjust', 'rightAdjust'}
    % Get the minimum number of samples
    nSamplesTarget = min(uniqueNSamples);

    % Create new endpoints
    switch alignMethod
    case 'leftAdjust'
        % Create new endpoints by repetition
        newEndPoints = match_format_vector_sets([1; nSamplesTarget], indices);
    case 'rightAdjust'
        % Get the new starting positions 
        newStarts = nSamples - nSamplesTarget + 1;

        % Create new endpoints
        newEndPoints = transpose([newStarts, nSamples]);
    end

    % Extract indices from new endpoints
    indices = extract_subvectors(indices, 'EndPoints', newEndPoints, ...
                                'AlignMethod', 'none');
case {'leftAdjustPad', 'rightAdjustPad'}
    % Get the maximum number of samples
    nSamplesTarget = max(uniqueNSamples);

    % Count the number of rows to pad for each vector
    nRowsToPad = nSamplesTarget - nSamples;

    % Create new endpoints
    switch alignMethod
    case 'leftAdjustPad'
        padDirection = 'post';
    case 'rightAdjustPad'
        padDirection = 'pre';
    end

    % Add NaNs to indices
    if iscell(indices)
        indices = cellfun(@(x, y) padarray_custom(x, y, NaN, padDirection), ...
                            indices, num2cell(nRowsToPad), ...
                            'UniformOutput', false);
    else
        indices = padarray_custom(indices, nRowsToPad, NaN, padDirection);
    end
case 'none'
    % Do nothing
    nSamplesTarget = NaN;
otherwise
    error_unrecognized(get_var_name(alignMethod), alignMethod, mfilename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function padded = padarray_custom(array, padSize, padValue, padDirection)
%% Pads array with padValue and creates array with padValue if empty

if isempty(array)
    padded = padValue * ones(padSize, 1);
else
    padded = padarray(array, padSize, padValue, padDirection);
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

% Force row vectors as column vectors
[endPoints, windows] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
            endPoints, windows);

%       cd/force_column_vector.m


%                   must be a numeric array or a cell array of numeric arrays

@(x) assert(isnumeric(x) || iscellnumeric(x), ...
            ['vecs must be either a numeric array ', ...
                'or a cell array of numeric arrays!']));

%   Note: If first argument is empty, 
%           the first and last indices will be returned

nRows = length(indices);

% If vecs is empty, or if all options are empty, return vecs

if iscellnumericvector(vecs)
[indices, vecs] = ...
        match_format_vector_sets(indices, vecs, 'ForceCellOutputs', true);

uniqueNSamples = unique(nSamples);

subVec(withoutNaNs, :) = vec(indices(withoutNaNs), :);
subVec(withoutNaNsThis, iCol) = vecThis(indicesThis(withoutNaNsThis));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
