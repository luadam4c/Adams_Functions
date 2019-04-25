function indices = create_indices (varargin)
%% Creates indices from endpoints (starting and ending indices)
% Usage: indices = create_indices (varargin)
% Explanation:
%       TODO
% Example(s):
%       create_indices([2, 5])
%       create_indices([5, 1])
%       create_indices([2; 5])
%       create_indices([[2; 2], [3; 3]])
%       create_indices({[2, 5], [2, 5]})
%       create_indices({[2; 5]; [2; 5]})
%       create_indices('IndexEnd', 5)
%       create_indices('IndexEnd', [2, 3])
%       create_indices([1, 50], 'MaxNum', 5)
%       create_indices([5; 1], 'MaxNum', 2)
% Outputs:
%       indices     - indices for each pair of idxStart and idxEnd
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
% Arguments:
%       endPoints   - (opt) the starting and ending indices
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%       varargin    - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ForceRowOutput': whether to force as row vector instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Vectors': vectors to create indices from
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%                   - 'IndexStart': first index
%                   must be empty or a numeric vector
%                   default == ones(nVectors, 1)
%                   - 'IndexEnd': last index
%                   must be empty or a numeric vector
%                   default == numel(vector) * ones(nVectors, 1)
%                   - 'MaxNum': maximum number of indices
%                   must be a positive integer scalar or Inf
%                   default == Inf
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/argfun.m
%       cd/compute_maximum_trace.m
%       cd/compute_minimum_trace.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/isemptycell.m
%       cd/isnumericvector.m
%       cd/match_and_combine_vectors.m
%       cd/match_format_vector_sets.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/plot_measures.m
%       cd/compute_peak_decay.m
%       cd/extract_columns.m
%       cd/extract_subvectors.m
%       cd/fit_2exp.m
%       cd/compute_combined_data.m
%       cd/parse_pulse_response.m
%       cd/plot_traces.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2018-12-18 Added 'ForceCellOutput' as an optional argument
% 2018-12-22 Now replaces out of range values with the actual endpoints
% 2018-12-24 Added 'IndexStart' and 'IndexEnd' as optional arguments
% 2019-01-13 Added 'ForceRowOutput' as an optional argument
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 2019-01-23 Now avoids putting indices together as a matrix if there is
%               only one index per vector
% 2019-02-24 Added 'MaxNum' as an optional parameter
% 2019-03-25 Now prevents negative indices by default
% 2019-04-24 Now allows indices to decrement
% TODO: Use argument 'ForcePositive' as false where necessary

%% Hard-coded parameters

%% Default values for optional arguments
endPointsDefault = [];          % no endpoint by default
forcePositiveDefault = true;    % force indices to be positive by default
forceCellOutputDefault = false; % don't force output as a cell array by default
forceRowOutputDefault = false;  % return column vectors by default
vectorsDefault = [];
indexEndDefault = [];           % set later
indexStartDefault = [];         % set later
maxNumDefault = Inf;            % no limit by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'endPoints', endPointsDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForcePositive', forcePositiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ForceRowOutput', forceRowOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Vectors', vectorsDefault, ...
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['Vectors must be either a numeric array', ...
                    'or a cell array!']));
addParameter(iP, 'IndexStart', indexStartDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexStart must be either empty or a numeric vector!'));
addParameter(iP, 'IndexEnd', indexEndDefault, ...
    @(x) assert(isnumericvector(x), ...
                'IndexEnd must be either empty or a numeric vector!'));
addParameter(iP, 'MaxNum', maxNumDefault, ...
    @(x) assert(isinf(x) || isaninteger(x), ...
                'MaxNum must be either Inf or a positive integer scalar!'));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
endPoints = iP.Results.endPoints;
forcePositive = iP.Results.ForcePositive;
forceCellOutput = iP.Results.ForceCellOutput;
forceRowOutput = iP.Results.ForceRowOutput;
vectors = iP.Results.Vectors;
indexStartUser = iP.Results.IndexStart;
indexEndUser = iP.Results.IndexEnd;
maxNum = iP.Results.MaxNum;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;

% TODO: warn if EndPoints provided and IndexStart and IndexEnd also provided

%% Preparation
% If nothing provided, return empty
if isempty(endPoints) && isempty(vectors) && ...
        isempty(indexStartUser) && isempty(indexEndUser)
    indices = [];
    return
end

% Decide whether to output row vectors
if forceRowOutput || isvector(endPoints) && ~iscolumn(endPoints) || ...
        isvector(indexStartUser) && ~iscolumn(indexStartUser) || ...
        isvector(indexEndUser) && ~iscolumn(indexEndUser)
    rowInstead = true;
else
    rowInstead = false;
end

% Make sure endPoints, indexStartUser, indexEndUser are in columns
[endPoints, indexStartUser, indexEndUser] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
            endPoints, indexStartUser, indexEndUser);

%% Do the job
% Extract the starting and ending indices from provided end points
[idxStartFromEndPoints, idxEndFromEndPoints] = ...
    argfun(@(x) extract_elements(endPoints, x), 'first', 'last');

% Replace with indexStartUser if provided
if ~isempty(indexStartUser)
    idxStart = indexStartUser;
else
    idxStart = idxStartFromEndPoints;
end

% Replace with indexEndUser if provided
if ~isempty(indexEndUser)
    idxEnd = indexEndUser;
else
    idxEnd = idxEndFromEndPoints;
end

% If one of idxStart and idxEnd is empty, create the other
if isempty(idxStart) && ~isempty(idxEnd)
    % Start the indices from 1
    idxStart = ones(size(idxEnd));
elseif ~isempty(idxStart) && isempty(idxEnd)
    % Use the length of the vectors
    idxEnd = count_samples(vectors);
end

% Make sure indices are columns
[idxStart, idxEnd] = argfun(@force_column_vector, idxStart, idxEnd);

% Match the formats of idxStart and idxEnd
[idxStart, idxEnd] = ...
    match_format_vector_sets(idxStart, idxEnd, 'MatchVectors', true);

% If original vectors are provided, 
%   fix indices if they are out of range
if ~(isempty(vectors) || iscell(vectors) && all(all(isemptycell(vectors))))
    % Count the number of samples in each vector
    nSamples = count_samples(vectors, 'TreatCellAsArray', treatCellAsArray, ...
                                    'TreatCellStrAsArray', treatCellStrAsArray);

    % Match the vector counts
    [idxStart, idxEnd] = ...
        argfun(@(x) match_format_vector_sets(x, nSamples, 'MatchVectors', true), ...
                idxStart, idxEnd);

    % Make sure endpoint indices are in range
    candidatesStart = match_and_combine_vectors(idxStart, 1);
    candidatesEnd = match_and_combine_vectors(idxEnd, nSamples);
    idxStart = compute_maximum_trace(candidatesStart, 'TreatRowAsMatrix', true);
    idxEnd = compute_minimum_trace(candidatesEnd, 'TreatRowAsMatrix', true);
end

% Create the indices
if iscell(idxStart) && iscell(idxEnd)
    indices = cellfun(@(x, y) create_indices_helper(x, y, maxNum, forcePositive), ...
                        idxStart, idxEnd, 'UniformOutput', false);
elseif isnumeric(idxStart) && isnumeric(idxEnd)
    indices = create_indices_helper(idxStart, idxEnd, maxNum, forcePositive);
else
    error('idxStart and idxEnd don''t match!!');
end

% Force as cell array output if requested
if forceCellOutput && ~iscell(indices)
    indices = force_column_cell(indices, 'RowInstead', rowInstead);
end

% Make the output consistent with the input
indices = force_column_vector(indices, 'RowInstead', rowInstead, ...
                                'IgnoreNonVectors', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = create_indices_helper (idxStart, idxEnd, maxNum, forcePositive)

% Match idxStart and idxEnd
[idxStart, idxEnd] = match_format_vectors(idxStart, idxEnd);

% Construct vectors of indices
if numel(idxStart) == 1 && numel(idxEnd) == 1
    % There is just one vector
    indices = create_one_indices(idxStart, idxEnd, maxNum, forcePositive);
else
    % There are multiple vectors
    indices = arrayfun(@(x, y) create_one_indices(x, y, maxNum, forcePositive), ...
                        idxStart, idxEnd, 'UniformOutput', false);

    % Count the number of samples in each indices vector
    nSamples = count_samples(indices);

    % Extract unique number of samples
    uniqueNSamples = unique(nSamples);
    
    % If the number of samples are all the same across all indices vectors
    %   unless all nSamples is one
    %   Note: 'AlignMethod' must be 'none' to prevent infinite loop
    if numel(uniqueNSamples) == 1 && uniqueNSamples ~= 1
        indices = force_matrix(indices, 'AlignMethod', 'none');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = create_one_indices (idxStart, idxEnd, maxNum, forcePositive)
%% Creates one set of indices

% Force the starting index to be positive if requested
if forcePositive && idxStart < 1
    idxStart = 1;
end

% Get the sign of the index increment
sgnIncr = sign(idxEnd - idxStart);

% Decide on the index increment
if ~isinf(maxNum)
    % Count the number of indices
    nIndices = abs(idxEnd - idxStart) + 1;

    % Decide on the index increment
    if nIndices > maxNum
        % Update index increment
        idxIncr = sgnIncr * ceil(nIndices / maxNum);

        % Update new number of indices
        nIndicesNew = ceil(nIndices / abs(idxIncr));

        % Update index start
        idxStart = idxEnd - idxIncr * (nIndicesNew - 1);
    else
        idxIncr = sgnIncr * 1;
    end
else
    idxIncr = sgnIncr * 1;
end

% Create the indices vector
indices = transpose(idxStart:idxIncr:idxEnd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

nSamples = match_dimensions(nSamples, size(idxStart));

% Create ones
firsts = ones(size(idxStart));

% Make sure endpoint indices are in range
idxStart = max([idxStart, firsts], [], 2);
idxEnd = min([idxEnd, nSamples], [], 2);

% Match the vector counts
[nSamples, idxStart, idxEnd] = ...
    match_format_vectors(nSamples, idxStart, idxEnd);

indices = {indices};

parse(iP, endPoints, varargin{:});

indices = arrayfun(@(x, y) transpose(x:y), idxStart, idxEnd, ...
                'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
