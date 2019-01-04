function indices = create_indices (endPoints, varargin)
%% Creates indices from endpoints (starting and ending indices)
% Usage: indices = create_indices (endPoints, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       indices     - indices for each pair of idxStart and idxEnd
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
% Arguments:
%       endPoints   - the starting and ending indices
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%       varargin    - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Vectors': vectors to create indices from
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array
%                   - 'IndexStart': first index to extract
%                   must be empty or a numeric vector
%                   default == ones(nVectors, 1)
%                   - 'IndexEnd': last index to extract
%                   must be empty or a numeric vector
%                   default == numel(vector) * ones(nVectors, 1)
%
% Requires:
%       cd/argfun.m
%       cd/compute_maximum_trace.m
%       cd/compute_minimum_trace.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/force_column_numeric.m
%       cd/force_matrix.m
%       cd/isnumericvector.m
%       cd/match_and_combine_vectors.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/extract_subvectors.m
%       cd/fit_2exp.m
%       cd/parse_pulse_response.m

% File History:
% 2018-12-17 Created by Adam Lu
% 2018-12-18 Added 'ForceCellOutput' as an optional argument
% 2018-12-22 Now replaces out of range values with the actual endpoints
% 2018-12-24 Added 'IndexStart' and 'IndexEnd' as optional arguments
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputDefault = false; % don't force output as a cell array by default
vectorsDefault = [];
indexEndDefault = [];           % set later
indexStartDefault = [];         % set later

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
addRequired(iP, 'endPoints', ...
    @(x) assert(isempty(x) || isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
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

% Read from the Input Parser
parse(iP, endPoints, varargin{:});
forceCellOutput = iP.Results.ForceCellOutput;
vectors = iP.Results.Vectors;
indexStartUser = iP.Results.IndexStart;
indexEndUser = iP.Results.IndexEnd;

% TODO: warn if EndPoints provided and IndexStart and IndexEnd provided

%% Preparation
% Make sure endPoints are in columns
endPoints = force_column_numeric(endPoints, 'IgnoreNonVectors', true);

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

% If original vectors are provided, 
%   fix indices if they are out of range
if ~isempty(vectors)
    % Make sure indices are columns
    [idxStart, idxEnd] = argfun(@force_column_numeric, idxStart, idxEnd);

    % Count the number of samples in each vector
    nSamples = count_samples(vectors);

    % Match the vector counts
    [idxStart, idxEnd] = ...
        argfun(@(x) match_format_vector_sets(nSamples, x, 'MatchVectors', true), ...
                idxStart, idxEnd);

    % Make sure endpoint indices are in range
    candidatesStart = match_and_combine_vectors(idxStart, 1);
    candidatesEnd = match_and_combine_vectors(idxEnd, nSamples);
    idxStart = compute_maximum_trace(candidatesStart);
    idxEnd = compute_minimum_trace(candidatesEnd);
end

% Create the indices
if iscell(idxStart) && iscell(idxEnd)
    indices = cellfun(@(x, y) create_indices_helper(x, y), ...
                        idxStart, idxEnd, 'UniformOutput', false);
elseif isnumeric(idxStart) && isnumeric(idxEnd)
    indices = create_indices_helper(idxStart, idxEnd);
else
    error('idxStart and idxEnd don''t match!!');
end

% Force as cell array output if requested
if forceCellOutput && ~iscell(indices)
    indices = {indices};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = create_indices_helper (idxStart, idxEnd)

% Match idxStart and idxEnd
[idxStart, idxEnd] = match_format_vectors(idxStart, idxEnd);

% Construct vectors of indices
if numel(idxStart) == 1 && numel(idxEnd) == 1
    % There is just one vector
    indices = transpose(idxStart:idxEnd);
else
    % There are multiple vectors
    indices = arrayfun(@(x, y) transpose(x:y), idxStart, idxEnd, ...
                    'UniformOutput', false);

    % Count the number of samples in each vector
    nSamples = count_samples(indices);

    % If the number of samples are all the same, combine the vectors
    %   Note: 'AlignMethod' must be 'none' to prevent infinite loop
    if numel(unique(nSamples)) == 1
        indices = force_matrix(indices, 'AlignMethod', 'none');
    end
end

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%