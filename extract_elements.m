function [elements, idxElement] = extract_elements (vecs, extractMode, varargin)
%% Extracts elements from vectors using a certain mode ('first', 'last', 'min', 'max')
% Usage: [elements, idxElement] = extract_elements (vecs, extractMode, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'first')
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'last')
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'min')
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'max')
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'firstdiff')
%       extract_elements({[3; 5; 4], [], [1, 2, -1]}, 'center')
%
% Outputs:
%       elements    - element(s) from each vector extracted
%                   specified as a numeric vector 
%                       or a cell array of column numeric vectors
%       idxElement  - indices(s) of elements extracted
%                   specified as a numeric vector 
%                       or a cell array of column numeric vectors
%
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array
%       extractMode - mode of extraction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - first element of each vector
%                       'last'  - last element of each vector
%                       'center'- center element of each vector
%                       'min'   - minimum-valued element of each vector
%                       'max'   - maximum-valued element of each vector
%                       'firstdiff' - first difference of each vector
%                       'specific'  - at a a specific index
%       varargin    - 'Index': index of the element from each vector
%                   must be a positive numeric vector
%                   default == []
%
% Requires:
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/isaninteger.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_dimensions.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_time_window.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/parse_lts.m
%       cd/parse_multiunit.m
%       cd/parse_pulse_response.m
%       cd/parse_stim.m
%       cd/plot_protocols.m
%       cd/select_similar_values.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-17 Now returns idxElement as well
% 2019-01-03 Now accepts cell arrays of non-vector arrays
% 2019-02-20 Fixed bug when using match_dimensions
% 2019-03-14 Now returns NaN if the vector is not long enough
% 2019-09-08 Added 'center' as a possible extract mode
% TODO: Add 'MaxNum' as an optional argument with default Inf
% TODO: Add 'Indices', 'Endpoints' and 'Windows' as optional arguments
%           and use extract_subvectors.m
% 

%% Hard-coded parameters
validExtractModes = {'first', 'last', 'center', 'min', 'max', ...
                        'firstdiff', 'specific'};

%% Default values for optional arguments
indexDefault = [];

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
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of vectors!']));
addRequired(iP, 'extractMode', ...
    @(x) any(validatestring(x, validExtractModes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Index', indexDefault, ...
    @(x) assert(isnumericvector(x), 'Index must be a numeric vector!'));

% Read from the Input Parser
parse(iP, vecs, extractMode, varargin{:});
index = iP.Results.Index;

% Validate extractMode
extractMode = validatestring(extractMode, validExtractModes);

%% Do the job
switch extractMode
case {'first', 'last', 'center', 'min', 'max', 'firstdiff'}
    % Extract from a position
    if iscellnumericvector(vecs)
        [elements, idxElement] = ...
            cellfun(@(x) extract_by_position(x, extractMode), vecs);
    elseif iscell(vecs)
        [elements, idxElement] = ...
            cellfun(@(x) extract_elements(x, extractMode), vecs, ...
                    'UniformOutput', false);
    else
        [elements, idxElement] = ...
            arrayfun(@(x) extract_by_position(vecs(:, x), extractMode), ...
                    transpose(1:size(vecs, 2)));
    end
case 'specific'
    % Check if index is provided
    if isempty(index)
        error(['The index for each vector must be provided ', ...
                'under the ''specific'' extract mode!!']);
    end

    % Extract from a specific index
    if iscell(vecs)
        % Force as a cell array
        if isnumeric(index)
            index = num2cell(index);
        end

        % Match the vector counts
        [vecs, index] = match_format_vector_sets(vecs, index);

        % Extract by index on each vector
        [elements, idxElement] = ...
            cellfun(@(x, y) extract_by_index(x, y), vecs, index);
    else
        % Count the number of vectors
        nVectors = count_vectors(vecs);

        % Match the number of indices
        index = match_dimensions(index, [nVectors, 1]);

        % Extract by index on each comlumn
        [elements, idxElement] = ...
            arrayfun(@(x, y) extract_by_index(vecs(:, x), y), ...
                    transpose(1:nVectors), index);
    end
otherwise
    error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element, idxElement] = extract_by_position (x, extractMode)

if numel(x) < 1
    element = NaN;
    idxElement = NaN;
    return
end

switch extractMode
    case 'first'
        element = x(1);
        idxElement = 1;
    case 'last'
        element = x(end);
        idxElement = numel(x);
    case 'center'
        nElements = numel(x);
        idxElement = (nElements + 1) / 2;
        if isaninteger(idxElement)
            element = x(idxElement);
        else
            element = nanmean([x(idxElement - 0.5), x(idxElement + 0.5)]);
        end
    case 'min'
        [element, idxElement] = min(x);
    case 'max'
        [element, idxElement] = max(x);
    case 'firstdiff'
        element = x(2) - x(1);
        idxElement = NaN;
    otherwise
        error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element, idxElement] = extract_by_index (x, index)

if numel(x) < 1
    element = NaN;
    idxElement = NaN;
    return
end

if isnan(index)
    element = NaN;
    idxElement = NaN;
elseif index == Inf
    element = x(end);
    idxElement = numel(x);
elseif index == -Inf
    element = x(1);
    idxElement = 1;
elseif index >= 1 && index <= numel(x)
    element = x(index);
    idxElement = index;
else
    error('The index %g is out of bounds!', index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) assert(isnumeric(x) || iscellnumeric(x), ...
            ['vecs must be either a numeric array', ...
                'or a cell array of numeric arrays!']));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
