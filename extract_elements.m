function [elements, idxElement] = extract_elements (vecs, extractMode, varargin)
%% Extracts elements from vectors using a certain mode ('first', 'last', 'min', 'max', 'all')
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
%       extract_elements({[3; 5], []; [], [1, 2]}, 'specific', 'Index', 1)
%       extract_elements([2, 3], 'first')
%       [elements, idxElement] = extract_elements({[1 2 3], [10 20], [100 200 300]}, 'all')
%       % ans{1} = [1 10 100]; ans{2} = [2 20 200]; ans{3} = [3 300]
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
%                       'first'     - first element of each vector
%                       'last'      - last element of each vector
%                       'center'    - center element of each vector
%                       'min'       - minimum-valued element of each vector
%                       'max'       - maximum-valued element of each vector
%                       'maxabs'    - maximum-absolute-valued element
%                       'firstdiff' - first difference of each vector
%                       'specific'  - at a a specific index
%                       'all'       - all elements by column
%       varargin    - 'Index': index of the element from each vector
%                   must be a positive numeric vector
%                   default == []
%                   - 'ReturnNan': Return NaN instead of empty
%                                       if nothing to extract
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/array_fun.m
%       cd/count_vectors.m
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/isaninteger.m
%       cd/iscellnumericvector.m
%       cd/isnumericvector.m
%       cd/match_dimensions.m
%       cd/match_format_vector_sets.m
%       cd/remove_empty.m
%
% Used by:
%       cd/adjust_edges.m
%       cd/align_subplots.m
%       cd/compute_peak_decay.m
%       cd/compute_peak_halfwidth.m
%       cd/compute_sampsizepwr.m
%       cd/compute_time_window.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/extract_columns.m
%       cd/find_closest.m
%       cd/is_overlapping.m
%       cd/m3ha_find_decision_point.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_simulate_population.ms
%       cd/parse_ipsc.m
%       cd/parse_lts.m
%       cd/parse_multiunit.m
%       cd/parse_pulse_response.m
%       cd/parse_repetitive_pulses.m
%       cd/parse_stim.m
%       cd/plot_autocorrelogram.m
%       cd/plot_repetitive_protocols.m
%       cd/resize_subplots_for_labels.m
%       cd/select_similar_values.m
%       cd/update_figure_for_corel.m
%       \Shared\Code\vIRt\virt_moore.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-17 Now returns idxElement as well
% 2019-01-03 Now accepts cell arrays of non-vector arrays
% 2019-02-20 Fixed bug when using match_dimensions
% 2019-03-14 Now returns NaN if the vector is not long enough
% 2019-09-08 Added 'center' as a possible extract mode
% 2019-10-04 Fixed bug
% 2019-12-01 Now allows vecs to be a matrix cell array
% 2020-02-26 Added 'ReturnNan' as a an optional arguement
% 2020-04-22 Added 'maxabs' as a possible extract mode
% TODO: Add 'TreatCellAsArray' as a parameter
% TODO: Add 'MaxNum' as an optional argument with default Inf
% TODO: Add 'Indices', 'Endpoints' and 'Windows' as optional arguments
%           and use extract_subvectors.m
% 2025-08-23 Added 'all' mode to extract elements by column.

%% Hard-coded parameters
validExtractModes = {'first', 'last', 'center', 'min', 'max', 'maxabs', ...
                        'firstdiff', 'specific', 'all'};

% TODO: Add as optional argument
treatCellAsArray = false;   % TODO: Not working yet. UniformOutput might need to be set according to this

%% Default values for optional arguments
indexDefault = [];
returnNanDefault = true;    % whether to return NaN instead of empty
                            %   if nothing to extract by default

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
addRequired(iP, 'vecs');
addRequired(iP, 'extractMode', ...
    @(x) any(validatestring(x, validExtractModes)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Index', indexDefault, ...
    @(x) assert(isnumericvector(x), 'Index must be a numeric vector!'));
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs, extractMode, varargin{:});
index = iP.Results.Index;
returnNan = iP.Results.ReturnNan;

% Validate extractMode
extractMode = validatestring(extractMode, validExtractModes);

%% Do the job
switch extractMode
case 'all'
    % This mode only works for cell arrays of vectors
    if ~iscell(vecs)
        error('The ''all'' extractMode only supports cell array inputs.');
    end

    % Determine which algorithm to use
    % Check if the script force_matrix.m exists in path
    if exist('force_matrix.m', 'file') == 2
        % Use custom method

        % First force all vectors as columns in a matrix padded by NaNs or
        %  empty cells
        matrix = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

        % Then transpose the matrix
        matrixTransposed = matrix';

        % Then extract each column as a cell array of vectors
        elementsWithEmpty = force_column_cell(matrixTransposed);

        % Remove empty
        elements = cellfun(@remove_empty, elementsWithEmpty, 'UniformOutput', false);

        % Set indices
        idxElement = (1:numel(elements))';
    else
        % Use Gemini-generated method

        % Find the length of the longest vector in the input cell array
        vectorLengths = cellfun(@numel, vecs);
        if isempty(vectorLengths)
            maxLength = 0;
        else
            maxLength = max(vectorLengths);
        end
    
        % Pre-allocate the output cell arrays for performance
        elements = cell(1, maxLength);
        idxElement = cell(1, maxLength);
    
        % Iterate through each "column" index, from 1 to the max length
        for i = 1:maxLength
            % For each column, build a new vector of elements and their indices
            currentElements = [];
            currentIndices = [];
    
            % Iterate through each of the original vectors
            for j = 1:numel(vecs)
                % Check if the current vector has an element at index 'i'
                if numel(vecs{j}) >= i
                    % If it does, append the element and its index
                    currentElements(end + 1) = vecs{j}(i);
                    currentIndices(end + 1) = i;
                end
            end
    
            % Assign the newly created vectors to the output cell arrays
            elements{i} = currentElements;
            idxElement{i} = currentIndices;
        end
    end

case {'first', 'last', 'center', 'min', 'max', 'maxabs', 'firstdiff'}
    % Extract from a position
    if iscellnumericvector(vecs) && ~treatCellAsArray
        [elements, idxElement] = ...
            array_fun(@(x) extract_by_position(x, extractMode, returnNan), ...
                        vecs);
    elseif iscell(vecs) && ~treatCellAsArray
        try
            [elements, idxElement] = ...
                array_fun(@(x) extract_elements(x, extractMode, ...
                                            'ReturnNan', returnNan), vecs);
        catch
            [elements, idxElement] = ...
                array_fun(@(x) extract_elements(x, extractMode, ...
                                            'ReturnNan', returnNan), vecs, ...
                        'UniformOutput', false);
        end
    else
        [elements, idxElement] = ...
            array_fun(@(x) extract_by_position(vecs(:, x), extractMode, ...
                                                returnNan), ...
                    transpose(1:size(vecs, 2)));
    end
case 'specific'
    % Check if index is provided
    if isempty(index)
        error(['The index for each vector must be provided ', ...
                'under the ''specific'' extract mode!!']);
    end

    % Extract from a specific index
    if iscell(vecs) && ~treatCellAsArray
        % Force as a cell array
        if isnumeric(index)
            index = num2cell(index);
        end

        % Match the vector counts
        [vecs, index] = match_format_vector_sets(vecs, index);

        % Extract by index on each vector
        [elements, idxElement] = ...
            array_fun(@(x, y) extract_by_index(x, y, returnNan), vecs, index);
    else
        % Count the number of vectors
        nVectors = count_vectors(vecs);

        % Force as a column vector
        vecs = force_column_vector(vecs, 'IgnoreNonVectors', true);
        
        % Match the number of indices
        index = match_dimensions(index, [nVectors, 1]);

        % Extract by index on each comlumn
        [elements, idxElement] = ...
            array_fun(@(x, y) extract_by_index(vecs(:, x), y, returnNan), ...
                    transpose(1:nVectors), index);
    end
otherwise
    error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element, idxElement] = extract_by_position (x, extractMode, returnNan)

if numel(x) < 1
    if returnNan
        element = NaN;
        idxElement = NaN;
    else
        element = [];
        idxElement = [];
    end
    return
end

switch extractMode
    case 'first'
        idxElement = 1;
        element = get_element(x, idxElement);
    case 'last'
        idxElement = numel(x);
        element = get_element(x, idxElement);
    case 'center'
        nElements = numel(x);
        idxElement = (nElements + 1) / 2;
        if ~isaninteger(idxElement)
            idxElement = idxElement - 0.5;
        end
        element = get_element(x, idxElement);
    case 'middle'
        nElements = numel(x);
        idxElement = (nElements + 1) / 2;
        if isaninteger(idxElement)
            element = get_element(x, idxElement);
        else
            element = nanmean([get_element(x, idxElement - 0.5); ...
                                get_element(x, idxElement + 0.5)]);
        end
    case 'min'
        [element, idxElement] = min(x);
    case 'max'
        [element, idxElement] = max(x);
    case 'maxabs'
        [element, idxElement] = max(abs(x));
    case 'firstdiff'
        element = x(2) - x(1);
        idxElement = NaN;
    otherwise
        error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element, idxElement] = extract_by_index (x, index, returnNan)

% Decide on the actual index
if numel(x) == 0
    if returnNan
        idxElement = NaN;
    else
        idxElement = [];
    end
elseif isnan(index)
    idxElement = NaN;
elseif index == Inf
    idxElement = numel(x);
elseif index == -Inf
    idxElement = 1;
elseif index >= 1 && index <= numel(x)
    idxElement = index;
else
    fprintf('Warning: The index %g is out of bounds!\n', index);
    if returnNan
        idxElement = NaN;
    else
        idxElement = [];
    end
end

% Get the element
element = get_element(x, idxElement, returnNan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function element = get_element (vector, index, returnNan)
%% Retrieves an element based on vector type

% Return NaN if index is NaN
if isempty(index)
    if returnNan
        element = NaN;
    else
        element = [];
    end
    return
elseif isnan(index)
    element = NaN;
    return
end

% Get element based on vector type
if iscell(vector)
    element = vector{index};
else
    element = vector(index);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) assert(isnumeric(x) || iscellnumeric(x), ...
            ['vecs must be either a numeric array', ...
                'or a cell array of numeric arrays!']));

addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of vectors!']));

        [elements, idxElement] = ...
            array_fun(@(x) extract_elements(x, extractMode), vecs, ...
                    'UniformOutput', false);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
