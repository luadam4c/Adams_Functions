function varargout = extract_columns (arrays, varargin)
%% Extracts columns from arrays or a cell array of arrays
% Usage: varargout = extract_columns (arrays, colNums (opt), varargin)
% Explanation:
%       TODO
%       By default, cell arrays are not considered as arrays
%
% Example(s):
%       extract_columns(magic(3), [1, 3], 'OutputMode', 'single')
%       [a, b] = extract_columns({magic(3); ones(4)}, [1:3])
%       [a, b] = extract_columns({magic(3); ones(3); zeros(4)}, ...
%                               {[2, 3], [1:3], [1, 3]})
%       c = extract_columns({magic(3); ones(4)}, [1:3], 'OutputMode', 'single')
%       d = extract_columns({{[1, 2]; [2, 1]}, {[4, 5], [3, 2]}}, 1:2, 'TreatCnvAsColumns', true, 'OutputMode', 'single')
%       e = extract_columns(repmat({{[1, 2]; [2, 1]}}, 2, 2), 1:2, 'TreatCnvAsColumns', true, 'OutputMode', 'single')
%       [a, b] = extract_columns({magic(3); ones(4)}, [1:3], 'AsRowVectors', true)
%
% Outputs:
%       varargout   - extracted column #1s, column #2s, etc.
%                       or extracted columns for each array 
%                           if 'OutputMode' is 'single'
%                   specified as a column vector 
%                       or a cell array of column vectors
%                       or a cell array of cell arrays of column vectors
%
% Arguments:    
%       arrays      - arrays to extract columns from
%       colNums     - (opt) column number(s) to extract instead of 1, 2, 3, ...
%                   must be either 'all' or a positive integer vector
%                       or a cell array of positive integer vectors
%                   default == 'all'
%       varargin    - 'OutputMode': mode of outputs
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'multiple' - all columns as separate outputs
%                       'single'   - all columns as one output
%                   default == 'multiple'
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCnvAsColumns': whether to treat a cell array
%                                           of numeric vectors as columns
%                                           of the same array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatRowVecAsOne': whether to treat row vectors
%                                           as a single vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'AsRowVectors': whether to extract as row vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/array_fun.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/create_indices.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/force_column_cell.m
%       cd/force_row_vector.m
%       cd/match_dimensions.m
%       cd/iscellnumeric.m
%       cd/iscellnumericvector.m
%       cd/ispositiveintegerarray.m
%
% Used by:    
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/compute_combined_data.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_network_plot_essential.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_example_jitter.m.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_simulate_population.m

% File History:
% 2018-10-24 Created by Adam Lu
% 2018-10-25 Updated usage of match_dimensions()
% 2018-12-18 Allow the option to treat cell arrays as arrays;
%               added 'TreatCellAsArray' (default == 'false')
% 2019-01-01 Fixed bugs if a cell array has only one numeric array
% 2019-01-12 Fixed bugs for single output and simplified 
% 2019-01-12 Now considers cell arrays of numeric vectors 
%               as a single array when 'TreatCnvAsColumns' is true
% 2019-01-13 Added 'AsRowVectors' as an optional argument
% 2019-01-22 Now uses iscellnumericvector instead of iscellvector
% 2019-12-01 Now allows arrays to be a matrix cell array
% 2020-01-01 Now uses array_fun.m
% TODO: Implement treatRowVecAsOne

%% Hard-coded parameters
validOutputModes = {'multiple', 'single'};

%% Default values for optional arguments
colNumberDefault  = 'all';      % extract all columns by default
outputModeDefault = 'multiple'; % separate outputs by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCnvAsColumnsDefault = false;   % treat cell arrays of numeric vectors as 
                                    % many arrays by default
treatRowVecAsOneDefault = true;     % treat row vectors as one vector by default
asRowVectorsDefault = false;    % extract as column vectors by default

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
addRequired(iP, 'arrays');

% Add optional inputs to the Input Parser
addOptional(iP, 'colNums', colNumberDefault, ...
    @(x) assert((ischar(x) || isstring(x)) && strcmpi(x, 'all') || ...
                ispositiveintegerarray(x) || ...
                iscell(x) && iscellnumeric(x), ...
                ['colNums must be either ''all'', "all" ', ...
                    'or a positive integer array ', ...
                    'or a cell array of positive integer vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputMode', outputModeDefault, ...
    @(x) any(validatestring(x, validOutputModes)));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCnvAsColumns', treatCnvAsColumnsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatRowVecAsOne', treatRowVecAsOneDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'AsRowVectors', asRowVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, arrays, varargin{:});
colNums = iP.Results.colNums;
outputMode = validatestring(iP.Results.OutputMode, validOutputModes);
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCnvAsColumns = iP.Results.TreatCnvAsColumns;
treatRowVecAsOne = iP.Results.TreatRowVecAsOne;
asRowVectors = iP.Results.AsRowVectors;

% Make sure colNums are column vectors
if ~(ischar(colNums) || isstring(colNums))
    colNums = force_column_vector(colNums);
end

% Check if colNums is compatible with arrays
if iscell(colNums)
    if ~iscell(arrays) || treatCellAsArray
        error(['colNums cannot be multiple vectors ', ...
                'if arrays is not a cell array!']);
    elseif numel(arrays) ~= numel(colNums)
        error(['If colNums are multiple vectors, it must have ', ...
                'the same number of vectors as the number of arrays!'])
    end
end

%% Preparation
% If arrays is empty, return empty outputs
if isempty(arrays)
    varargout = cell(1, nargout);
    return
end

% Count the number of columns for each array
if ~iscell(arrays) || treatCellAsArray
    nColumns = size(arrays, 2);
elseif treatCnvAsColumns && iscellnumericvector(arrays)
    nColumns = count_vectors(arrays);
else
    nColumns = array_fun(@count_vectors, arrays);
end

% Make sure provided column numbers are within range
if iscell(colNums)
    % Get the maximum column numbers
    maxColNumbers = array_fun(@max, colNums);

    % If any maximum column number is greater than 
    %   the corresponding number of columns, return error
    if any(maxColNumbers > nColumns)
        indProblematic = find(maxColNumbers > nColumns);
        for idx = 1:length(indProblematic)
            error(['The maximum column number requested exceeded ', ...
                    'the minimum number of columns for array number %d!'], ...
                    indProblematic(idx));
        end
    end    
elseif isnumeric(colNums)
    % Get the maximum column number
    maxColNumber = max(colNums);

    % Get the minimum number of columns
    minNCols = min(min(min(nColumns)));

    % If the maximum column number is greater than 
    %   the minimum number of columns, return error
    if maxColNumber > minNCols
        error(['The maximum column number cannot exceed ', ...
                'the minimum number of columns!']);
    end
end

% Modify colNums to be compatible with arrays
%   based on whether column numbers are provided
if iscell(colNums) || isnumeric(colNums)
    % Only need to match dimensions if there are multiple arrays
    if iscell(arrays) && ~treatCellAsArray && ...
            ~(treatCnvAsColumns && iscellnumericvector(arrays))
        % Force as a cell array if not already so
        colNums = force_column_cell(colNums);

        % Match the dimensions of the cell array with arrays
        colNums = match_dimensions(colNums, size(arrays));
    end
elseif ischar(colNums) || isstring(colNums)
    % Generate column numbers
    colNums = create_indices('IndexEnd', nColumns, ...
                            'ForceRowOutput', asRowVectors);
end

% Count the number of columns requested
nColNumbers = count_samples(colNums);

% Count the number of output arguments requested
% TODO: Pull this out to a function? Or just use nargoutchk?
%       [nOutputs, errorMessage] = check_nargout(nargout, outputMode, nColNumbers)
%       [nOutputs, errorMessage] = check_nargout(nargout, outputMode, maxNArgOut)
switch outputMode
    case 'multiple'
        % Return as many outputs as needed
        nOutputs = nargout;

        % If the number of outputs requested is greater than 
        %   the number of columns requested, return error
        if nOutputs > max(nColNumbers)
            error(['The number of outputs requested ', ...
                    'cannot exceed the number of columns requested!']);
        end
    case 'single'
        % Return a single output
        nOutputs = 1;

        % If more than one output requested, return error
        if nargout > 1
            error('There can only be one output under ''single'' mode!');
        end
    otherwise
        error('outputMode unrecognized!');
end

%% Extract columns
switch outputMode
    case 'multiple'
        if nargout > 0
            % Restrict the column numbers to nOutputs
            colNums = extract_subvectors(colNums, 'IndexEnd', nOutputs);
        end
        
        % Extract nOutputs columns
        varargout = extract_columns_helper(colNums, arrays, ...
                        treatCellAsArray, treatCnvAsColumns, asRowVectors);
    case 'single'
        % Extract all columns
        varargout{1} = extract_columns_helper(colNums, arrays, ...
                        treatCellAsArray, treatCnvAsColumns, asRowVectors);
    otherwise
        error('outputMode unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colExtracted = extract_columns_helper(colNums, arrays, ...
                            treatCellAsArray, treatCnvAsColumns, asRowVectors)
%% Extract the columns depending on the number of arrays

if treatCellAsArray || ~iscell(arrays) || ...
        (treatCnvAsColumns && iscellnumericvector(arrays))
    % Treat arrays as a single array
    colExtracted = array_fun(@(x) extract_from_one_array(x, arrays, ...
                                    treatCellAsArray, treatCnvAsColumns, ...
                                    asRowVectors), ...
                            colNums, 'UniformOutput', false);
else
    % Count the minimum number of columns requested
    nColNumbers = count_samples(colNums, 'ForceColumnOutput', false);
    nColsToExtract = min(min(min(nColNumbers)));

    % Create iColsToExtract
    if asRowVectors
        iColsToExtract = 1:nColsToExtract;
    else
        iColsToExtract = transpose(1:nColsToExtract);
    end

    % Treat each cell content as a different array
    colExtracted = array_fun(@(x) extract_specific_column(x, colNums, arrays, ...
                                    treatCellAsArray, treatCnvAsColumns, ...
                                    asRowVectors), ...
                            iColsToExtract, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colExtracted = extract_specific_column(iCol, colNums, arrays, ...
                            treatCellAsArray, treatCnvAsColumns, asRowVectors)
%% Extract a specific column (may be different for each array)

% Get the current column numbers to extract
colNumsThis = extract_elements(colNums, 'specific', 'Index', iCol);

% Match dimensions
colNumsThis = match_dimensions(colNumsThis, size(arrays));

% Extract this column from each array
colExtracted = ...
    array_fun(@(x, y) extract_from_one_array(x, y, ...
                    treatCellAsArray, treatCnvAsColumns, asRowVectors), ...
                num2cell(colNumsThis), arrays, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colExtracted = extract_from_one_array(colNum, array, ...
                            treatCellAsArray, treatCnvAsColumns, asRowVectors)
%% Extract a specific column from a specific array

if treatCellAsArray || ~iscell(array) || ...
        treatCnvAsColumns && iscellnumericvector(array)
    % Extract the column
    if treatCellAsArray || ~iscell(array)
        % Treat array as a non-cell array
        colExtracted = array(:, colNum);
    else
        % Treat each cell content as a different column of the same array
        colExtracted = array{colNum};
    end

    % Force as a row vector if requested
    if asRowVectors
        colExtracted = force_row_vector(colExtracted);
    end
else
    % Pass to extract_columns_helper again
    colExtracted = extract_columns_helper(colNum, array, ...
                            treatCellAsArray, treatCnvAsColumns, asRowVectors);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(colNums)
    % Match the dimensions of the cell array
    colNums = match_dimensions(colNums, size(arrays));
elseif isnumeric(colNums) || treatCellAsArray
    % Place colNums in a cell and repmat it as many times as 
    %   the number of arrays
    colNums = repmat({colNums}, size(arrays));
end

if isnumeric(colNums)
    nColNumbers = length(colNums);
elseif iscell(colNums)
    nColNumbers = cellfun(@length, colNums);
end

%                   must be a numeric array or a cell array of numeric arrays

@(x) assert(isnumeric(x) || iscellnumeric(x), ...
            ['arrays must be either a numeric array ', ...
                'or a cell array of numeric arrays!']));

if isnumeric(arrays) || treatCellAsArray
if iscell(arrays) && ~treatCellAsArray
if isnumeric(arrays) || treatCellAsArray
if isnumeric(arrays) || treatCellAsArray

%       cd/force_column_vector.m
% Transform arrays into cell arrays of column vectors
if nArrays == 1
    varargout{1} = force_column_vector(arrays(:, colNums));
else
    varargout{1} = ...
        cellfun(@(x, y) force_column_vector(x(:, y)), ...
                    arrays, colNums, ...
                    'UniformOutput', false);
end

if nArrays == 1
    varargout{1} = force_column_cell(arrays(:, colNums));
else
    varargout{1} = ...
        cellfun(@(x, y) force_column_cell(x(:, y)), ...
                    arrays, colNums, 'UniformOutput', false);
end

% If arrays is a cell array with one element, pull it out
%   Note: this will prevent cell arrays having numel == 1
%           when it is not considered an array
if ~treatCellAsArray && iscell(arrays) && numel(arrays) == 1
    arrays = arrays{1};
end

switch outputMode
    case 'multiple'
        % Extract columns
        varargout = cell(1, nOutputs);
        for iOutput = 1:nOutputs
            if isnumeric(arrays) || treatCellAsArray
                varargout{iOutput} = arrays(:, colNums(iOutput));
            else
                varargout{iOutput} = ...
                    cellfun(@(x, y) x(:, y(iOutput)), arrays, colNums, ...
                            'UniformOutput', false);
            end
        end
    case 'single'
        % Transform arrays into cell arrays of column vectors
        if isnumeric(arrays) || treatCellAsArray
            varargout{1} = extract_columns_helper(arrays, colNums);
        else
            varargout{1} = ...
                cellfun(@(x, y) extract_columns_helper(x, y), ...
                            arrays, colNums, 'UniformOutput', false);
        end
    otherwise
        error('outputMode unrecognized!');
end

function colExtracted = extract_columns_helper(array, colNums)
% Extract the columns for a single array

% Extract as a cell array
colExtracted = arrayfun(@(x) array(:, x), colNums, 'UniformOutput', false);

% Force as a column cell array
colExtracted{iCol} = force_column_cell(colExtractedThis);

% Count the number of arrays
if isnumeric(arrays) || treatCellAsArray
    nArrays = 1;
else
    nArrays = numel(arrays);
end

if isnumeric(arrays) || treatCellAsArray
    colNums = transpose(1:nColumns);
else
    colNums = arrayfun(@(x) transpose(1:x), nColumns, ...
                        'UniformOutput', false);
end

colExtracted = cell(nColsToExtract, 1);
for iCol = 1:nColsToExtract
     colExtracted{iCol} = ...
        cellfun(@(x, y) x(:, y(iCol)), arrays, colNums, ...
                'UniformOutput', false);
end

function colExtracted = extract_columns_helper(arrays, colNums, ...
                            nColsToExtract, treatCellAsArray, treatCnvAsColumns)


colExtracted = cell(nColsToExtract, 1);
for iCol = 1:nColsToExtract
    % Get the current column numbers to extract
    colNumsThis = extract_elements(colNums, 'specific', 'Index', iCol);

    % Extract this column from each array
    colExtracted{iCol} = ...
        cellfun(@(x) extract_columns_helper(x, y, ...
                        treatCellAsArray, treatCnvAsColumns)), ...
                arrays, colNumsThis, 'UniformOutput', false);
end

% Treat arrays as a non-cell array
colExtracted = arrayfun(@(x) arrays(:, x), colNums, ...
                        'UniformOutput', false);

addRequired(iP, 'arrays', ...
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['arrays must be either a numeric array ', ...
                    'or a cell array!']));

% Count the minimum number of columns requested
nColNumbers = count_samples(colNum);
nColsToExtract = min(nColNumbers);

% Create iColsToExtract
if asRowVectors
    iColsToExtract = 1:nColsToExtract;
else
    iColsToExtract = transpose(1:nColsToExtract);
end

% Treat each cell content as a different array
colExtracted = arrayfun(@(x) extract_specific_column(x, colNum, array, ...
                                treatCellAsArray, treatCnvAsColumns, ...
                                asRowVectors), ...
                        iColsToExtract, 'UniformOutput', false);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
