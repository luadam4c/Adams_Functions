function vectorsCell = force_column_cell (vectorsOrig, varargin)
%% Transforms a row cell array or a non-cell array to a column cell array of non-cell vectors
% Usage: vectorsCell = force_column_cell (vectorsOrig, varargin)
% Explanation:
%       This is an attempt to standardize the way multiple vectors are stored
%           -- always as column cell vectors
%       1. A row cell vector is converted to a column cell vector
%       2. A non-vector cell array is transformed to a column cell vector
%           of column cell arrays
%           However, if 'ToLinear' is true, it will simply be linearized 
%           as a column cell vector
%       3. An empty numeric array or a character array are placed in a cell array
%       4. A numeric vector is forced as a column vector
%           (force_column_vector.m is used), then placed in a cell array
%       5. A numeric non-vector array is transformed to a column cell vector
%           of column numeric vectors based on the first dimension
%
% Example(s):
%       load_examples;
%       force_column_cell(myCellNumeric2D)
%       force_column_cell(myCellRowVecs)
%       force_column_cell(myCellStr2D)
%       force_column_cell(myCellStr2D, 'ToLinear', true)
%       force_column_cell(myNumeric2D)
%       force_column_cell(myNumeric3D)
%
% Outputs:
%       vectorsCell - vectors as a column cell array
%                   specified as a column cell array
%
% Arguments:
%       vectorsOrig - original vectors
%                   Note: If an array, each column is a vector 
%                           to be placed in a cell
%                   must be a numeric array or a cell array 
%                       or a character vector
%       varargin    - 'ToLinearize': whether to linearize a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RowInstead': whether to force as row vector instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatVectorAsElement': whether to treat numeric vector
%                                           as an element
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_columns.m
%       cd/force_column_vector.m
%
% Used by:
%       cd/atfwrite.m
%       cd/all_dependent_functions.m
%       cd/combine_data_from_same_slice.m
%       cd/combine_variables_across_tables.m
%       cd/compute_rms_error.m
%       cd/compute_sampsizepwr.m
%       cd/combine_swd_sheets.m
%       cd/construct_fullpath.m
%       cd/create_average_time_vector.m
%       cd/create_indices.m
%       cd/create_subplots.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/extract_columns.m
%       cd/extract_fields.m
%       cd/filter_and_extract_pulse_response.m
%       cd/find_pulse_endpoints.m
%       cd/force_row_cell.m
%       cd/force_column_vector.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_individual_traces.m
%       cd/match_format_vector_sets.m
%       cd/parse_atf_swd.m
%       cd/parse_iox.m
%       cd/parse_multiunit.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_all_abfs.m
%       cd/plot_protocols.m
%       cd/plot_struct.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%       cd/plot_window_boundaries.m
%       cd/read_lines_from_file.m
%       cd/run_neuron.m
%       cd/save_all_zooms.m
%       cd/struct2arglist.m
%       cd/test_normality.m
%       cd/test_var_difference.m
%       cd/vertcat_spreadsheets.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-10-19 Now accepts character vectors
% 2018-10-27 Now places empty numeric arrays in a cell array
% 2018-12-18 Now defaults 2D cell arrays to be separated by columns
%               added 'ToLinearize' (default == 'false')
% 2018-12-19 Now uses extract_columns.m
% 2019-01-04 Fixed bug
% 2019-01-09 Now forces string arrays to become cell arrays of character vectors
% 2019-01-13 Added 'RowInstead' as an optional argument
% 2019-01-22 Now makes the vector format consistent
% 

%% Default values for optional arguments
toLinearizeDefault = false;     % whether to linearize a nonvector cell array
rowInsteadDefault = false;      % whether to force as row vector instead
treatVectorAsElementDefault = true;    % treat vectors as element by default

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
addRequired(iP, 'vectorsOrig');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ToLinearize', toLinearizeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RowInstead', rowInsteadDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatVectorAsElement', treatVectorAsElementDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectorsOrig, varargin{:});
toLinearize = iP.Results.ToLinearize;
rowInstead = iP.Results.RowInstead;
treatVectorAsElement = iP.Results.TreatVectorAsElement;

%% Do the job
if isempty(vectorsOrig) || iscell(vectorsOrig) && ...
        (rowInstead && isrow(vectorsOrig) || ...
        ~rowInstead && iscolumn(vectorsOrig))
    % Place in a cell array if not already so
    if ~iscell(vectorsOrig)
        vectorsCell = {vectorsOrig};
    else
        vectorsCell = vectorsOrig;
    end
elseif ischar(vectorsOrig) || isstring(vectorsOrig)
    % Convert to a cell array of character arrays
    vectorsCell = cellstr(vectorsOrig);

    % Pass to this function again
    vectorsCell = force_column_cell(vectorsCell, 'ToLinearize', toLinearize, ...
                                    'RowInstead', rowInstead);
elseif iscell(vectorsOrig) && (isvector(vectorsOrig) || toLinearize)
    % Reassign as a column
    vectorsCell = vectorsOrig(:);

    % Transpose to a row if requested
    if rowInstead
        vectorsCell = transpose(vectorsCell);
    end
elseif ~iscell(vectorsOrig) || ...
        iscell(vectorsOrig) && ~isvector(vectorsOrig) && ~toLinearize
    % Force any non-cell vector as a column vector
    if ~iscell(vectorsOrig)
        % If vectorsOrig is a row vector, columns will be extracted as row vectors
        if isrow(vectorsOrig)
            asRowVectors = true;
        else
            asRowVectors = false;
        end
        
        % Force as a column vector
        vectorsOrig = force_column_vector(vectorsOrig, 'IgnoreNonVectors', true);
    else
        % Columns will be extract as column vectors
        asRowVectors = false;
    end

    % Extract as a cell array or use num2cell
    if treatVectorAsElement
        vectorsCell = extract_columns(vectorsOrig, 'all', ...
                                'OutputMode', 'single', 'TreatCellAsArray', true, ...
                                'AsRowVectors', asRowVectors);
    else
        vectorsCell = num2cell(vectorsOrig);
    end

    % Pass to this function again
    vectorsCell = force_column_cell(vectorsCell, 'ToLinearize', toLinearize, ...
                                    'RowInstead', rowInstead);
else
    % Should not occur
    error('Code logic error!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

nVectors = size(vectorsOrig, 2);

% Reassign as a column
vectorsCell = vectorsCell(:);

if iscell(vectorsOrig) && ~iscolumn(vectorsOrig)

%       cd/count_vectors.m
% Count the number of vectors
nVectors = count_vectors(vectorsOrig);

% Extract as a cell array
vectorsCell = arrayfun(@(x) vectorsOrig(:, x), ...
                        transpose(1:nVectors), ...
                        'UniformOutput', false);

elseif ischar(vectorsOrig) || isnum(vectorsOrig) && isempty(vectorsOrig)
    % Place in a cell array
    vectorsCell = {vectorsOrig};

%       cd/isnum.m
addRequired(iP, 'vectorsOrig', ...
    @(x) isnum(x) || iscell(x) || ischar(x) || isstring(x));

% Extract as a cell array
vectorsCell = extract_columns(vectorsOrig, 'all', ...
                        'OutputMode', 'single', 'TreatCellAsArray', true, ...
                        'AsRowVectors', rowInstead);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
