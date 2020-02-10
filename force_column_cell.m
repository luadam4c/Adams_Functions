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
%       3. A non-matrix numeric array or a character array is 
%           simply placed in a cell array
%       4. A numeric non-vector array is transformed to a column cell vector
%           of column numeric vectors
%
% Example(s):
%       load_examples;
%       force_column_cell(myCellNumeric2D)
%       force_column_cell(myCellRowVecs)
%       force_column_cell(myCellStr2D)
%       force_column_cell(myCellStr2D, 'ToLinear', true)
%       force_column_cell(myNumeric2D)
%       force_column_cell(myNumeric3D)
%       force_column_cell(myCellStr2D)
%       force_column_cell(myCellStr2D, 'IgnoreNonVectors', true)
%       force_column_cell(1:5)
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
%                   - 'IgnoreNonVectors': whether to ignore non-vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatVectorAsElement': whether to treat numeric vector
%                                           as an element
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TreatRowVecAsOne': whether to treat row vectors
%                                           as a single vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_columns.m
%       cd/force_column_vector.m
%
% Used by:
%       cd/array_fun.m
%       cd/atfwrite.m
%       cd/all_dependent_functions.m
%       cd/combine_data_from_same_slice.m
%       cd/combine_variables_across_tables.m
%       cd/compute_sampsizepwr.m
%       cd/combine_swd_sheets.m
%       cd/construct_fullpath.m
%       cd/copy_into.m
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
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_decide_on_plot_vars.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_figure07.m
%       cd/match_format_vector_sets.m
%       cd/parse_atf_swd.m
%       cd/parse_current_family.m
%       cd/parse_iox.m
%       cd/parse_multiunit.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/plot_all_abfs.m
%       cd/plot_fitted_traces.m
%       cd/plot_repetitive_protocols.m
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
%       cd/update_figure_for_corel.m
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
% 2019-12-01 Added 'IgnoreNonVectors' as an optional argument
% 2020-01-02 Simplified code for non-matrix arrays
% 2020-01-02 Added 'TreatRowVecAsOne' as an optional argument
% 

%% Default values for optional arguments
toLinearizeDefault = false;     % whether to linearize a nonvector cell array
rowInsteadDefault = false;      % whether to force as row vector instead
ignoreNonVectorsDefault = false;	% don't ignore non-vectors by default
treatVectorAsElementDefault = true; % treat vectors as element by default
treatRowVecAsOneDefault = true;     % treat row vectors as one vector by default

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
addParameter(iP, 'IgnoreNonVectors', ignoreNonVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatVectorAsElement', treatVectorAsElementDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatRowVecAsOne', treatRowVecAsOneDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectorsOrig, varargin{:});
toLinearize = iP.Results.ToLinearize;
rowInstead = iP.Results.RowInstead;
ignoreNonVectors = iP.Results.IgnoreNonVectors;
treatVectorAsElement = iP.Results.TreatVectorAsElement;
treatRowVecAsOne = iP.Results.TreatRowVecAsOne;

%% Do the job
if isempty(vectorsOrig) || iscell(vectorsOrig) && ...
        (rowInstead && isrow(vectorsOrig) || ...
        ~rowInstead && iscolumn(vectorsOrig) || ...
        ignoreNonVectors && ~isvector(vectorsOrig))
    % Place in a cell array if not already so
    %   Otherwise, do nothing
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
                                    'RowInstead', rowInstead, ...
                                    'TreatRowVecAsOne', treatRowVecAsOne);
elseif iscell(vectorsOrig) && (isvector(vectorsOrig) || toLinearize)
    % Reassign as a column
    vectorsCell = vectorsOrig(:);

    % Transpose to a row if requested
    if rowInstead
        vectorsCell = transpose(vectorsCell);
    end
elseif ~iscell(vectorsOrig) || ...
        iscell(vectorsOrig) && ~isvector(vectorsOrig) && ~toLinearize
    % Convert to a cell array or a cell vector
    if ~iscell(vectorsOrig) && ~treatVectorAsElement
        % Convert to a cell array
        vectorsCell = num2cell(vectorsOrig);
    elseif ~iscell(vectorsOrig) && ...
            (iscolumn(vectorsOrig) || isrow(vectorsOrig) && treatRowVecAsOne)
        % Put in a cell array
        vectorsCell = {vectorsOrig};
    elseif ~iscell(vectorsOrig) && isrow(vectorsOrig) && ~treatRowVecAsOne
        % Convert to a cell array
        vectorsCell = num2cell(vectorsOrig);
    else
        % Either a non-cell non-vector array or a cell non-vector array
        %   Columns will be extract as column vectors and placed
        %   in a column cell array
        vectorsCell = extract_columns(vectorsOrig, 'all', ...
                        'OutputMode', 'single', 'TreatCellAsArray', true, ...
                        'TreatRowVecAsOne', treatRowVecAsOne, ...
                        'AsRowVectors', false);
    end

    % Pass to this function again (for resolving rowInstead and toLinearize)
    vectorsCell = force_column_cell(vectorsCell, 'ToLinearize', toLinearize, ...
                                    'RowInstead', rowInstead, ...
                                    'TreatRowVecAsOne', treatRowVecAsOne);
else
    % Should not occur
    error('Code logic error!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
