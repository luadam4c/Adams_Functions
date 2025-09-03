function [vecs1, vecs2] = match_format_vector_sets (vecs1, vecs2, varargin)
%% Matches two sets of vectors so that they are both cell arrays of the same number of column vectors
% Usage: [vecs1, vecs2] = match_format_vector_sets (vecs1, vecs2, varargin)
% Explanation:
%       This function takes two sets of vectors as input arguments
%       If one of the sets is a cell array of vectors, 
%           this function forces both sets as cell arrays of the same 
%           number of vectors so that cellfun() can be used
%       Otherwise, this function puts character vectors or numeric vectors 
%           in a cell array only if 'ForceCellOutputs' is set to true
%       cf. match_array_counts.m
%
% Example(s):
%       [a, b] = match_format_vector_sets((2:5)', 1:5)
%       [a, b] = match_format_vector_sets((2:5)', 1:5, 'MatchVectors', true)
%       [a, b] = match_format_vector_sets({1:5, 2:6}, 1:5)
%       [a, b] = match_format_vector_sets({1:5, 2:6}, 1:5, 'MatchVectors', true)
%       [a, b] = match_format_vector_sets({1:5, 2:6; [], 5:-1:1}, 1:5)
%       [a, b] = match_format_vector_sets({1:5, [2:6]'}, 1:5)
%       [a, b] = match_format_vector_sets([[1:5]', [2:6]'], [1:5]')
%       [a, b] = match_format_vector_sets({'yes'}, 1:5)
%       [a, b] = match_format_vector_sets({'yes'}, 1:5, 'TreatRowVecAsOne', false)
%       [a, b] = match_format_vector_sets('yes', magic(3))
%       [a, b] = match_format_vector_sets('apple', 1:5)
%       [a, b] = match_format_vector_sets('apple', 1:5, 'TreatRowVecAsOne', false)
%       [a, b] = match_format_vector_sets('apple', 'banana')
%       [a, b] = match_format_vector_sets('apple', '')
%       [a, b] = match_format_vector_sets(@lines, magic(3))
%       [a, b] = match_format_vector_sets({'apple'}, '')
%
% Outputs:
%       vecs1       - new first set of vectors
%                   must be a numeric vector, a character array or a cell array
%       vecs2       - new second set of vectors
%                   must be a numeric vector, a character array or a cell array
%
% Arguments:
%       vecs1       - first set of vectors
%                   must be a numeric array, a character array or a cell array
%       vecs2       - second set of vectors
%                   must be a numeric array, a character array or a cell array
%       varargin    - 'ForceCellOutputs': whether to force outputs as 
%                                           cell arrays
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchVectors': whether to match vectors within sets
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
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
%                   default == false
%                   - 'TreatCharAsScalar': whether to treat character arrays 
%                                           as scalars
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'TreatRowVecAsOne': whether to treat row vectors
%                                           as a single vector
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/apply_or_return.m
%       cd/argfun.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/force_string_end.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%       cd/isnumericvector.m
%       cd/match_column_count.m
%       cd/match_dimensions.m
%       cd/match_format_vectors.m
%       cd/match_row_count.m
%
% Used by:
%       cd/compute_activation_profile.m
%       cd/compute_default_sweep_info.m
%       cd/compute_relative_event_times.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/compute_sampsizepwr.m
%       cd/construct_fullpath.m
%       cd/create_indices.m
%       cd/create_time_vectors.m
%       cd/create_trace_plot_movie.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_pulse_response_endpoints.m
%       cd/find_window_endpoints.m
%       cd/read_neuron_outputs.m
%       cd/match_and_combine_vectors.m
%       cd/nan_except.m
%       cd/parse_phase_info.m
%       cd/parse_pleth_trace.m
%       cd/parse_pulse.m
%       cd/plot_fitted_traces.m
%       cd/plot_table_parallel.m
%       cd/plot_vertical_line.m
%       cd/plot_vertical_shade.m
%       cd/plot_traces.m
%       cd/test_var_difference.m
%       cd/transform_vectors.m

% File History:
% 2018-10-28 Adapted from code in find_window_endpoints.m 
%               and match_vector_counts.m
% 2018-10-31 Now uses isnumericvector.m and apply_or_return.m
% 2019-01-03 Added 'MatchVectors' as an optional argument
% 2019-01-04 Added 'TreatCellAsArray' (default == 'false')
% 2019-01-04 Added 'TreatCellStrAsArray' (default == 'true')
% 2019-01-04 Added 'TreatCharAsScalar' (default == 'true')
% 2019-01-04 Fixed bug
% 2019-10-03 Added 'TreatCellNumAsArray' as an optional argument
% 2019-12-01 Now does not necessarily force as column cell arrays
%               i.e., match 2D cell arrays
% 2019-12-23 Now lets the vectors be any type
% 2020-01-02 Added 'TreatRowVecAsOne' as an optional argument
% TODO: Add 'ForceColumnOutput' as an optional argument
% TODO: Accept more than two vector sets
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputsDefault = false;    % don't force as cell array by default
matchVectorsDefault = false;        % don't match vectors by default
treatCellAsArrayDefault = false;    % treat cell arrays as an array by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = false; % treat cell arrays of character arrays
                                    %   as an array by default
treatCharAsScalarDefault = true;% treat character arrays as scalars by default
treatRowVecAsOneDefault = true;     % treat row vectors as one vector by default

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
addRequired(iP, 'vecs1');
addRequired(iP, 'vecs2');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutputs', forceCellOutputsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchVectors', matchVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCharAsScalar', treatCharAsScalarDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatRowVecAsOne', treatRowVecAsOneDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs1, vecs2, varargin{:});
forceCellOutputs = iP.Results.ForceCellOutputs;
matchVectors = iP.Results.MatchVectors;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
treatCharAsScalar = iP.Results.TreatCharAsScalar;
treatRowVecAsOne = iP.Results.TreatRowVecAsOne;

%% Do the job
% If there are more than one vectors in either vecs1 or vecs2, 
%   put things in a format so cellfun can be used
if is_vector(vecs1, treatCellNumAsArray, treatCellStrAsArray, ...
                treatCellAsArray, treatRowVecAsOne) && ...
        is_vector(vecs2, treatCellNumAsArray, treatCellStrAsArray, ...
                    treatCellAsArray, treatRowVecAsOne)
    % Match vectors if requested
    if matchVectors
        % TODO: Pass in treatCharAsScalar && treatRowVecAsOne
        [vecs1, vecs2] = match_format_vectors(vecs1, vecs2);
    end

    % Force outputs to be cell arrays if requested
    if forceCellOutputs
        if ~iscell(vecs1)
            vecs1 = {vecs1};
        end
        if ~iscell(vecs2)
            vecs2 = {vecs2};
        end
    end
elseif ~forceCellOutputs && ...
            is_concatable(vecs1, treatCellNumAsArray, ...
                                treatCellStrAsArray, treatCellAsArray) && ...
            is_concatable(vecs2, treatCellNumAsArray, ...
                                treatCellStrAsArray, treatCellAsArray)
    % Force any row vector to be a column vector
    if treatRowVecAsOne
        [vecs1, vecs2] = ...
            argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
                     vecs1, vecs2);
    end

    % Find the maximum number of columns
    maxColumns = max(size(vecs1, 2), size(vecs2, 2));

    % Force vecs1/vecs2 to have the same number of columns
    [vecs1, vecs2] = ...
        argfun(@(x) match_column_count(x, maxColumns), vecs1, vecs2);

    % Match vectors if requested
    if matchVectors
        % Find the maximum number of rows
        maxRows = max(size(vecs1, 1), size(vecs2, 1));

        % Force vecs1/vecs2 to have the same number of rows
        [vecs1, vecs2] = argfun(@(x) match_row_count(x, maxRows), vecs1, vecs2);
    end
else
    % Force vecs1/vecs2 to become cell arrays of vectors
    [vecs1, vecs2] = ...
        argfun(@(x) force_column_cell(x, 'IgnoreNonVectors', true, ...
                                    'TreatRowVecAsOne', treatRowVecAsOne), ...
                vecs1, vecs2);

    % Find the maximum number of rows
    maxRows = max(size(vecs1, 1), size(vecs2, 1));

    % Find the maximum number of columns
    maxColumns = max(size(vecs1, 2), size(vecs2, 2));
    
    % Match up the vector counts
    % TODO: Incorporate comparison into match_dimensions.m
    [vecs1, vecs2] = ...
        argfun(@(x) match_dimensions(x, [maxRows, maxColumns]), vecs1, vecs2);
    
    % Match vectors if requested
    % TODO: Pass in treatCharAsScalar && treatRowVecAsOne
    if matchVectors
        [vecs1, vecs2] = cellfun(@(x, y) match_format_vectors(x, y), ...
                                vecs1, vecs2, 'UniformOutput', false);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isVector = is_vector(array, ...
                                treatCellNumAsArray, treatCellStrAsArray, ...
                                treatCellAsArray, treatRowVecAsOne)
%% Tests whether an array is considered a vector
% TODO: Use this in other functions?
% TODO: Test is_vector?

isVector = (isempty(array) || iscolumn(array) || ...
                    isrow(array) && treatRowVecAsOne) && ...
                is_non_cell(array, treatCellNumAsArray, ...
                            treatCellStrAsArray, treatCellAsArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isConcatable = is_concatable (array, ...
                                treatCellNumAsArray, treatCellStrAsArray, ...
                                treatCellAsArray)
%% Tests whether an array can be concatenated the usual way
%   Note: character arrays and function handles

isConcatable = ~ischar(array) && ~isa(array, 'function_handle') && ...
                        is_non_cell(array, treatCellNumAsArray, ...
                            treatCellStrAsArray, treatCellAsArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isNonCell = is_non_cell (array, treatCellNumAsArray, ...
                                    treatCellStrAsArray, treatCellAsArray)
%% Tests whether an array is considered a non-cell array
% TODO: Use this in other functions?

isNonCell = ~iscell(array) || ...
            iscellnumeric(array) && treatCellNumAsArray || ...
            iscellstr(array) && treatCellStrAsArray || ...
            iscell(array) && treatCellAsArray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
