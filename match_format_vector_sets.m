function [vecs1, vecs2] = match_format_vector_sets (vecs1, vecs2, varargin)
%% Matches two sets of vectors so that they are both cell arrays of the same number of column vectors
% Usage: [vecs1, vecs2] = match_format_vector_sets (vecs1, vecs2, varargin)
% Explanation:
%       This function takes two sets of vectors as input arguments
%       If there is more than one vector in one of the sets, 
%           this function forces both sets as column cell arrays of the same 
%           number of vectors so that cellfun() can be used
%       Otherwise, this function puts character vectors or numeric vectors 
%           in a cell array only if 'ForceCellOutputs' is set to true
%       cf. match_array_counts.m
%
% Example(s):
%       [a, b] = match_format_vector_sets((2:5)', 1:5)
%       [a, b] = match_format_vector_sets({1:5, 2:6}, 1:5)
%       [a, b] = match_format_vector_sets({1:5, [2:6]'}, 1:5)
%       [a, b] = match_format_vector_sets([[1:5]', [2:6]'], [1:5]')
%       [a, b] = match_format_vector_sets({'yes'}, 1:5)
%       [a, b] = match_format_vector_sets('yes', magic(3))
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
%       cd/match_dimensions.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/compute_default_sweep_info.m
%       cd/compute_relative_event_times.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/compute_sampsizepwr.m
%       cd/construct_fullpath.m
%       cd/create_indices.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/find_window_endpoints.m
%       cd/m3ha_plot_individual_traces.m
%       cd/match_and_combine_vectors.m
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
% TODO: Include the option to not force as column cell arrays
%           i.e., match 2D cell arrays
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs1', ...
    @(x) assert(isnum(x) || ischar(x) || isstring(x) || iscell(x), ...
                ['vecs1 must be either a numeric array, ', ...
                    'a character array, a string array or a cell array!']));
addRequired(iP, 'vecs2', ...
    @(x) assert(isnum(x) || ischar(x) || isstring(x) || iscell(x), ...
                ['vecs2 must be either a numeric array, ', ...
                    'a character array, a string array or a cell array!']));

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

% Read from the Input Parser
parse(iP, vecs1, vecs2, varargin{:});
forceCellOutputs = iP.Results.ForceCellOutputs;
matchVectors = iP.Results.MatchVectors;
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
treatCharAsScalar = iP.Results.TreatCharAsScalar;

%% Do the job
% If the vecs1 or vecs2 is a numeric vector, make sure it is a column vector
[vecs1, vecs2] = ...
    argfun(@(x) apply_or_return(isnumericvector(x), ...
                                @force_column_vector, x), ...
            vecs1, vecs2);

% If there are more than one vectors in either vecs1 or vecs2, 
%   put things in a format so cellfun can be used
if is_vector(vecs1, treatCellNumAsArray, treatCellStrAsArray, ...
                treatCellAsArray, treatCharAsScalar) && ...
        is_vector(vecs2, treatCellNumAsArray, treatCellStrAsArray, ...
                    treatCellAsArray, treatCharAsScalar)
    % Match vectors if requested
    if matchVectors
        [vecs1, vecs2] = match_format_vectors(vecs1, vecs2);
    end
else
    % Force vecs1/vecs2 to become column cell arrays of column vectors
    [vecs1, vecs2] = argfun(@force_column_cell, vecs1, vecs2);

    % Find the maximum number of rows
    maxVecs = max(numel(vecs1), numel(vecs2));

    % Match up the vector counts
    % TODO: Incorporate comparison into match_dimensions.m
    [vecs1, vecs2] = ...
        argfun(@(x) match_dimensions(x, [maxVecs, 1]), vecs1, vecs2);
    
    % Match vectors if requested
    if matchVectors
        [vecs1, vecs2] = cellfun(@(x, y) match_format_vectors(x, y), ...
                                vecs1, vecs2, 'UniformOutput', false);
    end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isVector = is_vector(vecs1, treatCellNumAsArray, ...
                                treatCellStrAsArray, treatCellAsArray, ...
                                treatCharAsScalar)
%% Tests whether
% TODO: Use this in other functions?

isVector = isnumericvector(vecs1) || ...
        ischar(vecs1) && treatCharAsScalar || ...
        iscellnumeric(vecs1) && treatCellNumAsArray || ...
        iscellstr(vecs1) && treatCellStrAsArray || ...
        iscell(vecs1) && treatCellAsArray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force vecs1/vecs2 to become column vectors
[vecs1, vecs2] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', false), ...
            vecs1, vecs2);

[vecs1, vecs2] = argfun(@force_column_cell, vecs1, vecs2);

% Count the vectors
[nVecs1, nVecs2] = ...
    argfun(@(x) count_vectors(x, 'TreatCellAsArray', treatCellAsArray, ...
                    'TreatCellStrAsArray', treatCellStrAsArray), ...
            vecs1, vecs2);

% Find the maximum number of vector counts
maxVecs = max(nVecs1, nVecs2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
