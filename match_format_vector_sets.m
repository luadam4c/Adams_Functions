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
%       [a, b] = match_format_vector_sets({1:5, 2:6}, 1:5)
%       [a, b] = match_format_vector_sets({1:5, [2:6]'}, 1:5)
%       [a, b] = match_format_vector_sets([[1:5]', [2:6]'], [1:5]')
%       [a, b] = match_format_vector_sets('yes', 1:5)
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
%
% Requires:
%       cd/apply_or_return.m
%       cd/argfun.m
%       cd/force_column_cell.m
%       cd/force_column_numeric.m
%       cd/iscellnumeric.m
%       cd/isnumericvector.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/compute_default_sweep_info.m
%       cd/compute_residuals.m
%       cd/compute_rms_error.m
%       cd/construct_fullpath.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/find_pulse_response_endpoints.m
%       cd/find_window_endpoints.m
%       cd/m3ha_plot_individual_traces.m
%       cd/plot_traces.m

% File History:
% 2018-10-28 Adapted from code in find_window_endpoints.m 
%               and match_vector_counts.m
% 2018-10-31 Now uses isnumericvector.m and apply_or_return.m
% TODO: Include the option to not force as column cell arrays
%           i.e., match 2D cell arrays
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputsDefault = false;    % don't force as cell array by default

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
    @(x) assert(isnumeric(x) || ischar(x) || iscell(x), ...
                ['vecs1 must be either a numeric array, ', ...
                    'a character array or a cell array!']));
addRequired(iP, 'vecs2', ...
    @(x) assert(isnumeric(x) || ischar(x) || iscell(x), ...
                ['vecs2 must be either a numeric array, ', ...
                    'a character array or a cell array!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutputs', forceCellOutputsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs1, vecs2, varargin{:});
forceCellOutputs = iP.Results.ForceCellOutputs;

%% Do the job
% If the vecs1 or vecs2 is a numeric vector, make sure it is a column vector
[vecs1, vecs2] = ...
    argfun(@(x) apply_or_return(isnumericvector(x), ...
                                @force_column_numeric, x), ...
            vecs1, vecs2);

% If there are more than one vectors in either vecs1 or vecs2, 
%   put things in a format so cellfun can be used
if (isnumericvector(vecs1) || ischar(vecs1)) && ...
    (isnumericvector(vecs2) || ischar(vecs2))
    % If both inputs are either a numeric vector (could be empty) 
    %   or a character array, force outputs to be cell arrays if requested
    if forceCellOutputs
        vecs1 = {vecs1};
        vecs2 = {vecs2};
    end
else
    % Force vecs1/vecs2 to become 
    %   column cell arrays of column numeric vectors
    [vecs1, vecs2] = argfun(@force_column_cell, vecs1, vecs2);

    % Find the maximum number of rows
    maxVecs = max(numel(vecs1), numel(vecs2));

    % Match up the vector counts
    % TODO: Incorporate comparison into match_dimensions.m
    [vecs1, vecs2] = ...
        argfun(@(x) match_dimensions(x, [maxVecs, 1]), vecs1, vecs2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[vecs1, vecs2] = ...
    argfun(@(x) force_column_numeric(x, 'IgnoreNonVectors', true), ...
            vecs1, vecs2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
