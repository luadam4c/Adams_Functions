function varargout = match_format_vectors (varargin)
%% Match the format of individual vectors by making them all column vectors with the same length
% Usage: varargout = match_format_vectors (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [a, b] = match_format_vectors(1:5, (2:6)')
%       [a, b] = match_format_vectors(1:5, (2:6)', 'RowInstead', true)
%       [a, b] = match_format_vectors(1:5, 2)
%       [a, b] = match_format_vectors(1:5, 2:3)
%
% Outputs:
%       varargout   - matched outputs
%
% Arguments:
%       varargin    - vectors to be matched
%                   - 'RowInstead': whether to force as row vector instead
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/compute_maximum_numel.m
%       cd/create_error_for_nargin.m
%       cd/extract_parameter_value_pairs.m
%       cd/first_matching_field.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%
% Used by:
%       cd/compute_gabab_conductance.m
%       cd/compute_time_window.m
%       cd/create_indices.m
%       cd/create_time_vectors.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/parse_ipsc.m

% File History:
% 2018-12-16 Created by Adam Lu

%% Default values for optional arguments
rowInsteadDefault = false;      % whether to force as row vector instead

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < nargout
    disp('Cannot request more outputs than provided inputs!!');
    error(create_error_for_nargin(mfilename));
end

%% Preparation
% % Look for 'RowInstead' in the argument list
[params, varargin] = extract_parameter_value_pairs(varargin);
if ~isempty(params)
    rowInstead = first_matching_field(params, {'RowInstead', 'rowInstead'});

    if isempty(rowInstead)
        rowInstead = rowInsteadDefault;
    end
else
    rowInstead = rowInsteadDefault;
end

% Decide on the dimension to match
if rowInstead
    dimToMatch = 2;
else
    dimToMatch = 1;
end
        
%% Do the job
% Force as column vectors
vararginTransformed = ...
    cellfun(@(x) force_column_vector(x, 'RowInstead', rowInstead), ...
            varargin, 'UniformOutput', false);

% Compute the maximum number of values over all vectors
maxNValues = compute_maximum_numel(vararginTransformed);

% Match the number of rows to maxNValues for each vector
varargout = cellfun(@(x) match_row_count(x, maxNValues, ...
                                        'DimToMatch', dimToMatch), ...
                    vararginTransformed, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of values in each vector
nValues = cellfun(@numel, vararginTransformed);

% Compute the maximum number of values
maxNValues = max(nValues);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
