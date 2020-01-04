function varargout = match_format_vectors (varargin)
%% Match the format of individual vectors by making them all column vectors with the same length
% Usage: varargout = match_format_vectors (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [a, b] = match_format_vectors(1:5, (2:6)')
%       [a, b] = match_format_vectors(1:5, 2)
%       [a, b] = match_format_vectors(1:5, 2:3)
%
% Outputs:
%       varargout   - matched outputs
%
% Arguments:
%       varargin    - vectors to be matched
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/match_row_count.m
%
% Used by:
%       cd/compute_time_window.m
%       cd/create_indices.m
%       cd/create_time_vectors.m

% File History:
% 2018-12-16 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < nargout
    disp('Cannot request more outputs than provided inputs!!');
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% Force as column vectors
vararginTransformed = cellfun(@force_column_vector, varargin, ...
                                'UniformOutput', false);

% Count the number of values in each vector
nValues = cellfun(@numel, varargin);

% Compute the maximum number of values
maxNValues = max(nValues);

% Match the number of rows to maxNValues for each vector
varargout = cellfun(@(x) match_row_count(x, maxNValues), ...
                    vararginTransformed, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
