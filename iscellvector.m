function isCellVector = iscellvector (x)
%% Returns whether an input is a cell array of non-cell vectors (may be empty)
% Usage: isCellVector = iscellvector (x)
% Explanation:
%       Tests whether the input is a cell array of non-cell vectors
% Example(s):
%       iscellvector({1:10, 2:20})
%       iscellvector({magic(3), 2:20})
%       iscellvector({'sets', 'lasts'})
%       iscellvector({{1:10, 2:20}, {1:10, 2:20}})
% Outputs:
%       isCellVector    - whether the input is a cell array of vectors
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires:
%
% Used by:
%       cd/compute_average_data.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/extract_columns.m

% File History:
% 2019-01-04 Adapted from iscellnumericvector.m
% 2019-01-10 Updated definition so that cell arrays are not vectors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isCellVector = iscell(x) && ...
                all(all(all(cellfun(@(x) isvector(x) && ~iscell(x), x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%