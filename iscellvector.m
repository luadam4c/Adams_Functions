function isCellVector = iscellvector (x)
%% Returns whether an input is a cell array of vectors (may be empty)
% Usage: isCellVector = iscellvector (x)
% Explanation:
%       Tests whether the input is a cell array of vectors
%
% Example(s):
%       iscellvector({1:10, 2:20})
%       iscellvector({magic(3), 2:20})
%       iscellvector({'sets', 'lasts'})
%       iscellvector({{1:10, 2:20}, {1:10, 2:20}})
%
% Outputs:
%       isCellVector    - whether the input is a cell array of vectors
%                       specified as a logical scalar
%
% Arguments:    
%       x               - an input to check
%
% Requires:
%
% Used by:
%       cd/all_ordered_pairs.m
%       cd/compute_combined_data.m
%       cd/vecfun.m

% File History:
% 2019-01-18 Adapted from iscellnonvector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isCellVector = iscell(x) && ...
                all(all(all(cellfun(@(x) isempty(x) || isvector(x), x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
