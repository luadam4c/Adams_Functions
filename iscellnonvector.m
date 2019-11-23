function isCellNonVector = iscellnonvector (x)
%% Returns whether an input is a cell array of non-cell vectors (may be empty)
% Usage: isCellNonVector = iscellnonvector (x)
% Explanation:
%       Tests whether the input is a cell array of non-cell vectors
% Example(s):
%       iscellnonvector({1:10, 2:20})
%       iscellnonvector({magic(3), 2:20})
%       iscellnonvector({'sets', 'lasts'})
%       iscellnonvector({{1:10, 2:20}, {1:10, 2:20}})
% Outputs:
%       isCellNonVector    - whether the input is a cell array of vectors
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires:
%
% Used by:

% File History:
% 2019-01-04 Adapted from iscellnumericvector.m
% 2019-01-10 Updated definition so that cell arrays are not vectors
% 2019-01-18 Renamed iscellvector -> iscellnonvector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isCellNonVector = iscell(x) && ...
                all(all(all(cellfun(@(x) isvector(x) && ~iscell(x), x))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
