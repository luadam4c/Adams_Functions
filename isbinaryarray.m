function isBinaryArray = isbinaryarray (x)
%% Returns whether the input is a binary array
% Usage: isBinaryArray = isbinaryarray (x)
% Explanation:
%       Tests whether each element of the input is a binary array
% Example(s):
%       isbinary([])
%       isbinary(false(10, 1))
%       isbinary(2)
%       isbinary(1)
%       isbinary([0 0 1])
% Outputs:
%       isbinaryarray   - whether the input is a binary array
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Used by:
%       cd/force_logical.m
%       cd/isbinaryscalar.m

% File History:
% 2018-12-12 Pulled from isbinaryscalar.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isBinaryArray = islogical(x) || isnumeric(x) && all(x == 0 | x == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%