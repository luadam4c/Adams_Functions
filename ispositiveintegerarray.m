function isPositiveIntegerArray = ispositiveintegerarray (x)
%% Returns whether an input is a positive integer array
% Usage: isPositiveIntegerArray = ispositiveintegerarray (x)
% Explanation:
%       Tests whether the input is a positive integer array
% Example(s):
%       ispositiveintegerarray(1:10)
%       ispositiveintegerarray(magic(3))
%       ispositiveintegerarray(-1:3)
% Outputs:
%       isPositiveIntegerArray
%                       - whether the input is a positive integer array
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%       cd/create_error_for_nargin.m
%       cd/isaninteger.m
%
% Used by:
%       cd/extract_columns.m

% File History:
% 2019-01-12 Modified from ispositiveintegervector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isPositiveIntegerArray = isnumeric(x) && all(all(all(isaninteger(x) & x > 0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%