function isInteger = isaninteger (array)
%% Returns whether each element of an array is an integer
% Usage: isInteger = isaninteger (array)
% Explanation:
%       Tests whether each element is an integer, regardless of type
%           Note: the built-in isinteger function will return true
%                   only if the element is of an int type
% Example(s):
%       isaninteger([Inf, 30, 30.5, -Inf, -100])
%       isaninteger({'sets'})
% Outputs:
%       isInteger   - whether each element is an integer
%                   specified as a logical array
% Arguments:    
%       array       - an array to check
%
% Used by:
%       cd/extract_elements.m
%       cd/ispositiveintegerarray.m
%       cd/ispositiveintegerscalar.m
%       cd/ispositiveintegervector.m
%       ~/Settings_Matlab/function_template.m

% File History:
% 2018-10-24 Adapted from https://www.mathworks.com/matlabcentral/answers/...
%                       326094-how-can-i-test-if-an-input-value-is-an-integer
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
% Check if the array is numeric
if ~isnumeric(array)
    isInteger = false;
    return
end

% Check if each element is an integer
%   Note: Since floor(Inf) == Inf and floor(-Inf) == -Inf, 
%           the finite check is necessary
isInteger = isfinite(array) & array == floor(array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%