function isPositiveVector = ispositivevector (x)
%% Returns whether an input is a positive vector
% Usage: isPositiveVector = ispositivevector (x)
% Explanation:
%       Tests whether the input is a positive vector
% Example(s):
%       ispositivevector(1:10)
%       ispositivevector(magic(3))
%       ispositivevector(-1:3)
% Outputs:
%       isPositiveVector
%                       - whether the input is a positive vector
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%
% Used by:
%       cd/match_time_info.m
%       cd/parse_pulse.m

% File History:
% 2018-12-17 Modified from ispositiveintegervector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isPositiveVector = isnumeric(x) && isvector(x) && all(x > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%