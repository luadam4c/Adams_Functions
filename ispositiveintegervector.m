function isPositiveIntegerVector = ispositiveintegervector (x)
%% Returns whether an input is a positive integer vector
% Usage: isPositiveIntegerVector = ispositiveintegervector (x)
% Explanation:
%       Tests whether the input is a positive integer vector
% Example(s):
%       ispositiveintegervector(1:10)
%       ispositiveintegervector(magic(3))
%       ispositiveintegervector(-1:3)
% Outputs:
%       isPositiveIntegerVector
%                       - whether the input is a positive integer vector
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%       cd/isaninteger.m
%
% Used by:
%       cd/compute_peak_halfwidth.m
%       cd/extract_columns.m
%       cd/match_dimensions.m
%       cd/m3ha_select_cells.m

% File History:
% 2018-10-24 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isPositiveIntegerVector = isnumeric(x) && isvector(x) && ...
                            all(isaninteger(x) & x > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%