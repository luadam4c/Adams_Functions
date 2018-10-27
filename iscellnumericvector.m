function isCellNumericVector = iscellnumericvector (x)
%% Returns whether an input is a cell array of numeric vectors
% Usage: isCellNumericVector = iscellnumericvector (x)
% Explanation:
%       Tests whether the input is a cell array of numeric vectors
% Example(s):
%       iscellnumericvector({1:10, 2:20})
%       iscellnumericvector({magic(3), 2:20})
%       iscellnumericvector({'sets', 'lasts'})
% Outputs:
%       isCellNumericVector   
%                       - whether the input is a cell array of numeric vectors
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/find_window_endpoints.m

% File History:
% 2018-10-25 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

%% Do the job
isCellNumericVector = iscell(x) && ...
                        all(cellfun(@(x) isnumeric(x) && isvector(x), x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%