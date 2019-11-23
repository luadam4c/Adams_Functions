function isPositiveScalar = ispositivescalar (x)
%% Returns whether an input is a positive scalar
% Usage: isPositiveScalar = ispositivescalar (x)
% Explanation:
%       Tests whether the input is a positive scalar
% Example(s):
%       ispositivescalar(1:10)
%       ispositivescalar(magic(3))
%       ispositivescalar(-1:3)
% Outputs:
%       isPositiveScalar
%                       - whether the input is a positive scalar
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Requires: 
%
% Used by:
%       cd/compute_and_plot_average_response.m
%       cd/parse_multiunit.m
%       cd/plot_traces.m

% File History:
% 2018-12-18 Modified from ispositivevector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isPositiveScalar = isnumeric(x) && isscalar(x) && x > 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%