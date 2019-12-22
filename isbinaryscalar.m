function isBinaryScalar = isbinaryscalar (x)
%% Returns whether an input is a binary scalar (may be empty)
% Usage: isBinaryScalar = isbinaryscalar (x)
% Explanation:
%       Tests whether the input is a binary scalar
% Example(s):
%       isbinaryscalar([])
%       isbinaryscalar(false)
%       isbinaryscalar(2)
%       isbinaryscalar(1)
%       isbinaryscalar([0 0 1])
% Outputs:
%       isBinaryScalar  - whether the input is a binary scalar (may be empty)
%                       specified as a logical scalar
% Arguments:    
%       x               - an input to check
%
% Required:
%       cd/isbinaryarray.m
%
% Used by:
%       cd/plot_fitted_traces.m

% File History:
% 2018-10-31 Modified from isnumericvector.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% TODO: Place in own function
isbinary = @(x) islogical(x) || isnumeric(x) && all(x == 0 | x == 1);

isBinaryScalar = isempty(x) || isscalar(x) && isbinary(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
