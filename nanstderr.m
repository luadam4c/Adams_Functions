function stdErr = nanstderr(X, varargin)
%% Calculate the standard error of the mean excluding NaN values
% Usage: stdErr = nanstderr(X, dim)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_stats.m
%		/media/adamX/Paula_IEIs/paula_iei4.m
%       cd/ZG_fit_IEI_distributions.m

% File History:
% 2017-12-14 Created
% 2019-03-14 Added 'dim' as an optional argument

%% Default values for optional arguments
dimDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'X', ...                   % vectors
    @(x) assert(isnum(x), 'X must be a numeric array!'));

% Add optional inputs to the Input Parser
addOptional(iP, 'dim', dimDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Read from the Input Parser
parse(iP, X, varargin{:});
dim = iP.Results.dim;

%% Do the job
if isempty(dim)
    stdErr = nanstd(X, 0) ./ sqrt(length(X(~isnan(X))));
else
    stdErr = nanstd(X, 0, dim) ./ sqrt(size(X(~isnan(X)), dim));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
