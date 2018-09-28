function rmsError = compute_rms_error(vec1, varargin)
%% Computes the root mean squared error given two vectors
% Usage: rmsError = compute_rms_error(vec1, varargin)
% Outputs:
%       rmsError    - root mean squared error
%                   specified as a numeric scalar
% Arguments:    
%       vec1        - the first vector
%                   must be a numeric vector
%       vec2        - (opt) the second vector
%                   must be a numeric vector with length equal to vec1
%                   default == nanmean(vec1) * ones(size(vec1))
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/import_rawtraces.m
%       /media/adamX/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m
%
% File History:
% 2018-07-09 Modified from the built-in rms.m
%   TODO: implement dim as in rms.m
% 

%% Default values for optional arguments
vec2Default = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vec1', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'vec2', [], ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, vec1, varargin{:});
vec2 = iP.Results.vec2;

% Set dependent argument defaults
if isempty(vec2)
    % Compare each sample point to the mean of the entire vector
    vec2 = nanmean(vec1) * ones(size(vec1));
end

% Check relationships between arguments
if length(vec1) ~= length(vec2)
    error('The two vectors must have the same length!');
end

%% Perform job
% Compute errors at every sample point
errors = vec1 - vec2;

% Compute the squared error
squaredError = (errors) .* conj(errors);

% Compute the mean-squared error
meanSquaredError = nanmean(squaredError);

% Compute the root-mean-squared error
rmsError = sqrt(meanSquaredError);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
