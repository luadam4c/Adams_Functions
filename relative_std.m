function rsd = relative_std (x)
%% Computes the relative standard deviation (%)
% Usage: rsd = relative_std (x)
% Explanation:
%       100 % * standard deviation / abs(mean), ignoring NaN values
% Example(s):
%       relative_std([29, NaN, 30, 35, 40])
% Outputs:
%       rsd     - relative standard deviation (%) 
%               specified as a nonnegative scalar
% Arguments:
%       x       - data array
%               must be a numeric array
%
% Used by:
%       cd/identify_repetitive_pulses.m

% File History:
% 2018-12-15 Created by Adam Lu
% 

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
addRequired(iP, 'x', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, x);

%% Do the job
rsd = 100 * nanstd(x) / abs(nanmean(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%