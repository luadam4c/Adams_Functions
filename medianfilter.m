function vecFilt = medianfilter (vec, varargin)
%% Applies a median filter to vectors
% Usage: vecFilt = medianfilter (vec, filtWindow, samplingInterval, varargin)
% Explanation:
%       Same as medfilt1() but with option of using a window in time units
% Example(s):
%       TODO
% Outputs:
%       vecFilt     - filtered vector(s)
%                   specified as a numeric array
% Arguments:
%       vec         - vector(s) to median filter
%                   must be a numeric array
%       filtWindow          - (opt) filtWindow in time units
%                           must be a 2-element numeric vector
%       samplingInterval    - (opt) sampling interval in time units
%                           must be a numeric scalar
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the medfilt1() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/find_nearest_odd.m
%
% Used by:
%       cd/find_passive_params.m

% File History:
% 2019-01-14 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
filtWindowDefault = 3;          % 3 sample points by default
samplingIntervalDefault = 1;    % treat filtWindow as samples by default
% param1Default = [];       % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'filtWindow', filtWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addOptional(iP, 'samplingInterval', samplingIntervalDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vec, varargin{:});
filtWindow = iP.Results.filtWindow;
samplingInterval = iP.Results.samplingInterval;
% param1 = iP.Results.param1;

% Keep unmatched arguments for the medfilt1() function
otherArguments = iP.Unmatched;
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Calculate the median filter window in samples
%   Note: Round down to the nearest odd integer to preserve values!!
%           However, must be >= 1
medianFilterWindowSamples = ...
    find_nearest_odd(filtWindow / samplingInterval, 'Direction', 'down');

% Median filter current vectors
vecsFilt = medfilt1(vecs, medianFilterWindowSamples, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%