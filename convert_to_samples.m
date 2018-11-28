function timeInSamples = ...
                convert_to_samples (timeLength, samplingInterval, varargin)
%% Converts time(s) from a time unit to samples based on a sampling interval in the same time unit
% Usage: timeInSamples = ...
%               convert_to_samples (timeLength, samplingInterval, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       timeInSamples - time(s) in samples
%                           specified as a numeric vector
% Arguments:
%       timeLength    - time(s) in time units
%                       must be a numeric vector
%       samplingInterval    - sampling interval(s) in time units
%                           must be a numeric vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/m3ha_import_raw_traces.m

% File History:
% 2018-11-28 Created by Adam Lu
% TODO: Add 'RoundMode' as an optional argument with 'round' as default
%       cf. nearest_odd.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'timeLength', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'samplingInterval', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% % Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, timeLength, samplingInterval, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% TODO: option to use floor and ceil instead
timeInSamples = round(timeLength ./ samplingInterval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%