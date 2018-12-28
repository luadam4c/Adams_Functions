function timeRel = compute_relative_time(ind, tVec, limits, varargin)
%% Computes time(s) relative to limits from indice(s)
% Usage: timeRel = compute_relative_time(ind, tVec, limits, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       timeRel     - time(s) relative to limits
%                   specified as a TODO
% Arguments:
%       ind         - indice(s) of the time vector
%                   must be a TODO
%       tVec        - time vector
%                   must be a TODO
%       limits      - time limits
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/compute_relative_value.m
%
% Used by:
%       cd/plot_pulse_response_with_stimulus.m

% File History:
% 2018-12-28 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'ind'); %, ...
    % TODO: validation function %);
addRequired(iP, 'tVec'); %, ...
    % TODO: validation function %);
addRequired(iP, 'limits'); %, ...
    % TODO: validation function %);

% Add parameter-tVec pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, ind, tVec, limits, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if isnan(ind)
    % If not a number, return not a number
    timeRel = NaN;
elseif all(isaninteger(ind) & ind > 0)
    % Compute the time relative to time axis start
    timeRel = compute_relative_value(tVec(ind), limits);
else
    error('Need to check ind in input parser!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%