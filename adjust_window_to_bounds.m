function winNew = adjust_window_to_bounds (winOrig, bounds, varargin)
%% Adjusts a time window so that it is within specific bounds
% Usage: winNew = adjust_window_to_bounds (winOrig, bounds, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       winNew      - TODO: Description of winNew
%                   specified as a TODO
%
% Arguments:
%       winOrig     - time window to adjust
%                   must be a TODO
%       bounds      - lower and upper bound
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_psth.m
%       cd/plot_psth.m

% File History:
% 2019-09-15 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
% TODO
addRequired(iP, 'winOrig');
addRequired(iP, 'bounds');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, winOrig, bounds, varargin{:});
% param1 = iP.Results.param1;

%% Preparation
% Initialize new window
winNew = winOrig;

%% Do the job
if winOrig(1) < bounds(1)
    winNew(1) = bounds(1);
end
if winOrig(2) > bounds(2)
    winNew(2) = bounds(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%