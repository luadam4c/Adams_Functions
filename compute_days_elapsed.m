function [daysElapsed] = compute_days_elapsed (date1, date2, varargin)
%% Computes the number of days elapsed from two dates in yyyyMMdd format
% Usage: [daysElapsed] = compute_days_elapsed (date1, date2, varargin)
% Explanation:
%       Computes the number of days elapsed between date1 and date2
%           and returns as numeric array
%
% Example(s):
%       compute_days_elapsed('20251104', '20251102')
%       compute_days_elapsed('20251104', '20251030')
%       compute_days_elapsed('20210125', '20200125')
%       compute_days_elapsed('20200125', '20190125')
%
% Outputs:
%       daysElapsed - days elapsed
%                   specified as a numeric array
%
% Arguments:
%       date1     - first date
%                   must be a numeric or datetime array
%       date2     - second date
%                   must be a numeric or datetime array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       \Shared\scAAV\analyze_reachr_motion.m

% File History:
% 2025-08-13 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'date1');
addRequired(iP, 'date2');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, date1, date2, varargin{:});

%% Preparation
% Convert inputs to date time array
date1 = datetime(num2str(date1), 'InputFormat', 'yyyyMMdd');
date2 = datetime(num2str(date2), 'InputFormat', 'yyyyMMdd');

%% Do the job
% Compute the number of days elapsed
daysElapsed = days(date2 - date1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%