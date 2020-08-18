function [nEventsBefore, nEventsAfter] = count_events (eventTimes, varargin)
%% Counts events before and after stimulation
% Usage: [nEventsBefore, nEventsAfter] = count_events (eventTimes, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [nEventsBefore, nEventsAfter] = count_events (-5:5)
%       [nEventsBefore, nEventsAfter] = count_events (-5:5, 'StimTime', 1)
%       [nEventsBefore, nEventsAfter] = count_events ({-5:5; -6:2}, 'StimTime', 1)
%
% Outputs:
%       nEventsBefore   - number of events before stimulation
%                       specified as a numeric array
%       nEventsAfter    - number of events after stimulation
%                       specified as a numeric array
%
% Arguments:
%       eventTimes  - an array or a cell array of event times
%                   must be a TODO
%       varargin    - 'StimTime': stimulation time
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_relative_events.m

% File History:
% 2020-08-18 Moved from plot_relative_events.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
stimTimeDefault = 0;

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
addRequired(iP, 'eventTimes');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'StimTime', stimTimeDefault);

% Read from the Input Parser
parse(iP, eventTimes, varargin{:});
stimTime = iP.Results.StimTime;

%% Do the job
if iscell(eventTimes)
    [nEventsBefore, nEventsAfter] = ...
        cellfun(@(x) count_events_helper(x, stimTime), eventTimes);
else
    [nEventsBefore, nEventsAfter] = count_events_helper(eventTimes, stimTime);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nEventsBefore, nEventsAfter] = ...
                count_events_helper (eventTimes, stimTime)

nEventsBefore = numel(eventTimes(eventTimes < stimTime));
nEventsAfter = numel(eventTimes(eventTimes >= stimTime));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%