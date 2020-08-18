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
%                   must be a numeric array
%                       or a cell array of numeric arrays
%       varargin    - 'StimTime': stimulation time
%                   must be a numeric scalar
%                   default == 0
%                   - 'BeforeWindow': time window before stim
%                   must be empty or a 2-element numeric vector
%                   default == [-Inf, stimTime]
%                   - 'AfterWindow': time window after stim
%                   must be empty or a 2-element numeric vector
%                   default == [stimTime, Inf]

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
beforeWindowDefault = [];
afterWindowDefault = [];

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
addParameter(iP, 'StimTime', stimTimeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'BeforeWindow', beforeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'AfterWindow', afterWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, eventTimes, varargin{:});
stimTime = iP.Results.StimTime;
beforeWindow = iP.Results.BeforeWindow;
afterWindow = iP.Results.AfterWindow;

%% Preparation
% Set default windows
if isempty(beforeWindow)
    beforeWindow = [-Inf, stimTime];
end

if isempty(afterWindow)
    afterWindow = [stimTime, Inf];
end


%% Do the job
if iscell(eventTimes)
    [nEventsBefore, nEventsAfter] = ...
        cellfun(@(x) count_events_helper(x, beforeWindow, afterWindow), ...
                eventTimes);
else
    [nEventsBefore, nEventsAfter] = ...
        count_events_helper(eventTimes, beforeWindow, afterWindow);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nEventsBefore, nEventsAfter] = ...
                count_events_helper (eventTimes, beforeWindow, afterWindow)

nEventsBefore = numel(eventTimes(eventTimes >= beforeWindow(1) & ...
                                    eventTimes < beforeWindow(2)));
nEventsAfter = numel(eventTimes(eventTimes >= afterWindow(1) & ...
                                    eventTimes < afterWindow(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%