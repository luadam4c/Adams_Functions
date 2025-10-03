function [phaseReset, phaseChange, relativeResetTime, T0, T1, eventTimeBefore, eventTimeBeforeBefore] = ...
    compute_phase_response (eventTimes, resetTime, toCompute)
%% Compute phase response from event times and a reset time
% Usage: [phaseReset, phaseChange, relativeResetTime, T0, T1, eventTimeBefore, eventTimeBeforeBefore] = ...
%   compute_phase_response (eventTimes, resetTime, toCompute)
% Explanation:
%       This function computes the phase reset and phase change given a series of
%       rhythmic event times and a single reset event time.
%
% Outputs:
%       phaseReset      - The phase of the reset time within the cycle.
%                       specified as a numeric scalar
%       phaseChange     - The change in phase caused by the reset.
%                       specified as a numeric scalar
%       relativeResetTime   - The time of the reset relative to the preceding event.
%                       specified as a numeric scalar
%       T0              - The unperturbed inter-event interval.
%                       specified as a numeric scalar
%       T1              - The perturbed inter-event interval.
%                       specified as a numeric scalar
%       eventTimeBefore - The time of the event immediately before the reset.
%                       specified as a numeric scalar
%       eventTimeBeforeBefore - The time of the event two cycles before the reset.
%                       specified as a numeric scalar
%
% Arguments:
%       eventTimes      - A vector of event times (e.g., whisk peaks).
%                       must be a numeric vector
%       resetTime       - The time of the resetting event (e.g., breath onset).
%                       must be a numeric scalar
%       toCompute       - A logical flag to indicate if computation is needed.
%                       must be a logical scalar
%
% Requires:
%       None
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m
%       \Shared\Code\vIRt\virt_analyze_whisk.m

% File History:
% 2025-09-24 Pulled code from virt_analyze_sniff_whisk.m
% 2025-10-02 Now returns empty if reset time is empty

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% No optional arguments to deal with

% Initialize outputs
if isempty(resetTime)
    % Return empty if resetTime is empty
    warning('Cannot compute phase response if reset time is empty!');
    eventTimeBefore = [];
    eventTimeBeforeBefore = [];
    T0 = [];
    T1 = [];
    relativeResetTime = [];
    phaseReset = [];
    phaseChange = [];
    return
else
    % Return NaN if resetTime is not empty
    eventTimeBefore = NaN;
    eventTimeBeforeBefore = NaN;
    T0 = NaN;
    T1 = NaN;
    relativeResetTime = NaN;
    phaseReset = NaN;
    phaseChange = NaN;

    % Display warning
    if numel(eventTimes) < 3
        warning('Not enough events to compute phase response. Will return NaN');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
% Only compute if desired and if there are enough events
if toCompute && numel(eventTimes) >= 3
    % Find indices of all events before the reset time
    eventNumsBefore = find(eventTimes < resetTime);

    % Find index of the event right after the reset time
    eventNumAfter = find(eventTimes >= resetTime, 1, 'first');

    % We need at least two events before and one event after to compute the phase response
    if numel(eventNumsBefore) >= 2 && ~isempty(eventNumAfter)
        % Get the event times of the last two events before reset
        eventTimeBefore = eventTimes(eventNumsBefore(end));
        eventTimeBeforeBefore = eventTimes(eventNumsBefore(end-1));

        % Get the event time of the first events after reset
        eventTimeAfter = eventTimes(eventNumAfter);        

        % Compute the first inter-event-interval (T0)
        % This is the unperturbed cycle duration
        T0 = eventTimeBefore - eventTimeBeforeBefore;

        % Compute the second inter-event-interval (T1)
        % This is the perturbed cycle duration
        T1 = eventTimeAfter - eventTimeBefore;

        % Compute the difference between reset time and event before
        % This is how far into the cycle the reset occurred
        relativeResetTime = resetTime - eventTimeBefore;

        % Compute the phase of reset within the unperturbed cycle (phase reset)
        % The phase is normalized to be between 0 and 2*pi
        phaseReset = 2*pi * mod(relativeResetTime / T0, 1);

        % Compute the change in event phase associated with the reset (delta phase)
        % This measures the change in period relative to the original period
        phaseChange = 2*pi * (T1 - T0) / T0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%