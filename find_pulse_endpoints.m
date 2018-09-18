function [idxStart, idxEnd] = find_pulse_endpoints (pulse)
%% Finds the indices of a current pulse's start and end
% Usage: [idxStart, idxEnd] = find_pulse_endpoints (pulse)
% Arguments:    
%       pulse       - pulse to find endpoints of
%                   must be a numeric vector
%
% Requires:
%       
% Used by:
%       /home/Matlab/Adams_Functions/find_pulse_response_endpoints.m
%
% File History:
% 2018-07-25 BT - Adapted from find_initial_slopes.m
% 2018-08-10 AL - Change the amplitude to take the value from pulseShifted
%                   rather than from pulse
% 2018-08-10 AL - Now checks number of arguments
% 2018-09-17 AL - Now returns empty indices if there is no pulse

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
addRequired(iP, 'pulse', ...                  % voltage vector
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, pulse);

%% Do the job
% Subtract the trace by the initial value
pulseShifted = pulse - pulse(1);

% Find the maximum absolute value and make that the pulse amplitude
[~, idxAbsMax] = max(abs(pulseShifted));
amplitude = pulseShifted(idxAbsMax);

% Find the sign of the amplitude
signAmplitude = sign(amplitude);

% Change search parameters based on positive or negative pulse
if signAmplitude == 1
    % Find the last point less than 1/4 of the amplitude
    idxStart = find(pulseShifted(1:idxAbsMax) < amplitude * 0.25, 1, 'last');

    % Find the last point greater than 3/4 of the amplitude
    idxEndRel = find(pulseShifted(idxAbsMax:end) > amplitude * 0.75, 1, 'last');
elseif signAmplitude == -1
    % Find the last point greater than 1/4 of the amplitude
    idxStart = find(pulseShifted(1:idxAbsMax) > amplitude * 0.25, 1, 'last');

    % Find the last point less than 3/4 of the amplitude
    idxEndRel = find(pulseShifted(idxAbsMax:end) < amplitude * 0.75, 1, 'last');
else
    idxStart = [];
    idxEnd = [];
    return;
end

% Shift the indices to correspond to entire ivecCpr
idxEnd = idxEndRel + idxAbsMax - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Find the pulse amplitude at that point
amplitude = pulse(idxAbsMax);

%}
