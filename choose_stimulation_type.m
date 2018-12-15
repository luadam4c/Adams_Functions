function stimType = choose_stimulation_type (responseType)
%% Chooses the stimulation type based on the response type
% Usage: stimType = choose_stimulation_type (responseType)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       stimType    - the channel type for the stimulation response
%                   specified as a one of: 
%                       'Voltage'       - voltage
%                       'Current'       - current
% Arguments:
%       responseType- the channel type for the pulse response
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'Voltage'       - voltage
%                       'Current'       - current
%                       'Conductance'   - conductance
%                       'Other'         - other un-identified types
%
% Used by:
%       cd/compute_all_pulse_responses.m
%       cd/compute_average_pulse_response.m

% File History:
% 2018-12-15 Created by Adam Lu
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

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
addRequired(iP, 'responseType', ...
    @(x) any(validatestring(x, validChannelTypes)));

% Read from the Input Parser
parse(iP, responseType);

%% Do the job
switch responseType
    case 'Voltage'
        stimType = 'Current';
    case {'Current', 'Conductance', 'Other'}
        stimType = 'Voltage';
    otherwise
        error('Code logic error!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%