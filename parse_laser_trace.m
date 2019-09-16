function varargout = parse_laser_trace (vectors, siMs, varargin)
%% Parses laser traces
% Usage: [parsedParams, parsedData] = parse_laser_trace (vectors, siMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data_Not_Used/code_testing/test2AtNight_200Hz.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       laserVec = channelValues{strcmp(channelNames, 'Sound')};
%       siMs = spike2Table{strcmp(channelNames, 'Sound'), 'siSeconds'} * 1000;
%       [parsedParams, parsedData] = parse_laser_trace(laserVec, siMs, 'TraceFileName', spike2MatPath);
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       vectors     - vectors containing many pulse responses
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vectors
%       siMs        - sampling interval in ms
%                   must be a positive vector
%       varargin    - Any other parameter-value pair for parse_repetitive_pulses()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/parse_repetitive_pulses.m
%
% Used by:
%       cd/parse_spike2_mat.m

% File History:
% 2019-09-12 Modified from parse_gas_trace.m
% 2019-09-13 Now uses parse_repetitive_pulses.m
% 

%% Hard-coded parameters
% Note: Must be consistent with plot_relative_events.m
pulseTableSuffix = '_laser_pulses';
pulseShape = 'square';
pulseDirection = 'upward';
minInterPulseIntervalMs = 10000;

%% Default values for optional arguments

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'vector'}));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, vectors, siMs, varargin{:});

% Keep unmatched arguments for the parse_repetitive_pulses() function
otherArguments = iP.Unmatched;

%% Do the job
% Output variably
if nargout >= 2
    [varargout{1}, varargout{2}] = parse_repetitive_pulses(vectors, siMs, ... 
                            'PulseTableSuffix', pulseTableSuffix, ...
                            'PulseShape', pulseShape, ...
                            'PulseDirection', pulseDirection, ...
                            'MinInterPulseIntervalMs', minInterPulseIntervalMs, ...
                            otherArguments);
else
    varargout{1} = parse_repetitive_pulses(vectors, siMs, ... 
                            'PulseTableSuffix', pulseTableSuffix, ...
                            'PulseShape', pulseShape, ...
                            'PulseDirection', pulseDirection, ...
                            'MinInterPulseIntervalMs', minInterPulseIntervalMs, ...
                            otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
