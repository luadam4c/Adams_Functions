function varargout = parse_gas_trace (vectors, siMs, varargin)
%% Parses gas traces
% Usage: [parsedParams, parsedData] = parse_gas_trace (vectors, siMs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data_Not_Used/code_testing/test2AtNight_200Hz.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       channelValues = spike2Table.channelValues;
%       channelNames = spike2Table.channelNames;
%       gasVec = channelValues{strcmp(channelNames, 'O2')};
%       siMs = spike2Table{strcmp(channelNames, 'O2'), 'siSeconds'} * 1000;
%       [parsedParams, parsedData] = parse_gas_trace(gasVec, siMs, 'PulseDirection', 'downward', 'TraceFileName', spike2MatPath);
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
% 2019-09-09 Created by Adam Lu
% 2019-09-10 Added 'PulseDirection' as an optional argument
% 2019-09-13 Now uses parse_repetitive_pulses.m
% 

%% Hard-coded parameters
% Note: Must be consistent with parse_iox.m and plot_relative_events.m
pulseTableSuffix = '_gas_pulses';
pulseShape = 'first-order';

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
                                    otherArguments);
else
    varargout{1} = parse_repetitive_pulses(vectors, siMs, ...
                                    'PulseTableSuffix', pulseTableSuffix, ...
                                    'PulseShape', pulseShape, ...
                                    otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
