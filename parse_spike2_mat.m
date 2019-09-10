function parsedDataTable = parse_spike2_mat (spike2MatPath, varargin)
%% Parses a Spike2-exported MATLAB file
% Usage: parsedDataTable = parse_spike2_mat (spike2MatPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data/new_pleth_data/test2AtNight_200Hz.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%
% Outputs:
%       parsedDataTable     - parsed data
%                           specified as a table
%
% Arguments:
%       spike2MatPath   - Spike2-exported MATLAB path
%                       must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_fields.m
%       cd/convert_to_char.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/force_string_end.m
%       cd/parse_gas_trace.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-08 Moved from spike2Mat2Text.m
% 

%% Hard-coded parameters
MS_PER_S = 1000;
isTrace = @(x) isfield(x, 'values') && isfield(x, 'interval');

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'spike2MatPath', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, spike2MatPath, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Make sure the MATLAB path ends with .mat
spike2MatPath = force_string_end(spike2MatPath, '.mat');

%% Do the job
%   Note: this is a struct of structs
spike2FileContents = load(spike2MatPath);

%% Reorganize the data as a matrix
% Find all fields with a 'values' field
channelStructsTable = all_fields(spike2FileContents, 'ValueFunc', isTrace);

% Get all the channel data
channelStructs = channelStructsTable.Value;

% Get all the channel names
channelNames = extract_fields(channelStructs, 'title', 'UniformOutput', false);

% Get all the channel units
channelUnits = extract_fields(channelStructs, 'units', 'UniformOutput', false);
channelUnits = cellfun(@convert_to_char, channelUnits, 'UniformOutput', false);

% Get all the channel values
channelValues = extract_fields(channelStructs, 'values', 'UniformOutput', false);

% Get all the channel scales
channelScales = extract_fields(channelStructs, 'scale', 'UniformOutput', true);

% Get all the channel offsets
channelOffsets = extract_fields(channelStructs, 'offset', 'UniformOutput', true);

% Get all the channel starting times
channelStarts = extract_fields(channelStructs, 'start', 'UniformOutput', true);

% Get all the sampling intervals
siSeconds = extract_fields(channelStructs, 'interval', 'UniformOutput', true);

% Count all the number of samples
nSamples = count_samples(channelValues);

% Compute the total durations in seconds
totalDurationsSec = siSeconds .* nSamples;

% Make sure all data vectors are the same length
if numel(unique(nSamples)) > 1 || numel(unique(siSeconds)) > 1
    % Resample to align the data
    % TODO
    disp('Not implemented yet!');
    parsedDataTable = table.empty;
    return
end

% Adjust channel start times so that they are all the same
channelStarts = adjust_start_times(channelStarts, siSeconds);

%% Parse gas trace if it exists
% Test if a gas trace exists
isGasTrace = strcmp(channelNames, 'O2');

if any(isGasTrace)
    % Get gas vector(s)
    gasVec = channelValues{isGasTrace};

    % Get the sampling interval in ms
    siMs = siSeconds(isGasTrace) * MS_PER_S;

    % TODO: Change this one CO2 is read as well
    pulseDirection = 'downward';

    % Parse gas vectors and create pulse tables
    parse_gas_trace(gasVec, siMs, 'TraceFileName', spike2MatPath, ...
                    'PulseDirection', pulseDirection);
end

%% Output results
% Place in a table
parsedDataTable = ...
    table(channelNames, channelUnits, channelValues, channelScales, ...
            channelOffsets, channelStarts, ...
            siSeconds, nSamples, totalDurationsSec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function channelStarts = adjust_start_times (channelStarts, siSeconds)

% Make sure it is a multiple of siSeconds
channelStarts = floor(channelStarts ./ siSeconds) .* siSeconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%