function parsedDataTable = parse_spike2_mat (spike2MatPath, varargin)
%% Parses a Spike2-exported MATLAB file
% Usage: parsedDataTable = parse_spike2_mat (spike2MatPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberR01/Pleth/Data/new_pleth_data/test2AtNight_200Hz.mat';
%       spike2Table = parse_spike2_mat(spike2MatPath);
%       parse_spike2_mat(spike2MatPath, 'ParseGas', true);
%       spike2Table = parse_spike2_mat(spike2MatPath, 'ChannelNames', {'Sound', 'O2', 'Pleth 2', 'WIC#2', 'WIC#1'});
%       [~, matPaths] = all_files('Ext', 'mat');    
%       for i = 1:numel(matPaths); parse_spike2_mat(matPaths{i}); end
%
% Outputs:
%       parsedDataTable     - parsed data
%                           specified as a table
%
% Arguments:
%       spike2MatPath   - path to Spike2-exported .mat file
%                       must be a string scalar or a character vector
%       varargin    - 'ParseText': whether to parse text marks
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ParseGas': whether to parse pleth pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ParseLaser': whether to parse laser pulses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ChannelNames': channel names to restrict to
%                   must be a string array or a cell array of character vectors
%                   default == all channels with a 'values' field and 
%                                   an 'interval' field
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/all_fields.m
%       cd/convert_to_char.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_fields.m
%       cd/extract_subvectors.m
%       cd/find_in_strings.m
%       cd/force_string_end.m
%       cd/match_row_count.m
%       cd/parse_gas_trace.m
%       cd/parse_laser_trace.m
%       cd/plot_traces_spike2_mat.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/spike2Mat2Text.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-09-08 Moved from spike2Mat2Text.m
% 2019-09-11 Added 'ParseGas' as an optional argument
% 2019-09-20 Added 'ParseText' as an optional argument
% 2019-09-29 Added 'ChannelNames' as an optional argument with default {}
% 

%% Hard-coded parameters
MS_PER_S = 1000;
isTrace = @(x) isfield(x, 'values') && isfield(x, 'interval');

%% Default values for optional arguments
parseTextDefault = false;
parseGasDefault = false;
parseLaserDefault = false;
channelNamesDefault = {};          % set later

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
addParameter(iP, 'ParseText', parseTextDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ParseGas', parseGasDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ParseLaser', parseLaserDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ChannelNames', channelNamesDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));

% Read from the Input Parser
parse(iP, spike2MatPath, varargin{:});
parseText = iP.Results.ParseText;
parseGas = iP.Results.ParseGas;
parseLaser = iP.Results.ParseLaser;
channelNamesUser = iP.Results.ChannelNames;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Make sure the MATLAB path ends with .mat
spike2MatPath = force_string_end(spike2MatPath, '.mat');

%% Reorganize the data as a matrix
% Load everything
%   Note: this is a struct of structs
spike2FileContents = load(spike2MatPath);

% Load data
if isempty(channelNamesUser)
    % Find all fields with a 'values' field
    channelStructsTable = all_fields(spike2FileContents, 'ValueFunc', isTrace, ...
                                        'OutputType', 'table');

    % Get all the channel data
    channelStructs = channelStructsTable.Value;

    % Remove the channel structs table
    clear channelStructsTable
else
    % Find all structures
    allChannelStructs = all_fields(spike2FileContents, 'OutputType', 'value');

    % Extract all channel names
    allChannelNames = extract_fields(allChannelStructs, 'title', 'UniformOutput', false);

    % Find the channels of interest
    indOfInterest = cellfun(@(x) find_in_strings(x, allChannelNames, ...
                                    'SearchMode', 'exact', 'MaxNum', 1, ...
                                    'ReturnNaN', true), channelNamesUser);
    
    % Restrict to the channels of interest
    channelStructs = allChannelStructs(indOfInterest);
    
    % Remove the channel structs table
    clear spike2FileContents allChannelStructs
end

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
    if max(nSamples) - min(nSamples) <= 1
        % Just eliminate the last data point
        nSamplesCommon = min(nSamples);
        nSamples = repmat(nSamplesCommon, size(nSamples));
        siSeconds = repmat(mean(siSeconds), size(siSeconds));
        channelValues = extract_subvectors(channelValues, ...
                                        'Endpoints', [1, nSamplesCommon]);
    else
        % Resample to align the data
        % TODO
        disp('Not implemented yet!');
        parsedDataTable = table.empty;
    %     return
    end
end

% Adjust channel start times so that they are all the same
channelStarts = adjust_start_times(channelStarts, siSeconds);

%% Parse text marks it they exist
if parseText
    textTable = parse_text_spike2_mat(spike2MatPath, spike2FileContents);
end

%% Parse gas trace if it exists
% Test if a gas trace exists
isGasTrace = strcmp(channelNames, 'O2');

if any(isGasTrace) && parseGas
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

%% Parse laser trace if it exists
% Test if a laser trace exists
isLaserTrace = strcmp(channelNames, 'Sound');

if any(isLaserTrace) && parseLaser
    % Get laser vector(s)
    laserVec = channelValues{isLaserTrace};

    % Get the sampling interval in ms
    siMs = siSeconds(isLaserTrace) * MS_PER_S;

    % Parse gas vectors and create pulse tables
    parse_laser_trace(laserVec, siMs, 'TraceFileName', spike2MatPath);
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

function textTable = parse_text_spike2_mat (spike2MatPath, spike2FileContents)
% TODO: parse_text_spike2_mat.m

% Hard-coded parameters
textTableSuffix = '_texts';
isText = @(x) isfield(x, 'text');
figSuffix = strcat(textTableSuffix, '_detection');

% Extract file parts
fileBase = extract_fileparts(spike2MatPath, 'pathbase');

% Make sure the sheet base ends with pulseTableSuffix
textTableBase = force_string_end(fileBase, textTableSuffix);
figBase = force_string_end(fileBase, figSuffix);

% Create paths for the pulse table
textTablePath = force_string_end(textTableBase, '.csv');

% Find all fields with a 'text' field
textStructs = all_fields(spike2FileContents, 'ValueFunc', isText, ...
                            'OutputType', 'value');

% Extract just the first structure with the 'text' field
% TODO: What if there are multiple 'text' fields
textStruct = textStructs{1};

% Extract the text content
%   Note: this is a character matrix
textContent = textStruct.text;

% Extract the text times
%   Note: this is a numeric column vector
textTimes = textStruct.times;

% Extract the text codes
%   Note: this is a numeric matrix
% TODO: What is this?
textCodes = textStruct.codes;

% Convert to a cell array
textStrings = cellstr(textContent);


% Record trace path and whether path exists
tracePath = {spike2MatPath};

% Determine if the trace file exists
[tracePath, pathExists] = construct_and_check_fullpath(tracePath);

% Match row counts
[tracePath, pathExists] = ...
    argfun(@(x) match_row_count(x, numel(textStrings)), tracePath, pathExists);

% TODO: Plot detection result

% Place into a table
textTable = table(textStrings, textTimes, textCodes, tracePath, pathExists, ...
                    'VariableNames', {'String', 'Time', 'Code', ...
                                        'tracePath', 'pathExists'});

% Save the table
writetable(textTable, textTablePath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the .mat file
spike2MatFile = matfile(spike2MatPath);

% Extract all structure names
allStructNames = all_fields(spike2MatFile, 'OutputType', 'name');

% Number of channels
nChannels = numel(channelNames);

% Content of matfile
spike2MatContent = whos(spike2MatFile);

% All structure names
structNames = extract_fields(spike2MatContent, 'name');

% Find the actual structure names
[~, actualChannelNames] = find_in_strings(channelNames, structNames);

% Load only data needed
for iChannel = 1:nChannels
    % Get the current channel name
    channelNameThis = channelNames{iChannel};

    % Find the current channel struct
    spike2MatFile
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
