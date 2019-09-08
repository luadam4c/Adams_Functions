function varargout = parse_spike2_mat (spike2MatPath, varargin)
%% Parses a Spike2-exported MATLAB file
% Usage: [parsedParamsStruct, parsedDataTable] = ...
%               parse_spike2_mat (spike2MatPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       spike2MatPath = '/media/shareX/2019octoberRO1/Pleth/Data/matFiles/test2AtNight_200Hz.mat';
%       [params, data] = parse_spike2_mat('test2AtNight_200Hz');
%
% Outputs:
%       parsedParamsStruct  - parsed params structure
%                           specified as a structure
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
%       cd/struct2arglist.m
%
% Used by:
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-08 Moved from spike2Mat2Text.m
% 

%% Hard-coded parameters
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

% Get all the channel values
channelValues = extract_fields(channelStructs, 'values', 'UniformOutput', false);

% Get all the channel units
channelUnits = extract_fields(channelStructs, 'units', 'UniformOutput', false);
channelUnits = cellfun(@convert_to_char, channelUnits, 'UniformOutput', false);

% Get all the sampling intervals
siSecondsAll = extract_fields(channelStructs, 'interval', 'UniformOutput', true);

% Count all the number of samples
nSamplesAll = count_samples(channelValues);

% Make sure all data vectors are the same length
if numel(unique(nSamplesAll)) > 1
    % Resample to align the data
    % TODO
    disp('Not implemented yet!');
    varargout{1} = struct.empty;
    varargout{2} = table.empty;
    return
end

% Make sure sampling intervals are the same
if numel(unique(siSecondsAll)) > 1
    % Need to resample the vectors
    % TODO
    disp('Not implemented yet!');
    varargout{1} = struct.empty;
    varargout{2} = table.empty;
    return
else
    % Just choose one sampling interval
    siSeconds = siSecondsAll(1);
end

%% Output results
% Parsed params 
parsedParamsStruct.siSeconds = siSeconds;

% Parsed data
if nargout >= 2
    parsedDataTable = ...
        table(channelNames, channelUnits, channelValues, ...
            siSecondsAll, nSamplesAll, ...
            'VariableNames', {'channelNames', 'channelUnits', ...
                            'channelValues', 'siSeconds', 'nSamples'});
end

% Variable arguments out
varargout{1} = parsedParamsStruct;
if nargout >= 2
    varargout{2} = parsedDataTable;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%