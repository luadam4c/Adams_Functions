function varargout = spike2Mat2Text (spike2MatPath, varargin)
%% Converts a Spike2-exported .mat file to a text file (atf, txt or csv)
% Usage: [textPath, dataStruct] = spike2Mat2Text (spike2MatPath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       textPath     - path to output text file
%                   specified as a character vector
%       dataStruct   - a structure with fields:
%                       spike2FileContents
%                       spike2MatPath
%                       textPath
%                       channelMatrix
%                       channelNames
%                       channelUnits
%                       siSeconds
%                       timeStart
%                       comment
%                   specified as a scalar structure
%
% Arguments:
%       spike2MatPath   - path to Spike2-exported .mat file
%                       must be a string scalar or a character vector
%       varargin    - 'TextType': type of text file
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'atf'   - Axon Text File
%                       'txt'   - Plain Text File
%                       'csv'   - comma-separated value file
%                   default == 'atf'
%                   - 'TextPath': path to text file
%                   must be a string scalar or a character vector
%                   default == replace(spike2MatPath, '.mat', ['.', textType])
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_fields.m
%       cd/atfwrite.m
%       cd/convert_to_char.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/dlmwrite_with_header.m
%       cd/extract_fields.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%
% Used by:
%       /home/Matlab/plethRO1/spike2Loader.m

% File History:
% 2019-09-06 Moved from spike2Loader.m
% 2019-09-06 Added 'TextType' and 'TextPath' as optional arguments
% TODO: Use parse_spike2_mat.m
% 

%% Hard-coded parameters
validTextTypes = {'atf', 'txt', 'csv'};
isTrace = @(x) isfield(x, 'values') && isfield(x, 'interval');

% TODO: Make optional argument
timeStart = 0;

%% Default values for optional arguments
textTypeDefault  = 'atf';
textPathDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'spike2MatPath');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TextType', textTypeDefault, ...
    @(x) any(validatestring(x, validTextTypes)));
addParameter(iP, 'TextPath', textPathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, spike2MatPath, varargin{:});
textType = validatestring(iP.Results.TextType, validTextTypes);
textPath = iP.Results.TextPath;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Create an output Axon Text File path
if isempty(textPath)
    textPath = replace(spike2MatPath, '.mat', ['.', textType]);
end

%% Load the data
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
    varargout{1} = '';
    varargout{2} = struct;
    return
else
    % Force as a matrix
    channelMatrix = force_matrix(channelValues);
end

% Make sure sampling intervals are the same
if numel(unique(siSecondsAll)) > 1
    % Need to resample the vectors
    % TODO
    disp('Not implemented yet!');
    varargout{1} = '';
    varargout{2} = struct;
    return
else
    % Just choose one sampling interval
    siSeconds = siSecondsAll(1);
end

%% Create the text file
% Count the number of samples
nSamples = size(channelMatrix, 1);

switch textType
    case 'atf'
        % Create a comment
        comment = sprintf('Data from Spike2-exported MATLAB file %s', ...
                    spike2MatPath);

        % Write to ATF file(s)
        atfwrite(channelMatrix, 'SignalNames', channelNames, ...
              'SignalUnits', channelUnits, ...
              'SamplingIntervalSeconds', siSeconds, ...
              'TimeStart', timeStart, 'FileName', textPath, ...
              'Comment', comment);
    case 'txt'
        % Create an output text file path
        textPath = replace(spike2MatPath, '.mat', '.txt');

        % Create a time vector in ms
        timeVectorMs = create_time_vectors(nSamples, 'TimeUnits', 'ms', ...
                                    'SamplingIntervalSeconds', siSeconds, ...
                                    'TimeStart', timeStart);

        % Count the number of significant figures needed
        nSigFig = ceil(log10(nSamples));

        % Print data only as an Axon Plain Text File
        dlmwrite_with_header(textPath, [timeVectorMs, channelMatrix], ...
                                'Delimiter', '\t', 'Precision', nSigFig);
    case 'csv'
        % Create an output spreadsheet path
        textPath = replace(spike2MatPath, '.mat', '.csv');

        % Convert to valid variable names
        channelNamesValid = matlab.lang.makeValidName(channelNames);

        % Create a table with each channel as a separate column
        channelTable = table(channelValues{:}, 'VariableNames', ...
                            channelNamesValid);

        % Write the table as a csv file
        writetable(channelTable, textPath);
    otherwise
        error('textType unrecognized!');
end

%% Output info
varargout{1} = textPath;
if nargout >= 2    
    dataStruct.spike2FileContents = spike2FileContents;
    dataStruct.spike2MatPath = spike2MatPath;
    dataStruct.textPath = textPath;
    dataStruct.channelMatrix = channelMatrix;
    dataStruct.channelNames = channelNames;
    dataStruct.channelUnits = channelUnits;
    dataStruct.siSeconds = siSeconds;
    dataStruct.timeStart = timeStart;
    dataStruct.comment = comment;

    varargout{2} = dataStruct;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%