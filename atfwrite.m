function atfwrite (dataMatrix, varargin)
%% Writes a data matrix to an Axon Text File formatted text file (.atf)
% Usage: atfwrite (dataMatrix, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       atfwrite(myRandomSignals1);
%
% Arguments:
%       dataMatrix  - data matrix where each column is a vector
%                   must be a numeric 2D array
%       varargin    - 'FileName': full or relative path to .atf file
%                   must be a string scalar or a character vector
%                   default == [create_time_stamp, '_', inputname(1), '.atf'];
%                   - 'OutFolder': directory to place .atf file
%                   must be a string scalar or a character vector
%                   default == set in construct_fullpath.m
%                   - 'SamplingIntervalSeconds': sampling interval in seconds
%                   must be a positive vector
%                   default == 1 second
%                   - 'SignalNames': signal names for each signal
%                   must be a string vector or a cell array of character vectors
%                   default == Signal #1, Signal #2, ...
%
% Requires:
%       cd/construct_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_stamp.m
%       cd/create_time_vectors.m
%       cd/force_column_cell.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-03 Created by Adam Lu
% 

%% Hard-coded parameters
nRowsToSkip = 4;

%% Default values for optional arguments
fileNameDefault = '';
outFolderDefault = '';      % set in construct_fullpath.m
signalNamesDefault = {};    % set later
samplingIntervalSecDefault = 0.01;    % sampling at 100 Hz by default

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
addRequired(iP, 'dataMatrix', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileName', fileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    
addParameter(iP, 'SignalNames', signalNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['SignalNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'SamplingIntervalSeconds', samplingIntervalSecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));

% Read from the Input Parser
parse(iP, dataMatrix, varargin{:});
fileName = iP.Results.FileName;
outFolder = iP.Results.OutFolder;
signalNames = iP.Results.SignalNames;
siSeconds = iP.Results.SamplingIntervalSeconds;

%% Preparation
% Count the number of signals
nSignals = size(dataMatrix, 2);

% Count the number of samples
nSamples = size(dataMatrix, 1);

% Create default file name
if isempty(fileName)
    fileName = [create_time_stamp, '_', inputname(1), '.atf'];
end

% Create default file name
filePath = construct_fullpath(fileName, 'Directory', outFolder);

% Decide on signal names
if ~isempty(signalNames)
    % Count the number of signal names
    nNames = numel(signalNames);

    % Make sure they are the same
    if nNames ~= nSignals
        fprintf('Number of signal names and data columns don''t match!\n');
        return
    end
else
    % Create default signal names
    signalNames = create_labels_from_numbers(1:nSignals, 'Prefix', 'Signal #');
end

%% Deal with time
% Create a time vector
timeVectorSec = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                                'SamplingIntervalSeconds', siSeconds);

% Count the number of significant figures needed
nSigFig = ceil(log10(nSamples));

%% Prepare for the ATF header
% Force signal names as a column cell array
signalNames = force_column_cell(signalNames);

% Create a string of signal names separated by commas
signalNamesStr = print_cellstr(signalNames, 'Delimiter', ',', ...
                                'OmitQuotes', true, 'OmitNewline', true, ...
                                'OmitBraces', true, 'ToPrint', false);

% Create the SignalExported line
signalsExportedLine = ['"SignalsExported=', signalNamesStr, '"\n'];

% Create Signals headers
signalsHeader = ['"Signals="'; strcat('"', signalNames, '"')];

% Create the Signals line
signalsLine = print_cellstr(signalsHeader, 'Delimiter', '\t', ...
                            'OmitQuotes', true, 'OmitNewline', false, ...
                            'OmitBraces', true, 'ToPrint', false);

% Print the trace labels with signal units
traceLabels = create_labels_from_numbers(1:nSignals, 'Prefix', 'Trace #', ...
                                        'Suffix', ' ()');

% Create trace headers
traceHeader = ['"Time (s)"'; strcat('"', traceLabels, '"')];

% Create the header info line
headerInfo = print_cellstr(traceHeader, 'Delimiter', '\t', ...
                            'OmitQuotes', true, 'OmitNewline', false, ...
                            'OmitBraces', true, 'ToPrint', false);

%% Create the ATF header
% Open the file
fid = fopen(filePath, 'w');

% Print the version info
fprintf(fid, 'ATF\t1.0\n');

% Print nRowsToSkip and nSignals
fprintf(fid, '%d\t%d\n', nRowsToSkip, nSignals);

% Print the acquisition mode
fprintf(fid, '"AcquisitionMode=Fixed-Length Event-Driven"\n');

% Print sweep start time
fprintf(fid, '"SweepStartTimesMS=0.000"\n');

% Print the SignalsExported line
fprintf(fid, signalsExportedLine);

% Print the Signals line
fprintf(fid, signalsLine);

% Print the header info line
fprintf(fid, headerInfo);

% Close the file
fclose(fid);

%% Append the data
dlmwrite(filePath, [timeVectorSec, dataMatrix], '-append', ...
                    'Delimiter', '\t', 'Precision', nSigFig);

% dlmwrite_with_header(filePath, [timeVectorSec, dataMatrix], ...
%                     'AppendToFile', true, 'Delimiter', '\t', ...
%                     'Precision', nSigFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%