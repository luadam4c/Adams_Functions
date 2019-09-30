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
%                   - 'SignalUnits': signal units for each signal
%                   must be a string vector or a cell array of character vectors
%                   default == '', '', ...
%                   - 'TimeStart': start time in seconds
%                       Note: this is the time a sampling interval before the 
%                               first data time
%                   must be a numeric vector
%                   default == 0 seconds
%                   - 'Comment': comment to put in ATF header
%                   must be a string scalar or a character vector
%                   default == sprintf('Data from %s', inputname(1))
%
% Requires:
%       cd/construct_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_time_stamp.m
%       cd/create_time_vectors.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/print_cellstr.m
%       cd/compute_sigfig.m
%
% Used by:
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-03 Created by Adam Lu
% 2019-09-06 Added 'TimeStart' as an optional argument
% 2019-09-06 Added 'SignalUnits' as an optional argument
% 2019-09-08 Now computes nDecimals instead of nSigFig
% TODO: Add 'TextMarks' somehow
% 

%% Hard-coded parameters
maxNSamplesAtf = 1e6;

% The following must be consistent with combine_swd_sheets.m
pieceStr = '_piece';            % string in file names that separate pieces

%% Default values for optional arguments
fileNameDefault = '';
outFolderDefault = '';              % set in construct_fullpath.m
signalNamesDefault = {};            % set later
signalUnitsDefault = {};            % set later
samplingIntervalSecDefault = 0.01;  % sampling at 100 Hz by default
timeStartDefault = 0;               % start at 0 seconds by default
commentDefault = '';                % set later

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
addParameter(iP, 'SignalUnits', signalUnitsDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['SignalUnits must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'SamplingIntervalSeconds', samplingIntervalSecDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'Comment', commentDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, dataMatrix, varargin{:});
fileName = iP.Results.FileName;
outFolder = iP.Results.OutFolder;
signalNames = iP.Results.SignalNames;
signalUnits = iP.Results.SignalUnits;
siSeconds = iP.Results.SamplingIntervalSeconds;
timeStart = iP.Results.TimeStart;
comment = iP.Results.Comment;

%% Preparation
% Count the number of samples
nSamples = size(dataMatrix, 1);

% Count the number of signals
nSignals = size(dataMatrix, 2);

% Count the number of significant figures needed
%   Note: print at least 5 significant figures
nSigFig = max(5, compute_nsigfig_for_time(nSamples, siSeconds));

% Create a precision string
precision = ['%.', num2str(nSigFig), 'g'];

% Create default file name
if isempty(fileName)
    fileName = [create_time_stamp, '_', inputname(1), '.atf'];
end

% Construct file path to use
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

% Decide on signal units
if ~isempty(signalUnits)
    % Count the number of signal names
    nUnits = numel(signalUnits);

    % Make sure they are the same
    if nUnits ~= nSignals
        fprintf('Number of signal units and data columns don''t match!\n');
        return
    end
else
    % Create default signal units
    signalUnits = repmat({''}, size(signalNames));
end

% Create a default comment
if isempty(comment)
    comment = sprintf('Data from %s', inputname(1));
end

%% Do the job
% Split files if more than maxNSamples samples
if nSamples > maxNSamplesAtf
    % Compute the number of files
    nFiles = ceil(nSamples / maxNSamplesAtf);

    % Extract the path base
    filePathBase = extract_fileparts(filePath, 'pathbase');

    % Create file paths
    filePaths = create_labels_from_numbers(1:nFiles, ...
                'Prefix', [filePathBase, pieceStr], 'Suffix', '.atf');

    for iFile = 1:nFiles
        % TODO: Use create_indices.m somehow
        % Get the starting row
        rowStart = (iFile - 1) * maxNSamplesAtf + 1;

        % Get the ending row
        rowEnd = min(nSamples, iFile * maxNSamplesAtf);

        % Create indices
        rowsThis = rowStart:rowEnd;

        % Compute time start for this file
        timeStartThis = timeStart + (rowStart - 1) * siSeconds;

        % Modify the comment
        commentThis = [comment, ', piece', num2str(iFile)];

        % Create the .atf file
        atfwrite_helper(dataMatrix(rowsThis, :), siSeconds, ...
                        signalNames, signalUnits, ...
                        timeStartThis, commentThis, ...
                        precision, filePaths{iFile});
    end
else
    % Create the .atf file
    atfwrite_helper(dataMatrix, siSeconds, signalNames, signalUnits, ...
                    timeStart, comment, precision, filePath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function atfwrite_helper(dataMatrix, siSeconds, signalNames, signalUnits, ...
                            timeStart, comment, precision, filePath)
%% Writes one .atf file
% TODO: Construct the header text beforehand and pass it in

%% Hard-coded constants
MS_PER_S = 1000;

%% Deal with time
% Count the number of samples
nSamples = size(dataMatrix, 1);

% Create a time vector for this file
timeVectorSec = create_time_vectors(nSamples, 'TimeUnits', 's', ...
                                'SamplingIntervalSeconds', siSeconds, ...
                                'TimeStart', timeStart);
                          
% Compute the time start in milliseconds
timeStartMs = timeStart * MS_PER_S;

%% Construct the ATF header
atfHeaderLines = construct_atf_header(signalNames, signalUnits, ...
                                        timeStartMs, comment);

%% Print the ATF header
% Open the file
fid = fopen(filePath, 'w');
fprintf('Creating file at %s ... \n', filePath);

% TODO: fprintf_custom.m?
% Replace any '%' with '%%'
atfHeaderLines = replace(atfHeaderLines, '%', '%%');

% Print the ATF file header
fprintf(fid, atfHeaderLines);

% Close the file
fclose(fid);

%% Append the data
% Append the data
dlmwrite(filePath, [timeVectorSec, dataMatrix], '-append', ...
                    'Delimiter', '\t', 'Precision', precision);

% dlmwrite_with_header(filePath, [timeVectorSec, dataMatrix], ...
%                     'AppendToFile', true, 'Delimiter', '\t', ...
%                     'Precision', nSigFig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function atfHeaderLines = construct_atf_header (signalNames, signalUnits, ...
                                                timeStartMs, comment)
%% Constructs an Axon Text File header

%% Hard-coded parameters
nRowsToSkip = 5;

%% Preparation
% Count the number of signals
nSignals = numel(signalNames);

%% Do the job
% First line is version info
versionInfoLine = 'ATF\t1.0\n';

% Second line is nRowsToSkip and nSignals
secondLine = sprintf('%d\\t%d\\n', nRowsToSkip, nSignals);

% Acquisition Mode
acquisitionModeLine = '"AcquisitionMode=Fixed-Length Event-Driven"\n';

% Comment Line
commentLine = sprintf('"Comment=%s"\n', comment);

% Sweep Start Time
sweepStartTimesLine = sprintf('"SweepStartTimesMS=%.3f"\n', timeStartMs);

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

% Print the trace labels with no signal units
traceLabels = create_labels_from_numbers(1:nSignals, 'Prefix', 'Trace #', ...
                                        'Suffix', ' ()');

% Replace with correct signal units
traceLabels = cellfun(@(x, y) replace(x, '()', ['(', y, ')']), ...
                        traceLabels, signalUnits, 'UniformOutput', false);

% Create trace headers
traceHeader = ['"Time (s)"'; strcat('"', traceLabels, '"')];

% Create the header info line
headerInfo = print_cellstr(traceHeader, 'Delimiter', '\t', ...
                            'OmitQuotes', true, 'OmitNewline', false, ...
                            'OmitBraces', true, 'ToPrint', false);

% Combine into the ATF file header
atfHeaderLines = [versionInfoLine, secondLine, acquisitionModeLine, ...
                commentLine, sweepStartTimesLine, signalsExportedLine, ...
                signalsLine, headerInfo];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nDecimals = compute_ndecimals (number)
%% Computes the number of decimal places needed
% TODO: Pull out to its own function

% Compute the place 
nDecimalsFirstDigit = -floor(log10(number));

% Compute the number of significant figures
nSigFig = compute_sigfig(number);

if nDecimalsFirstDigit <= 0
    nDecimals = 0;
elseif nSigFig > 1
    nDecimals = nDecimalsFirstDigit + (nSigFig - 1);
else
    nDecimals = nDecimalsFirstDigit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nSigFig = compute_nsigfig_for_time (nSamples, samplingInterval)

% Compute the number of places to the right of the decimal point needed
nPlacesRight = compute_ndecimals(samplingInterval);

% Compute the total duration
totalDuration = nSamples * samplingInterval;

% Compute the number of places to the left of the decimal point needed
nPlacesLeft = ceil(log10(totalDuration));

% Compute the number of significant figures needed
nSigFig = nPlacesRight + nPlacesLeft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
