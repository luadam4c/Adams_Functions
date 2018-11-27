function [swdAssystTable, swdAssystCsvFile] = ...
                parse_assyst_swd (assystFile, varargin)
%% Parse spike-wave-discharge (SWD) event info from an Assyst.txt file
% Usage: [swdAssystTable, swdAssystCsvFile] = ...
%               parse_assyst_swd (assystFile, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       assystFile  - Assyst output file, must have the ending 'Assyst.txt'
%                   must be a string scalar or a character vector
%       varargin    - 'AbfFileName': Name of the corresponding .abf file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == [fileBase, '.abf']
%                   - 'OutFolder': directory to output swd table file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of originalEventFile
%
% Requires:
%       cd/argfun.m
%       cd/construct_and_check_abfpath.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/plot_traces_EEG.m

% File History:
% 2018-11-22 Created by Adam Lu
% 2018-11-27 Now reads durationOrig as a duration array 
%               rather than a datetime array
% 

%% Hard-coded parameters
varNamesOrig = {'swdNumber', 'startTimeOrig', 'endTimeOrig', ...
                'durationOrig', 'comment'};
varNames = {'startTime', 'endTime', 'duration', 'abfPath', 'pathExists'};
startTimePattern = '^Record Start time';
samplingRatePattern = '^Sampling rate';
timeStepPattern = '^Time step';
swdPattern1 = '^Event Category: SWD';
swdPattern2 = '^Event Category: Seizures';
nEventsPattern = '^Number of events';
headerPattern = '^#';

% Must be consistent with read_swd_sheet.m
dateTimePattern = 'M/d/yyyy HH:mm:ss.SSS';
durationPattern = 'hh:mm:ss.SSS';

%% Default values for optional arguments
abfFileNameDefault = '';        % set later
outFolderDefault = '';          % set later

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
addRequired(iP, 'assystFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AbfFileName', abfFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, assystFile, varargin{:});
abfFileName = iP.Results.AbfFileName;
outFolder = iP.Results.OutFolder;

%% Preparation
% Get the fileDir and fileBaseAssyst
[fileDir, fileBaseAssyst, ~] = fileparts(assystFile);

% Replace '_Assyst' with '' to get fileBase
fileBase = replace(fileBaseAssyst, '_Assyst', '');

% Decide on the output folder
if isempty(outFolder)
    outFolder = fileDir;
end

% Construct Assyst SWD table csv file
swdAssystCsvFile = fullfile(outFolder, [fileBase, '_Assyst_SWDs.csv']);

% Set default .abf file name
if isempty(abfFileName)
    abfFileName = [fileBase, '.abf'];
end

% Construct full path to abf file
[abfPath, pathExists] = ...
    construct_and_check_abfpath(abfFileName, ...
                                'Directory', pwd, 'Verbose', false);

% Try again one level up if path doesn't exist
if ~pathExists
    [abfPath, pathExists] = ...
        construct_and_check_abfpath(abfFileName, ...
                                'Directory', fileparts(pwd), 'Verbose', false);
end

% Construct the date time format specifier for textscan()
dateTimeFormatSpec = ['%{', dateTimePattern, '}D'];
timeFormatSpec = ['%{', durationPattern, '}T'];

%% Do the job
% Open the file
fid = fopen(assystFile);

% Initialize an empty string
thisLine = '';

% Read line(s) excluding newline character until reaching startTimePattern
while isempty(regexpi(thisLine, startTimePattern))
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        fprintf('%s is not found!', startTimePattern);
        return
    end
end

% Store the recording start time
tempCell1 = textscan(thisLine, ['%s ', dateTimeFormatSpec], 'Delimiter', '=');
recordingStartTime = tempCell1{2};

% Read line(s) excluding newline character until reaching samplingRatePattern
while isempty(regexpi(thisLine, samplingRatePattern))
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        fprintf('%s is not found!', samplingRatePattern);
        return
    end
end

% Store the sampling rate in Hz
tempCell2 = textscan(thisLine, '%s %f', 'Delimiter', '=');
samplingRateHz = tempCell2{2};

% Read line(s) excluding newline character until reaching timeStepPattern
while isempty(regexpi(thisLine, timeStepPattern))
    thisLine = fgetl(fid);
end

% Store the sampling rate in seconds
tempCell3 = textscan(thisLine, '%s %f', 'Delimiter', '=');
siSeconds = tempCell3{2};

% Read line(s) excluding newline character until reaching swdPattern
while isempty(regexpi(thisLine, swdPattern1)) && ...
        isempty(regexpi(thisLine, swdPattern2))
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        fprintf('%s or %s is not found!', swdPattern1, swdPattern2);
        return
    end
end

% Read line(s) excluding newline character until reaching nEventsPattern
while isempty(regexpi(thisLine, nEventsPattern))
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        fprintf('%s is not found!', nEventsPattern);
        return
    end
end

% Store the number of SWDs
tempCell4 = textscan(thisLine, '%s %d', 'Delimiter', ':');
nSWDs = tempCell4{2};

% Read in the header
while isempty(regexpi(thisLine, headerPattern))
    thisLine = fgetl(fid);
    if ~ischar(thisLine)
        fprintf('%s is not found!', headerPattern);
        return
    end
end

% Read in the SWDs
swdAssystCell = textscan(fid, ['%d ', dateTimeFormatSpec, ' ', ...
                            dateTimeFormatSpec, ' ', ...
                            timeFormatSpec, ' %s'], nSWDs, ...
                            'Delimiter', '\t');

% Convert to a table
swdAssystTable = table(swdAssystCell{:}, 'VariableNames', varNamesOrig);

% Convert the start times to seconds
startTime = seconds(swdAssystTable.startTimeOrig - recordingStartTime);

% Convert the end times to seconds
endTime = seconds(swdAssystTable.endTimeOrig - recordingStartTime);

% Convert the duration to seconds
duration = seconds(swdAssystTable.durationOrig);

% Close the file
fclose(fid);

% Make sure the dimensions match up
[abfPath, pathExists] = ...
    argfun(@(x) match_dimensions(x, size(startTime)), abfPath, pathExists);

%% Output results
% Create analyzed variables
swdAssystTable = addvars(swdAssystTable, ...
                        startTime, endTime, duration, abfPath, pathExists, ...
                        'Before', 'swdNumber', 'NewVariableNames', varNames);

% Write the table to a file
writetable(swdAssystTable, swdAssystCsvFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

abfFileName = fullfile(fileDir, [fileBase, '.abf']);

timePattern = 'HH:mm:ss.SSS';
timeFormatSpec = ['%{', timePattern, '}D'];

duration = second(swdAssystTable.durationOrig);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%