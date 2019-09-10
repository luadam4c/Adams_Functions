function pulseTables = parse_iox (varargin)
%% Parses a .iox.txt file and return a pulse table
% Usage: pulseTables = parse_iox (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       pulseTables = parse_iox;
%       pulseTable = parse_iox('2015_09_28')
%
% Outputs:
%       pulseTable  - TODO: Description of pulseTable
%                   specified as a TODO
%       endTimes    - TODO: Description of pulseTable
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the .iox files, e.g. '2015_09_28'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .iox files to detect
%                   must be a character vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from directory
%                   - 'TraceFileNames': Name of the corresponding trace file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/argfun.m
%       cd/construct_and_check_fullpath.m
%       cd/construct_fullpath.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/find_in_strings.m
%       cd/force_column_cell.m
%       cd/match_dimensions.m
%       cd/read_lines_from_file.m
%       ~/Adams_Functions/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-10 Created by Adam Lu
% 

%% Hard-coded parameters
ioxSuffixExtension = '.rf_1.iox.txt';

% Note: Must be consistent with parse_gas_trace.m
pulseTableSuffix = '_gas_pulses';

%% Default values for optional arguments
directoryDefault = pwd;         % look for .abf files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from directory by default
traceFileNamesDefault = '';      % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'TraceFileNames', traceFileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
traceFileNames = iP.Results.TraceFileNames;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Initialize output
pulseTables = table.empty;

% Decide on the files to use
if isempty(fileNames)
    % Find all iox files in the directory
    [~, fileNames] = all_files('Directory', directory, ...
                                'Extension', ioxSuffixExtension);
else
    % Extract the common parent directory
    directory = extract_common_directory(fileNames);
end

% Return usage message if no .abf files found
if isempty(fileNames)
    fprintf('No .iox files were found in %s!\n', directory);
    return
end

% Force file name(s) with the iox text file ending
fileNames = force_string_end(fileNames, ioxSuffixExtension);

% Extract everything before ioxSuffixExtension
filePathBases = extractBefore(fileNames, ioxSuffixExtension);

% Extract just the file bases
%   Note: After removing the extension, the paths look like directories
fileBases = extract_fileparts(filePathBases, 'dirbase');

% Construct trace paths
if isempty(traceFileNames)
    traceFileNames = fullfile(directory, strcat(fileBases, '.txt'));
else
    traceFileNames = construct_fullpath(traceFileNames, 'Directory', directory);
end

% Make sure the pulse table base ends with pulseTableSuffix
pulseTableNames = force_string_end(fileBases, [pulseTableSuffix, '.csv']);

% Create path(s) to the pulse table file(s)
pulseTablePaths = fullfile(directory, pulseTableNames);

%% Do the job
if iscell(fileNames)
    pulseTables = cellfun(@(x, y, z) parse_one_iox(x, y, z), ...
                            fileNames, traceFileNames, pulseTablePaths, ...
                            'UniformOutput', false);
else
    pulseTables = parse_one_iox(fileNames, traceFileNames, pulseTablePaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pulseTable = parse_one_iox (fileName, traceFileName, pulseTablePath)

%% Hard-coded parameters
eventStr = 'hypoxia';
regexpDateTime = '\d{2}\:\d{2}\:\d{2}\.\d{3}';

%% Read in start and end times
% Read in all lines with the eventStr
lineStrs = read_lines_from_file(fileName, 'Keyword', eventStr);

% Find the lines with start times
[~, startLineStrs] = find_in_strings('period-start', lineStrs, ...
                                    'SearchMode', 'substrings');

% Find the lines with end times
[~, endLineStrs] = find_in_strings('period-stop', lineStrs, ...
                                    'SearchMode', 'substrings');

% Extract the start and end times
[startTimeStrs, endTimeStrs] = ...
    argfun(@(x) regexp(x, regexpDateTime, 'match'), startLineStrs, endLineStrs);

% Convert to the number of seconds
[startTime, endTime] = ...
    argfun(@(x) seconds(cellfun(@duration, x)), ...
            startTimeStrs, endTimeStrs);

%% Check trace paths
% Determine if the trace file exists
[tracePath, pathExists] = construct_and_check_fullpath(traceFileName);

% Force as a cell array
tracePath = force_column_cell(tracePath);

% Match the row count
[tracePath, pathExists] = ...
    argfun(@(x) match_dimensions(x, size(startTime)), tracePath, pathExists);

%% Create the pulse table
% Compute the duration
durationVar = endTime - startTime;

% Create a pulse table
pulseTable = table(startTime, endTime, durationVar, tracePath, pathExists, ...
                    'VariableNames', {'startTime', 'endTime', 'duration', ...
                                        'tracePath', 'pathExists'});

% Write to spreadsheet files
writetable(pulseTable, pulseTablePath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%