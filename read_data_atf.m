function [dataTable, atfParams] = read_data_atf (varargin)
%% Reads the data table from an Axon Text File formatted text file (.atf)
% Usage: [dataTable, atfParams] = read_data_atf (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       filePaths = write_data_atf(myRandomSignals1);
%       dataTable = read_data_atf('FilePaths', filePaths);
%
% Outputs:
%       dataTable  - data table where each column is a vector
%                   specified as a 2D table or a cell array of them
%       atfParams    - structure of metadata with fields:
%                       samplingIntervalSeconds
%                       signalNames
%                       signalUnits
%                       timeStart
%                       comment
%                   specified as a structure array
%
% Arguments:
%       varargin    - 'Recursive': whether to search recursively
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FilePaths': names of '.atf' text files to load
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from pwd
%
% Requires:
%       cd/all_files.m
%       cd/apply_over_cells.m
%       cd/compute_sampling_interval.m
%       cd/construct_and_check_fullpath.m
%       cd/convert_units.m
%       cd/find_first_match.m
%       cd/read_lines_from_file.m
%       cd/sscanf_full.m
%
% Used by:
%       cd/combine_swd_resp_data.m
%       cd/parse_atf_swd.m

% File History:
% 2020-06-28 Modified from write_data_atf.m
% 2020-08-16 Added timeEndSec
% TODO: Don't read ATF scored files
% TODO: Utilize pieceStr to combine data across files as an option

%% Hard-coded parameters
% The following must be consistent with write_data_atf.m
pieceStr = '_piece';            % string in file names that separate pieces

%% Default values for optional arguments
recursiveDefault = false;       % don't search recursively by default
directoryDefault = '';          % set later
filePathsDefault = {};          % detect from pwd by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Recursive', recursiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FilePaths', filePathsDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, varargin{:});
recursive = iP.Results.Recursive;
directory = iP.Results.Directory;
filePaths = iP.Results.FilePaths;

%% Preparation
% Initialize output
dataTable = [];
atfParams = struct;

% Decide on the files to use
if isempty(filePaths)
    % Decide on the directory if not provided
    if isempty(directory)
        % Use the present working directory
        directory = pwd;
    end

    % Find all '.atf' files in the directory
    [~, filePaths] = all_files('Recursive', recursive, ...
                                'Directory', directory, 'Extension', '.atf');

    % Return usage message if no files found
    if isempty(filePaths)
        fprintf('Type ''help %s'' for usage\n', mfilename);
        return
    end
end

% Check if each path exists
[filePaths, pathExists] = ...
    construct_and_check_fullpath(filePaths, 'Directory', directory);

% Return if not all paths exist
if ~all(pathExists)
    return
end

%% Do the job
% Read in each file
if iscell(filePaths)
    [dataTable, atfParams] = ...
        cellfun(@read_data_atf_helper, filePaths, 'UniformOutput', false);

    % Concatenate the structures
    atfParams = apply_over_cells(@vertcat, atfParams);
else
    [dataTable, atfParams] = read_data_atf_helper(filePaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataTable, atfParams] = read_data_atf_helper(filePath)
%% Reads one .atf file

%% Hard-coded parameters
% Must be consistent with write_data_atf.m
nHeaderLines = 7;
delimiter = '\t';

% Extract the acquisition mode
acquisitionModeLine = read_lines_from_file(filePath, 'LineNumber', 3);
acquisitionMode = extractBefore(extractAfter(acquisitionModeLine, '='), '"');

% Extract the comment
commentLine = read_lines_from_file(filePath, 'LineNumber', 4);
comment = extractBefore(extractAfter(commentLine, '='), '"');

% Extract the start time in ms
timeStartLine = read_lines_from_file(filePath, 'LineNumber', 5);
timeStartMs = sscanf_full(timeStartLine, '%g');

% Convert the start time to seconds
timeStartSec = convert_units(timeStartMs, 'ms', 's');

% Read in the signal names
signalsLine = read_lines_from_file(filePath, 'LineNumber', 7);
results1 = textscan(signalsLine, '%q');
signalNames = results1{1};
signalNames{1} = 'Time';

% Read in the signal units
headerLine = read_lines_from_file(filePath, 'LineNumber', 8);
results2 = textscan(headerLine, '%q');
headerNames = results2{1};
signalUnits = extractBefore(extractAfter(headerNames, '('), ')');

% Read in the data
dataTable = readtable(filePath, 'FileType', 'text', ...
                        'HeaderLines', nHeaderLines);

% Fix the column names
dataTable.Properties.VariableNames = signalNames;

% Extract time column if exists
if any(contains(signalNames, 'Time'))
    [~, timeStr] = find_first_match('Time', signalNames);

    % Extract time column if exists
    timeVec = dataTable.(timeStr);

    % Compute sampling interval in seconds
    siSeconds = compute_sampling_interval(timeVec);

    % Store time end
    timeEndSec = timeVec(end);
else
    siSeconds = NaN;
    timeEndSec = NaN;
end

% Save parameters
atfParams.acquisitionMode = acquisitionMode;
atfParams.comment = comment;
atfParams.timeStartSec = timeStartSec;
atfParams.timeEndSec = timeEndSec;
atfParams.siSeconds = siSeconds;
atfParams.signalNames = signalNames;
atfParams.signalUnits = signalUnits;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
