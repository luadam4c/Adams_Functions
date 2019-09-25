function [swdManualTable, swdManualCsvFile] = ...
                parse_atf_swd (originalEventFile, varargin)
%% Parse spike-wave-discharge (SWD) event info from .atf file or converted .csv file
% Usage: [swdManualTable, swdManualCsvFile] = ...
%               parse_atf_swd (originalEventFile, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       swdManualTable = parse_atf_swd('WAGS04_30_2018_cage3_Manual_SWDs.atf');
%
% Outputs:
%       swdManualTable      - a table of spike-wave discharge event info
%                           specified as a 2D table
%
% Arguments:
%       originalEventFile   - original event file, could be .atf or 
%                               converted .csv
%                           must be a string scalar or a character vector
%       varargin    - 'TraceFileName': Name of the corresponding trace file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == extracted from the .atf file
%                   - 'OutFolder': directory to output swd table file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of originalEventFile
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%
% Requires:
%       cd/argfun.m
%       cd/atf2sheet.m
%       cd/construct_and_check_fullpath.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/is_overlapping.m
%       cd/issheettype.m
%       cd/match_dimensions.m
%       cd/read_lines_from_file.m
%       cd/sscanf_full.m
%
% Used by:
%       cd/parse_all_swds.m
%       cd/plot_traces_EEG.m

% File History:
% 2018-11-21 Created by Adam Lu
% 2018-12-26 Added 'SheetType' as an optional argument
% 2019-09-08 Now uses the trace file name as the basis 
%               for constructing sheet file name
% 2019-09-09 Updated the construction of trace file paths
% 2019-09-24 Added check for overlapping windows
% 

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
varNames = {'startTime', 'endTime', 'duration', 'tracePath', 'pathExists'};

%% Default values for optional arguments
traceFileNameDefault = '';      % set later
outFolderDefault = '';          % set later
sheetTypeDefault = 'csv';       % default spreadsheet type

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
addRequired(iP, 'originalEventFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TraceFileName', traceFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, originalEventFile, varargin{:});
traceFileName = iP.Results.TraceFileName;
outFolder = iP.Results.OutFolder;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Decide what file type the first input is
if regexpi(originalEventFile, '.atf$')
    atfFile = originalEventFile;
    atfCsvFile = '';
elseif regexpi(originalEventFile, '_atf.csv$')
    atfFile = '';
    atfCsvFile = originalEventFile;
else
    atfFile = '';
    atfCsvFile = '';
end

% Get the file directory
fileDir = fileparts(originalEventFile);

% Decide on the output folder
if isempty(outFolder)
    outFolder = fileDir;
end

%% Do the job
% Initialize outputs
swdManualTable = table.empty;
swdManualCsvFile = '';

% Read the table from the file
if ~isfile(atfFile) && ~isfile(atfCsvFile)
    % Do nothing
    return
elseif isfile(atfCsvFile)
    % Display warning if atf file also provided
    if isfile(atfFile)
        fprintf(['Table with be read from the csv file ', ...
                'instead of the atf file!\n']);
    end

    % Read in the SWD manual table from the converted csv file
    atfTable = readtable(atfCsvFile);
elseif isfile(atfFile)
    % Make sure it is not an ATF trace file
    isScoredAtf = is_scored_atf(atfFile);
    if ~isScoredAtf
        fprintf('%s is not a scored ATF file!\n', atfFile);
        return;
    end

    % Read in the SWD manual table and print to a csv file
    [atfTable, atfCsvFile] = atf2sheet(atfFile, 'SheetType', sheetType);
    fprintf('%s created!\n', atfCsvFile);
end

% Make sure there is an event recorded
if height(atfTable) == 0
    return
end

% Get the first channel name
firstSignalName = atfTable{1, 'Signal'};

% Check whether each row is the same as firstSignalName
isFirstSignal = strcmp(atfTable.Signal, firstSignalName);

% Restrict to the entries for the first channel only
swdManualTableOfInterest = atfTable(isFirstSignal, :);

% Get the start and end times in ms
startTimesMs = swdManualTableOfInterest.Time1_ms_;
endTimesMs = swdManualTableOfInterest.Time2_ms_;

% Convert to seconds
startTime = startTimesMs / MS_PER_S;
endTime = endTimesMs / MS_PER_S;

% Make sure none of the windows overlap
isOverlapping = is_overlapping(transpose([startTime, endTime]));
if isOverlapping
    fprintf('The file %s cannot be parsed because windows overlap!\n', ...
            originalEventFile);
    return
end

% Compute duration
duration = endTime - startTime;

% If not provided, read in the trace file names
if isempty(traceFileName)
    % Get the .abf file name for each SWD
    traceFileName = swdManualTableOfInterest.FileName;
end

% Construct full path to original data file
%   Note: Must be (or copied to) the same directory as the original event file
[tracePath, pathExists] = ...
    construct_and_check_fullpath(traceFileName, 'Directory', fileDir);

% Force as a column cell array
tracePath = force_column_cell(tracePath);

% Make sure the dimensions match up
[tracePath, pathExists] = ...
    argfun(@(x) match_dimensions(x, size(startTime)), tracePath, pathExists);

% Extract the file base of the trace file
traceFileBase = extract_fileparts(traceFileName{1}, 'base');

% Extract the file extension of the trace file
traceFileExt = extract_fileparts(traceFileName{1}, 'ext');

% Construct manual SWD table csv file
swdManualCsvFile = ...
    fullfile(outFolder, [traceFileBase, '_Manual_SWDs.', sheetType]);

%% Correct the start and end times if the data comes from an ATF file
if strcmpi(traceFileExt, '.atf')
    % Note: The scored atf file itself also has the info, 
    %   but at a weird location

    % Read the sweep start time line
    headerLine = read_lines_from_file(tracePath{1}, 'MaxNum', 1, ...
                                'Keyword', 'SweepStartTimesMS');
    if isempty(headerLine)
        error('headerLine not found in %s', tracePath{1});
    end
                            
    % Read in the start time of the trace
    traceStartTimeMs = sscanf_full(headerLine, '%g');

    % Extract the start time of the trace
    traceStartTimeSeconds = traceStartTimeMs / MS_PER_S;

    % Modify the start and end times
    startTime = traceStartTimeSeconds + startTime;
    endTime = traceStartTimeSeconds + endTime;
end

%% Output results
% Create a table for the parsed SWDs
swdManualTable = table(startTime, endTime, duration, ...
                        tracePath, pathExists, 'VariableNames', varNames);

% Write the table to a file
writetable(swdManualTable, swdManualCsvFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isScoredAtf = is_scored_atf (atfFile)
%% Determine whether the ATF file is scored

% Read in the second line of the file
line2 = read_lines_from_file(atfFile, 'LineNumber', 2);

% Read in the number of lines to skip following this line before the header
%   is reached
nLinesToSkip = sscanf(line2, '%d', 1);

% This number should be zero for a scored ATF file
isScoredAtf = nLinesToSkip == 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%