function [tablesAll, sheetFullFileNames] = atf2sheet (varargin)
%% Converts .atf text file(s) to a spreadsheet file(s) (type specified by the 'SheetType' argument)
% Usage: [tablesAll, sheetFullFileNames] = atf2sheet (atfFileOrDir (opt), varargin)
% Explanation:
%       Converts .atf text file(s) to a spreadsheet file(s) 
%           Default file type to convert is csv, 
%               but this can be changed by the 'SheetType' argument)
%
% Example(s):
%       atf2sheet(pwd);
%       atf2sheet(pwd, 'SheetType', 'xlsx');
% Outputs:
%       tablesAll           - data from the atf file(s)
%                           specified as a table or a cell array of tables
%       sheetFullFileNames  - spreadsheet file name(s)
%                           specified as a string scalar or a character vector
%                               or a cell array of them
% Side Effects:
%       Creates a spreadsheet file
% Arguments:    
%       atfFileOrDir    - (opt) .atf file name or directory containing .atf files
%                       must be a string scalar or a character vector
%                       default == pwd
%       varargin    - 'SheetType': sheet type; 
%                       e.g., 'csv', 'xlsx', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   - 'OutFolder': directory to output csv file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of atf file (atfDir)
%                   - 'NLinesToSkip': number of lines to skip
%                   must be a positive integer scalar
%                   default == 2
%                   - 'Delimiter': Delimiter in .atf file
%                   must be a string scalar or a character vector
%                   default == '\t'
%                   - 'Encoding': Encoding in .atf file
%                   must be a string scalar or a character vector
%                   default == 'ISO-8859-15'
%
% Requires:
%       cd/all_files.m
%       cd/check_dir.m
%       cd/issheettype.m
%       cd/parse_file_or_directory.m
%
% Used by:
%       cd/parse_atf_swd.m
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m

% File History:
% 2018-05-16 Created by Adam Lu, some code from abf2mat.m
% 2018-05-17 Changed default OutFolder to atfDir
% 2018-05-23 Now has tablesAll as the first output
% 2018-11-21 Changed default sheetType xlsx -> csv
% 2018-11-29 Now uses parse_file_or_directory.m, all_files.m and check_dir.m
% TODO: isdelimiter.m
% 

%% Hard-coded parameters
N_LINES_TO_SKIP = 2;            % atf files all seem to have two irrelevant lines

%% Default values for optional arguments
atfFileOrDirDefault = pwd;      % convert all .atf files in the present working 
                                %   directory by default
sheetTypeDefault = 'csv';       % save as a comma-separated value file by default
outFolderDefault = '';          % default directory to output spreadsheet file
nLinesToSkipDefault = N_LINES_TO_SKIP;  % default number of lines to skip
delimiterDefault = '\t';        % default delimiter
encodingDefault = 'ISO-8859-15';% default encoding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'atfFileOrDir', atfFileOrDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'NLinesToSkip', nLinesToSkipDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Encoding', encodingDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
atfFileOrDir = iP.Results.atfFileOrDir;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
outFolder = iP.Results.OutFolder;
nLinesToSkip = iP.Results.NLinesToSkip;
delimiter = iP.Results.Delimiter;
encoding = iP.Results.Encoding;

%% Preparation
% Parse first argument
[atfDir, atfFileName, multipleFiles] = parse_file_or_directory(atfFileOrDir);

% Set default output directory
if isempty(outFolder)
    outFolder = atfDir;
end

% Make sure the output directory exists
check_dir(outFolder);

%% Do the job
% Find all .atf files to convert
if multipleFiles
    % List all the .atf files in the directory
    [~, allAtfPaths] = ...
        all_files('Directory', atfDir, 'Extension', '.atf');

    % Separate the file parts
    [allAtfDirs, allAtfBases, allAtfExts] = ...
        cellfun(@fileparts, allAtfPaths, 'UniformOutput', false);

    % Put together the file names
    allAtfFileNames = ...
        cellfun(@(x, y) [x, y], allAtfBases, allAtfExts, 'UniformOutput', false);

    % Convert all files
    [tablesAll, sheetFullFileNames]  = ...
        cellfun(@(x, y) convert_atf2sheet(x, y, outFolder, sheetType, ...
                                        delimiter, nLinesToSkip, encoding), ...
                allAtfDirs, allAtfFileNames, 'UniformOutput', false);
else
    % Convert this file
    [tablesAll, sheetFullFileNames]  = ...
        convert_atf2sheet(atfDir, atfFileName, outFolder, ...
                                sheetType, delimiter, nLinesToSkip, encoding);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [table, sheetFullFileName] = convert_atf2sheet (atfDir, atfFileName, outFolder, sheetType, delimiter, nLinesToSkip, encoding)
%% Convert an .aft file to a spreadsheet file

% Get the full file name of the .atf file
atfFullFileName = fullfile(atfDir, atfFileName);

% Construct the full file name for the sheet file
sheetFullFileName = fullfile(outFolder, ...
                            strrep(atfFileName, '.atf', ['.', sheetType]));

% Read in the .atf file, ignoring the first two lines
table = readtable(atfFullFileName, 'FileType', 'text', ...
                  'Delimiter', delimiter, 'HeaderLines', nLinesToSkip, ...
                  'Encoding', encoding);

% Write to a spreadsheet file
writetable(table, sheetFullFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

sheetTypeDefault = 'xlsx';      % default spreadsheet type

if exist(atfFileOrDir, 'file') == 2                         % it's a file
    % Store the directory containing the file
    [atfDir, atfFileBase, atfFileExt] = fileparts(atfFileOrDir);

    % If the directory is empty, it is the present working directory
    if isempty(atfDir)
        atfDir = pwd;
    end

    % Get the relative path for the file name
    atfFileName = [atfFileBase, atfFileExt];

    % Set flag
    multipleFiles = false;
elseif exist(atfFileOrDir, 'dir') == 7                      % it's a directory
    % The argument is already the full directory
    atfDir = atfFileOrDir;

    % Set atfFileName to empty (will be vary for each file)
    atfFileName = '';

    % Set flag
    multipleFiles = true;
elseif exist(fullfile(pwd, atfFileOrDir), 'file') == 2      % it's a file
    % The first argument is just the file name in current directory
    atfDir = pwd;

    % Get the relative path for the file name
    atfFileName = atfFileOrDir;

    % Set flag
    multipleFiles = false;
elseif exist(fullfile(pwd, atfFileOrDir), 'dir') == 7       % it's a directory
    % The first argument is a subdirectory in current directory
    %   Get full path to directory
    atfDir = fullfile(pwd, atfFileOrDir);

    % Set atfFileName to empty (will be vary for each file)
    atfFileName = '';

    % Set flag
    multipleFiles = true;
else
    message = sprintf('The .atf file or directory %s does not exist!', ...
                        atfFileOrDir);
    mTitle = 'File or Directory Not Found';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                          'MessageMode', 'show', 'Verbose', true, ...
                          'CreateMode', 'replace');
    return;
end

if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

allAtfFiles = dir(fullfile(atfDir, '*.atf'));

% Compute the number of .atf files in the directory
nAtfFiles = length(allAtfFiles);

% Loop through all files
tablesAll = cell(nAtfFiles, 1);
sheetFullFileNames = cell(nAtfFiles, 1);
parfor iFile = 1:nAtfFiles
    atfFileName = allAtfFiles(iFile).name;
    [tablesAll{iFile}, sheetFullFileNames{iFile}] = ...
        convert_atf2sheet(atfDir, atfFileName, outFolder, ...
                            sheetType, delimiter, nLinesToSkip, encoding);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%