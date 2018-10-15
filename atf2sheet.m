function [tablesAll, sheetFullFileNames] = atf2sheet (atfFileOrDir, varargin)
%% Converts .atf text file(s) to a spreadsheet file(s) (type specified by the 'SheetType' argument)
% Usage: [tablesAll, sheetFullFileNames] = atf2sheet (atfFileOrDir, varargin)
% Explanation:
%       Converts .atf text file(s) to a spreadsheet file(s) 
%           Default file type to convert is xlsx, 
%               but this can be changed by the 'SheetType' argument)
%
% Example(s):
%       atf2sheet(pwd);
%       atf2sheet(pwd, 'SheetType', 'csv');
% Outputs:
%       tablesAll           - data from the atf file(s)
%                           specified as a table or a cell array of tables
%       sheetFullFileNames  - spreadsheet file name(s)
%                           specified as a string scalar or a character vector
%                               or a cell array of them
% Side Effects:
%       Creates a spreadsheet file
% Arguments:    
%       atfFileOrDir    - .atf file name or directory name
%                       must be a string scalar or a character vector
%       varargin    - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
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
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m

% File History:
% 2018-05-16 Created by Adam Lu, some code from abf2mat.m
% 2018-05-17 Changed default OutFolder to atfDir
% 2018-05-23 Now has tablesAll as the first output
% TODO: isdelimiter.m
% 

%% Hard-coded parameters
N_LINES_TO_SKIP = 2;            % atf files all seem to have two irrelevant lines

%% Default values for optional arguments
sheetTypeDefault = 'xlsx';      % default spreadsheet type
outFolderDefault = '';          % default directory to output spreadsheet file
nLinesToSkipDefault = N_LINES_TO_SKIP;  % default number of lines to skip
delimiterDefault = '\t';        % default delimiter
encodingDefault = 'ISO-8859-15';% default encoding

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
addRequired(iP, 'atfFileOrDir', ...                  % .atf file name
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

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
parse(iP, atfFileOrDir, varargin{:});
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
outFolder = iP.Results.OutFolder;
nLinesToSkip = iP.Results.NLinesToSkip;
delimiter = iP.Results.Delimiter;
encoding = iP.Results.Encoding;

% Parse first argument
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

% Set dependent argument defaults
if isempty(outFolder)
    % Default output directory is atfDir
    outFolder = atfDir;
end

% Make sure the output directory exists
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

% Find all .atf files to convert
if multipleFiles
    % List all the .atf files in the directory
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
else
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%