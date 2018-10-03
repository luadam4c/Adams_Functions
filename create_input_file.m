function create_input_file (dataDir, defParams, varargin)
%% Create an input spreadsheet file from data file names in a directory based on default parameters
% Usage: create_input_file (dataDir, defParams, varargin)
% Explanation:
%       TODO
% Side Effects:
%       Generates an Excel file in the present working directory
% Arguments:    
%       dataDir     - a data directory
%                   must be a directory
%       defParams   - default parameters
%                   must be a structure with fields containing scalar values
%       varargin    - 'DataType': data type
%                   must be a scalartext
%                   default == 'mat'
%                   - 'SheetType': file type for input spreadsheet file
%                   must be a string scalar or a character vector
%                   default == 'xlsx'
%                   - 'OutFolder': directory to output csv file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'show'
%
% Requires:
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:    
%       /home/Matlab/EEG_gui/EEG_gui.m
%
% File History:
% 2018-04-30 Created by Adam Lu
% 2018-05-02 Added outFolder
% 2018-05-02 Changed to using tables and added SheetType as a parameter
% 2018-05-04 Added prompt to overwrite and message boxes

%% Hard-coded parameters
dataFileHeader = 'Files';       % header for the data file names
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
dataTypeDefault = 'mat';        % detect matfiles by default
sheetTypeDefault = 'xlsx';      % default file type for input file
outFolderDefault = '';          % default directory to output csv file
messageModeDefault = 'show';    % default: display message box without pausing

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
addRequired(iP, 'dataDir', ...                  % a data directory
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'defParams', ...                % default parameters
    @(x) validateattributes(x, {'struct'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));

% Read from the Input Parser
parse(iP, dataDir, defParams, varargin{:});
dataType = iP.Results.DataType;
sheetType = iP.Results.SheetType;
outFolder = iP.Results.OutFolder;

% Match possibly ambiguous strings to valid strings
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);

% Set dependent argument defaults
if isempty(outFolder)
    outFolder = pwd;
end

% Make sure the output directory exists
check_dir(outFolder, 'MessageMode', 'show');


%% Create and check full input file name
% Split the data directory path into parts
if isunix
    dataDirParts = strsplit(dataDir, '/');
else
    dataDirParts = strsplit(dataDir, '\');
end

% The "data directory name" is the last part that is:
%   1. nonempty
%   2. does not match matfile
check_criteria = @(x) any(length(x)) && ~any(strfind(x, 'matfile'));
indMatchCriteria = cellfun(@(x) check_criteria(x), dataDirParts);
dataDirName = dataDirParts{find(indMatchCriteria, 1, 'last')};

% Create input spreadsheet file based on data directory name and data type
inputFileName = [dataDirName, '_all', dataType, 'files.', sheetType];

% Create full input file name
inputFullFileName = fullfile(outFolder, inputFileName);

% Check if already exists
if exist(inputFullFileName, 'file') == 2
    qString = {sprintf('An input file of the name %s already exists!'...
                        , inputFullFileName), ...
                'Would you like to replace it?'};
    qTitle = 'Prompt to Deal With Existing Input File';
    choice1 = 'Yes';
    choice2 = 'No';
    answer = questdlg(qString, qTitle, ...
                        choice1, choice2, choice2);
    switch answer
    case choice1
        % Do nothing and continue
    case choice2
        % Exit function
        return;
    otherwise
        % Exit function
        return;        
    end
end

%% Deal with files
% Check if the data directory exists and exit the program
if ~isdir(dataDir)
    msg = sprintf('The data directory %s does not exist!!\n', dataDir);
    mtitle = 'Input File Creation Unsuccessful';
    print_or_show_message(msg, 'Icon', 'error', ...
                         'MTitle', mtitle, 'MessageMode', messageMode);
    return;
end

% Look for files of the given data type in the data directory
%   This returns a row structure array
dataFiles = dir(fullfile(dataDir, ['*.', dataType]));

% Extract the file names as a row cell array, sorted alphabetically
dataFileNames = sort({dataFiles.name});

% Convert the file names to full file names
dataFullFileNames = cellfun(@(x) fullfile(dataDir, x), dataFileNames, ...
                            'UniformOutput', false);

% Display a warning if no files of the given data type are found
if isempty(dataFiles)
    msg = sprintf('Warning: There are no .%s files in %s!!\n', dataType, dataDir);
    mtitle = 'No Data Files';
    print_or_show_message(msg, 'Icon', 'warn', ...
                         'MTitle', mtitle, 'MessageMode', messageMode);
end

% Count the number of files
nFiles = numel(dataFiles);

%% Deal with parameters
% Extract the parameter names
%   This returns a column cell array
paramNames = fieldnames(defParams);

% Count the number of parameters
nParams = numel(paramNames);

% Make sure there is a parameter
if nParams == 0
    msg = sprintf('There are no parameters!!\n');
    mtitle = 'Input File Creation Unsuccessful';
    print_or_show_message(msg, 'Icon', 'error', ...
                         'MTitle', mtitle, 'MessageMode', messageMode);
    return;
end

% Check if all fields in the default parameters structure are scalar values
if ~all(structfun(@isscalar, defParams))
    % Determine which fields are not scalar
    idxfields = find(~structfun(@isscalar, defParams));

    % Print all field names that are not scalars and exit the program
    for idx = idxfields
        msg = sprintf('The parameter %s is not a scalar!!\n', paramNames{idx});
        mtitle = 'Input File Creation Unsuccessful';
        print_or_show_message(msg, 'Icon', 'error', ...
                             'MTitle', mtitle, 'MessageMode', messageMode);
    end
    return;
end

% Repeat the same parameter values for nFiles rows
inputParams = structfun(@(x) x * ones(nFiles, 1), defParams, ...
                        'UniformOutput', false);

%% Create table
% Create a table for the parameter values
paramsTable = struct2table(inputParams);

% Place the file names in a column cell array
fileNames = dataFullFileNames';

% Convert the cell array to a table
namesTable = cell2table(fileNames, 'VariableNames', {dataFileHeader});

% Add the file names to the table at the very left
inputTable = [namesTable, paramsTable];
% inputTable = addVars(paramsTable, fileNames, 'Before', 1);  
                                                % addVars introduced in R2018a

%% Write table
% Write the input table to the spreadsheet file
writetable(inputTable, inputFullFileName);

% Display message box
msg = {sprintf('New input file %s created!', inputFullFileName), 'Well done!'};
mtitle = 'Input File Creation Successful';
print_or_show_message(msg, 'MTitle', mtitle, 'MessageMode', messageMode);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Default parameters
% defParams.WL = 0.5;         % Window size for standard deviation computations 
%                             %     (in seconds)
% defParams.SNR = sqrt(2);    % Signal to noise ratio
% defParams.minD = 2;         % Minimum duration (in seconds)
% defParams.minF = 6.0;       % Minimum frequency (in Hz)
% defParams.maxF = 11;        % Maximum frequency (in Hz)
% defParams.FP = 2;           % Frequency proximity (in Hz) 
%                             %   i.e., how close the harmonics are to each other

% Extract the parameter values
% paramValues = struct2cell(defParams);
paramValues = cellfun(@num2str, struct2cell(defParams), ...
                        'UniformOutput', false);

%% Create cell array
% Allocate a cell array with files as rows and parameters as columns
cellArray = cell(nFiles + 1, nParams + 1);

% Write the header for the data files
cellArray{1, 1} = dataFileHeader;

% Place the parameter names in the cell array as the header
cellArray(1, 2:(nParams + 1)) = paramNames';

% Place the file names in the first column of the cell array
cellArray(2:(nFiles + 1), 1) = dataFullFileNames';

% Place the values in each row of the cell array
cellArray(2:(nFiles + 1), 2:(nParams + 1)) = repmat(paramValues', nFiles, 1);

%% Write cell array
% Split the data directory path into parts
dataDirParts = strsplit(dataDir, '/');

% The "data directory name" is the last part that is:
%   1. nonempty
%   2. does not match matfile
check_criteria = @(x) any(length(x)) && ~any(strfind(x, 'matfile'));
indMatchCriteria = cellfun(@(x) check_criteria(x), dataDirParts);
dataDirName = dataDirParts{find(indMatchCriteria, 1, 'last')};

% Create input Excel file based on data directory name and data type
% inputFileName = [dataDirName, '_all', dataType, 'files.xls'];

% Write cell array to an Excel file in the present working directory
% xlswrite(fullfile(pwd, inputFileName), cellArray);

% Write cell array to a csv file in the output directory
inputFileName = [dataDirName, '_all', dataType, 'files.csv'];
fid = fopen(fullfile(outFolder, inputFileName), 'w');
for iRow = 1:size(cellArray, 1)
    % Print the row
    fprintf(fid, [strjoin(cellArray(iRow, :), ', '), '\n']);
end
fclose(fid);

% Create the header with spaces (doesn't work!)
inputTable.Properties.VariableNames = [dataFileHeader; paramNames];

if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    msg = sprintf('New directory is made: %s\n\n', outFolder);
    mtitle = 'New Directory Made';
    print_or_show_message(msg, 'MTitle', mtitle, 'MessageMode', messageMode);
end

%}
