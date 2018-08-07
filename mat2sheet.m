function [sheetFullFileNameAll, tableAll, varsAll, varNamesAll] = mat2sheet (matFileOrDir, varargin)
%% Converts .mat files to a spreadsheet file(s) (type specified by the 'SheetType' argument)
% Usage: [sheetFullFileNameAll, tableAll, varsAll, varNamesAll] = mat2sheet (matFileOrDir, varargin)
% Explanation:
%       Converts .mat file(s) to a spreadsheet file(s) 
%           Default file type to convert is xlsx, 
%               but this can be changed by the 'SheetType' argument)
%
% Example(s):
%       mat2sheet(pwd);
%       mat2sheet(pwd, 'SheetType', 'csv');
% Outputs:
%       sheetFullFileNames  - spreadsheet file name(s)
%                           specified as a string scalar or a character vector
% Side Effects:
%       Creates a spreadsheet file
% Arguments:    
%       matFileOrDir    - .mat file name or directory name
%                       must be a string scalar or a character vector
%       varargin    - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
%                   - 'OutFolder': directory to output csv file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of mat file (matDir)
%
% Requires:
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/print_or_show_message.m
%
% Used by:
%       /home/Matlab/EEG_gui/plot_EEG_event_raster.m
%       /home/Matlab/Adams_Functions/ZG_extract_all_IEIs.m

% File History:
% 2018-05-17 Modified from atf2sheet.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
sheetTypeDefault = 'xlsx';      % default spreadsheet type
outFolderDefault = '';          % default directory to output spreadsheet file

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
addRequired(iP, 'matFileOrDir', ...                  % .mat file name
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, matFileOrDir, varargin{:});
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
outFolder = iP.Results.OutFolder;

% Parse first argument
if exist(matFileOrDir, 'file') == 2                         % it's a file
    % Store the directory containing the file
    [matDir, matFileBase, matFileExt] = fileparts(matFileOrDir);

    % If the directory is empty, it is the present working directory
    if isempty(matDir)
        matDir = pwd;
    end

    % Get the relative path for the file name
    matFileName = [matFileBase, matFileExt];

    % Set flag
    multipleFiles = false;
elseif exist(matFileOrDir, 'dir') == 7                      % it's a directory
    % The argument is already the full directory
    matDir = matFileOrDir;

    % Set flag
    multipleFiles = true;
elseif exist(fullfile(pwd, matFileOrDir), 'file') == 2      % it's a file
    % The first argument is just the file name in current directory
    matDir = pwd;

    % Get the relative path for the file name
    matFileName = matFileOrDir;

    % Set flag
    multipleFiles = false;
elseif exist(fullfile(pwd, matFileOrDir), 'dir') == 7       % it's a directory
    % The first argument is a subdirectory in current directory
    %   Get full path to directory
    matDir = fullfile(pwd, matFileOrDir);

    % Set flag
    multipleFiles = true;
else
    message = sprintf('The .mat file or directory %s does not exist!', ...
                        matFileOrDir);
    mTitle = 'File or Directory Not Found';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                          'MessageMode', 'show', 'Verbose', true, ...
                          'CreateMode', 'replace');
    return;
end

% Set dependent argument defaults
if isempty(outFolder)
    % Default output directory is matDir
    outFolder = matDir;
end

% Make sure the output directory exists
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

% Find all .mat files to convert
if multipleFiles
    % List all the .mat files in the directory
    allMatFiles = dir(fullfile(matDir, '*.mat'));

    % Compute the number of .mat files in the directory
    nMatFiles = length(allMatFiles);

    % Loop through all files
    sheetFullFileNameAll = cell(nMatFiles, 1);
    tableAll = cell(nMatFiles, 1);
    varsAll = cell(nMatFiles, 1);
    varNamesAll = cell(nMatFiles, 1);
    parfor iFile = 1:nMatFiles
        matFileName = allMatFiles(iFile).name;
        [sheetFullFileNameAll{iFile}, tableAll{iFile}, ...
            varsAll{iFile}, varNamesAll{iFile}] = ...
            convert_mat2sheet(matDir, matFileName, outFolder, sheetType);
    end
else
    [sheetFullFileNameAll, tableAll, varsAll, varNamesAll] = ...
        convert_mat2sheet(matDir, matFileName, outFolder, sheetType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sheetFullFileName, table, vars, varNames] = convert_mat2sheet (matDir, matFileName, outFolder, sheetType)
%% Convert an .mat file to a spreadsheet file

% Get the full file name of the .mat file
matFullFileName = fullfile(matDir, matFileName);

% Construct the full file name for the sheet file
sheetFullFileName = fullfile(outFolder, ...
                            strrep(matFileName, '.mat', ['.', sheetType]));

% Load the .mat file into a structure
m = load(matFullFileName);

% Get all variable names from the matfile
%   Note: this is a column cell array
varNames = fieldnames(m);

% Store all variables from the matfile in a cell array
%   Note: this is a column cell array
vars = struct2cell(m);

% Iteratively reduce any structure into its components
varNo = 1;                  % start with the first variable
while varNo <= numel(varNames)
    % Get the current variable name and variable
    varName = varNames{varNo};
    varThis = vars{varNo};

    % If it is a character array, put it in a cell array
    if ischar(varThis)
        varThis = {varThis};
    end

    % If it is a row vector, force it to be a column
    if isrow(varThis)
        varThis = varThis(:);
    end

    % Compute the number of columns
    nCols = size(varThis, 2);

    % Either break the variable apart or continue to the next variable
    if nCols > 1                        % it has more than one column
        [vars, varNames] = ...
            break_variable_apart('ManyColumns', nCols, varNo, varName, ...
                                    varThis, varNames, vars);
    elseif isstruct(varThis)            % it is a structure
        [vars, varNames] = ...
            break_variable_apart('IsStruct', nCols, varNo, varName, ...
                                    varThis, varNames, vars);
    else
        % Make sure the variables are updated if not replaced
        varNames{varNo} = varName;
        vars{varNo} = varThis;

        % Increment the variable number
        varNo = varNo + 1;
    end
end

% Get the final variable count
nVars = numel(vars);

% Find the maximum length of the components of the vars cell array
nRows = max(cellfun(@length, vars));

% Pad the components of the vars cell array to length
%parfor iVar = 1:nVars
for iVar = 1:nVars
    % Get the current variable
    varThis = vars{iVar};

    % Get the variable length
    origLength = length(varThis);

    % If the length is smaller than maximum length, pad 
    if iscellstr(varThis)       % if it is a cell array of strings/chars
        % Create new cell string padded with empty strings
        varNew = cell(nRows, 1);
        for iRow = 1:nRows
            if iRow <= origLength
                varNew{iRow} = varThis{iRow};
            else
                varNew{iRow} = '';
            end
        end

        % Save to vars
        vars{iVar} = varNew;
    elseif isnumeric(varThis)   % if it is numeric
        % Create new cell string padded with empty strings
        varNew = NaN * ones(nRows, 1);
        varNew(1:origLength) = varThis;

        % Save to vars
        vars{iVar} = varNew;
    else
        error('Variable type not recognized for %s', varNames{iVar});
    end
end

% Convert the vars cell array into a table
table = struct2table(cell2struct(vars, varNames, 1));

% Write the table to a spreadsheet file
writetable(table, sheetFullFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vars, varNames] = break_variable_apart(reason, nCols, varNo, varName, varThis, varNames, vars);
% Break the variable apart according to given reason

% Save the other variable names and parts
varNamesBefore = varNames(1:varNo-1);
varsBefore = vars(1:varNo-1);
varNamesAfter = varNames(varNo+1:end);
varsAfter = vars(varNo+1:end);

% Create new variable names and parts
varNamesNew = cell(nCols, 1);
varsNew = cell(nCols, 1);
switch reason
case 'ManyColumns'
    parfor iCol = 1:nCols
        % The new variable names just append the column number
        varNamesNew{iCol} = [varName, '_', num2str(iCol)];

        % The new variables are the corresponding columns
        varsNew{iCol} = varThis(:, iCol);
    end
case 'IsStruct'
    % Get all field names
    fields = fieldnames(varThis);

    % Concatenate the variable name with the field names
    varNamesNew = strcat(varName, '_', fields);

    % Get the new variables
    varsNew = struct2cell(varThis);
otherwise
    error('Reason unrecognized!');
end

% Concatenate the variable names and variables back together
varNames = [varNamesBefore; varNamesNew; varNamesAfter];
vars = [varsBefore; varsNew; varsAfter];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the current variable
varThis = vars.(varName);

%}
