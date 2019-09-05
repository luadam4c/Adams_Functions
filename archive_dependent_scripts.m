function varargout = archive_dependent_scripts (mFileName, varargin)
%% Archive all dependent scripts of a function
% Usage: functionListTable = archive_dependent_scripts (mFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       archive_dependent_scripts('Glucose_analyze', 'FileExt', 'zip')
%       archive_dependent_scripts('Glucose_analyze', 'FileExt', 'tar')
%       archive_dependent_scripts('Glucose_analyze', 'FileExt', 'gz')
%
% Outputs:
%       functionListTable       - see all_dependent_functions.m
%                               specified as a table
% Arguments:
%       mFileName   - .m file name
%                   must be a string scalar or a character vector
%       varargin    - 'OutFolder': directory to place archive file
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'OutFilePath': full path to archive file
%                   must be a string scalar or a character vector
%                   default == fullfile(outFolder, outFileName)
%                   - 'OutFileName': name of archive file
%                       Note: If provided, OutFilePath will override this
%                   must be a string scalar or a character vector
%                   default == [mFileName, '_dependent_files_', create_time_stamp]
%                   - 'FileExt': file extension for the archive
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'zip'   - Windows zip
%                       'tar'   - tar ball
%                       'gz'    - GNU zip
%                   default == 'zip'
%                   - Any other parameter-value pair for 
%                           all_dependent_functions.m
%
% Requires:
%       cd/all_dependent_functions.m
%       cd/check_dir.m
%       cd/create_time_stamp.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/force_string_end.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/clc2_analyze.m

% File History:
% 2019-08-11 Created by Adam Lu
% 2019-08-16 Added 'OutFilePath' and 'OutFileName' as optional arguments
% 2019-08-20 Now returns appropriate error when .m file cannot be found
% 

%% Hard-coded parameters
validFileExts = {'', 'zip', 'tar', 'gz'};

% TODO: Make these optional arguments
saveListFlag = true;
printListFlag = false;

%% Default values for optional arguments
outFolderDefault = pwd;
outFilePathDefault = '';
outFileNameDefault = '';
fileExtDefault = '';            % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'mFileName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['mFileName must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFilePath', outFilePathDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFileName', outFileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileExt', fileExtDefault, ...
    @(x) any(validatestring(x, validFileExts)));

% Read from the Input Parser
parse(iP, mFileName, varargin{:});
outFolder = iP.Results.OutFolder;
outFilePath = iP.Results.OutFilePath;
outFileName = iP.Results.OutFileName;
fileExt = validatestring(iP.Results.FileExt, validFileExts);

% Keep unmatched arguments for the all_dependent_functions() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% If an empty file name is provided, return error
if exist(mFileName) ~= 2
    mFileName = force_string_end(mFileName, '.m');
    fprintf('The file %s cannot be found!\n', mFileName);
    return
end

% If the file extension for the archive is provided, extract it
%   otherwise, set default
if isempty(fileExt)
    if ~isempty(outFilePath)
        fileExt = extract_fileparts(outFilePath, 'ext');
    elseif ~isempty(outFileName)
        fileExt = extract_fileparts(outFileName, 'ext');
    else
        fileExt = 'zip';
    end
end

% Remove the first '.' from extensions
% TODO: Make this a function remove_string_start.m
%   fileExt = remove_string_start(fileExt, '.');
if regexp(fileExt, '^\.')
    fileExt = extractAfter(fileExt, '.');
end

% Create a file path for the archive
if isempty(outFilePath)
    if ~isempty(outFileName)
        % Extract just the file base
        [~, outFileBase] = fileparts(outFileName);

        % Create a default file path
        outFilePathBase = fullfile(outFolder, outFileBase);
    else
        % Create a default file path
        outFilePathBase = fullfile(outFolder, ...
                        [mFileName, '_dependent_files_', create_time_stamp]);
    end
else
    % Extract just the part without the extension
    outFilePathBase = extract_fileparts(outFilePath, 'pathbase');
end

% Create a temporary folder for copying files
outFolderTemp = outFilePathBase;
check_dir(outFolderTemp);

%% Do the job
% Retrieve a cell array of function paths
fprintf('Retrieving a list of all files dependent on %s ... \n', mFileName);
functionListTable = all_dependent_functions(mFileName, ...
                        'OutFolder', outFolderTemp, ...
                        'OriginalOutput', false, ...
                        'SaveFlag', saveListFlag, ...
                        'PrintFlag', printListFlag, ...
                        otherArguments{:});

% Extract the full paths
pathList = functionListTable.fullPath;

% Count the number of files
nFiles = numel(pathList);

% Copy all of them 
fprintf('Copying files to %s ... \n', outFolderTemp);
parfor iFile = 1:nFiles
    copyfile(pathList{iFile}, outFolderTemp);
end

% Archive the files in the temporary folder
fprintf('Archiving files to %s ... \n', outFilePathBase);
switch fileExt
    case 'zip'
        zip(outFilePathBase, outFolderTemp);
    case 'gz'
        % Create a tar ball
        tar(outFilePathBase, outFolderTemp);

        % Compress the tar ball
        gzip([outFilePathBase, '.tar']);

        % Delete the tar ball
        delete([outFilePathBase, '.tar']);
    case 'tar'
        tar(outFilePathBase, outFolderTemp);
    otherwise
        error('The archive file extension %s is unrecognized!', fileExt);
end

% Remove the temporary folder and all its contents
%   Note: the 's' flag makes it recursive
fprintf('Removing temporary directory %s ... \n', outFolderTemp);
rmdir(outFolderTemp, 's');

%% Output results
varargout{1} = functionListTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
