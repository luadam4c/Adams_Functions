function success = backup_folders (sourceParentDir, destParentDir, varargin)
%% Backups a specific list of folders (or all) from a source to a destination
% Usage: success = backup_folders (sourceParentDir, destParentDir, foldersToBackup)
% Explanation:
%       This function orchestrates a backup process. It takes a source
%       directory and a destination directory. It can optionally take a
%       specific list of subfolders to backup. If no list is provided, it
%       scans the source directory using all_files.m and backs up all 
%       subdirectories found.
%       It uses all_files.m recursively to find files and only copies 
%       files that are new or have a newer modification timestamp than 
%       the destination version.
%       
%       It generates a timestamped .txt log file of the standard output 
%       in the destination directory.
%
% Example(s):
%       backup_folders('C:\Data', 'D:\Backup');
%       backup_folders('C:\Data', 'D:\Backup', {'Experiment1', 'Notes'});
%
% Outputs:
%       success     - Boolean indicating overall status
%                   returns true if all requested folders processed without error
%
% Side Effects:
%       Creates directories and copies files to destParentDir.
%       Prints status messages to the Command Window.
%       Creates a log file in destParentDir.
%
% Arguments:
%       sourceParentDir - The root directory containing folders to backup
%                   must be a string scalar or character vector
%       destParentDir   - The root directory where backups will be saved
%                   must be a string scalar or character vector
%       varargin    - (opt) foldersToBackup: Cell array of folder names
%                   must be a cell array of character vectors or string array
%                   default == {} (which triggers scanning for all subfolders)
%
% Requires:
%       cd/all_files.m
%       cd/check_dir.m
%
% Used by:
%       cd/backup_shared_Adam.m

% File History:
% 2026-01-19 Created by Gemini
% 2026-01-19 Updated to use all_files.m
% 2026-01-19 Updated to use check_dir.m
% 2026-01-19 Added diary logging and detailed print output

%% Default values for optional arguments
foldersToBackupDefault = {};            % default to empty (implies all)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error('Not enough input arguments. Source and Destination are required.');
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = false;

% Add required inputs
addRequired(iP, 'sourceParentDir', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'destParentDir', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add optional inputs
addOptional(iP, 'foldersToBackup', foldersToBackupDefault, ...
    @(x) assert(iscellstr(x) || isstring(x), ...
    'foldersToBackup must be a cell array of strings or a string array!'));

% Parse inputs
parse(iP, sourceParentDir, destParentDir, varargin{:});
sourceParentDir = iP.Results.sourceParentDir;
destParentDir = iP.Results.destParentDir;
foldersToBackupUser = iP.Results.foldersToBackup;

% Check if source directory exists
if ~exist(sourceParentDir, 'dir')
    error('Source directory does not exist: %s', sourceParentDir);
end

% Handle "Default All" case
% If foldersToBackup is empty, scan sourceParentDir for all directories
if isempty(foldersToBackupUser)
    % Use all_files to find subdirectories
    subDirsStruct = all_files('Directory', sourceParentDir, ...
        'SubDirInstead', true, 'Recursive', false, 'Verbose', false);
    
    if isempty(subDirsStruct)
        fprintf('No subfolders found in %s.\n', sourceParentDir);
        success = true;
        return;
    end

    % Create a cell array of subdirectory names
    foldersToBackup = {subDirsStruct.name};
else
    foldersToBackup = foldersToBackupUser;
end

%% Do the job

% 1. Setup Logging
% Ensure destination exists so we can write the log there
check_dir(destParentDir, 'Verbose', false);

timeStamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
logFileName = sprintf('backup_log_%s.txt', timeStamp);
logPath = fullfile(destParentDir, logFileName);

% Start the diary
if exist(logPath, 'file')
    delete(logPath); 
end
diary(logPath);

% Use onCleanup to ensure diary turns off even if the script crashes
cleanupObj = onCleanup(@() diary('off'));

fprintf('Starting backup process at %s\n', datestr(now));
fprintf('Source: %s\n', sourceParentDir);
fprintf('Destination: %s\n', destParentDir);
if isempty(foldersToBackupUser)
     fprintf('Mode: Auto-scan (all subfolders)\n');
else
     fprintf('Mode: Specific folders list (%d folders)\n', length(foldersToBackup));
end
fprintf('--------------------------------------------------\n');

% 2. Call the main backup logic function
% The run_backup_loop function prints the names of files as they are copied
success = run_backup_loop(sourceParentDir, destParentDir, foldersToBackup);

% 3. Report Final Status
if success
    fprintf('\nBackup completed successfully.\n');
else
    fprintf('\nBackup completed with errors.\n');
end

fprintf('Log saved to: %s\n', logPath);
fprintf('--------------------------------------------------\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function success = run_backup_loop(sourceParent, destParent, folderList)
% Orchestrates the backup for the list of high-level folders using all_files.

success = true;

for i = 1:length(folderList)
    folderName = folderList{i};
    
    % Construct full paths for the specific project/experiment folder
    currentSourcePath = fullfile(sourceParent, folderName);
    currentDestPath = fullfile(destParent, folderName);
    
    % Check if source folder exists
    if exist(currentSourcePath, 'dir')
        try
            fprintf('Processing "%s"...\n', folderName);
            
            % Create the root destination folder for this item if it doesn't exist
            if ~exist(currentDestPath, 'dir')
                fprintf('  [NEW ROOT] %s\n', folderName);
            end
            check_dir(currentDestPath, 'Verbose', false);
            
            % Use all_files recursively to get a flat list of all files
            filesStruct = all_files('Directory', currentSourcePath, ...
                                    'Recursive', true, 'Verbose', false);
            
            if isempty(filesStruct)
                fprintf('  No files found in source.\n');
                continue;
            end
            
            % Iterate through every file found
            for k = 1:length(filesStruct)
                fileInfo = filesStruct(k);
                
                % Determine relative path from the current source root
                % e.g. if currentSourcePath is 'C:\Data\Exp1' 
                % and file is at 'C:\Data\Exp1\Sub\file.txt'
                % relativePath is '\Sub'
                relativePath = strrep(fileInfo.folder, currentSourcePath, '');
                
                % Define destination directory for this specific file
                destFolderPath = fullfile(currentDestPath, relativePath);
                
                % Create destination subdirectory if missing
                check_dir(destFolderPath, 'Verbose', false);
                
                % Define full source and destination filenames
                srcFileFull = fullfile(fileInfo.folder, fileInfo.name);
                destFileFull = fullfile(destFolderPath, fileInfo.name);
                
                % Logic to decide if we copy
                shouldCopy = false;
                statusMsg = '';
                
                if ~exist(destFileFull, 'file')
                    shouldCopy = true;
                    statusMsg = '  [NEW FILE]';
                else
                    % Check if source is newer than destination
                    destData = dir(destFileFull);
                    
                    % Compare datenums (source > dest means source is newer)
                    % precise comparison sometimes needs a small tolerance, 
                    % but strict inequality is usually fine for backups
                    if fileInfo.datenum > destData.datenum
                        shouldCopy = true;
                        statusMsg = '  [UPDATE]';
                    end
                end
                
                % Execute Copy and Print Filename
                if shouldCopy
                    % Print the specific file action to Standard Output (and Diary)
                    fprintf('%s %s\n', statusMsg, fullfile(relativePath, fileInfo.name));
                    
                    % Copy the file
                    [copyStatus, copyMsg] = copyfile(srcFileFull, destFileFull, 'f');
                    if ~copyStatus
                        fprintf('  [ERROR] Could not copy %s: %s\n', fileInfo.name, copyMsg);
                        success = false; 
                    end
                end
            end
            
            fprintf('  Done.\n');
            
        catch ME
            fprintf('[FAILED]\n');
            fprintf('  Error: %s\n', ME.message);
            success = false;
        end
    else
        fprintf('Skipping "%s": Source folder not found.\n', folderName);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%