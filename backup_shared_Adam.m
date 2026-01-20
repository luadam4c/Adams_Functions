% backup_shared_Adam.m
% A script to backup a specific list of folders from a source to a destination.
% Requires: 
%       cd/backup_folders.m

%% Hard-coded parameters
% Define the parent directory where your source folders are located
sourceParentDir = '\\moorelaboratory.dts.usc.edu\Shared\Adam';

% Define the parent directory where you want the backups to go
destParentDir = '\\moorelaboratory.dts.usc.edu\User Directories\al_627\Shared-Adam-Backup';

% List of specific folder names (as strings) you want to back up
foldersToBackup = { ...
    'Meeting Notes', ...
    'Miscellaneous', ...
    'Onboarding', ...
    'scAAV', ...
    'Suckling', ...
    'Whisking', ...
};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
backupStatus = backup_folders(sourceParentDir, destParentDir, foldersToBackup);

if ~backupStatus
    fprintf('\nBackup completed with errors.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%