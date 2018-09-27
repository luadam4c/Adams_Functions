function [subDirs, fullPaths] = subdirs(directory)
%% Get all the subdirectories from a given directory
% Usage: [subDirs, fullPaths] = subdirs(directory)
%
% Requires:
%       cd/extract_fullpaths.m
%
% Used by: 
%       /media/ashleyX/Recordings/analyze_recordings.m
%       

% File History:
% 2018-09-27 Adapted from the web by Adam Lu
%               https://www.mathworks.com/matlabcentral/answers/166629-is-there-any-way-to-list-all-folders-only-in-the-level-directly-below-a-selected-directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a list of all files and folders in this folder
files = dir(directory);

% Get a logical vector that tells which entries are directories
isDir = [files.isdir];

% Get a logical vector that tells which entries are irrelevant ('.' or '..')
isIrrelevant = cellfun(@(x) any(strcmp(x, {'.', '..'})), {files.name});

% Keep only those that are relevant directories
subDirs = files(isDir & ~isIrrelevant);

% Extract the full paths
fullPaths = extract_fullpaths(subDirs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%