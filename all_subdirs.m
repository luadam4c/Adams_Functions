function [subDirs, fullPaths] = all_subdirs(varargin)
%% Returns all the subdirectories in a given directory
% Usage: [subDirs, fullPaths] = all_subdirs(varargin)
%
% Outputs:
%       subDirs     - a files structure for the subdirectories
%                   specified as a structure array with fields:
%                       name
%                       folder
%                       date
%                       bytes
%                       isdir
%                       datenum
%       fullPaths   - the full paths to the subdirectories
%                   specified as a column cell array of character vectors
% Arguments:    
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/construct_and_check_fullpath.m
%       cd/extract_fullpaths.m
%
% Used by: 
%       /media/ashleyX/Recordings/analyze_recordings.m TODO: Update this
%       /home/zhongxiao/SCIDmiceLTP/Code/analyze_SCIDmiceLTP.m
%       

% File History:
% 2018-09-27 Adapted from the web by Adam Lu
%               https://www.mathworks.com/matlabcentral/answers/
%                   166629-is-there-any-way-to-list-all-folders-only-in-the-level-directly-below-a-selected-directory
% 2018-10-04 Renamed subdirs() -> all_subdirs()

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
directoryDefault = '';              % construct_and_check_fullpath('') == pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;

% Make sure the directory is an existing full path
[directory, dirExists] = construct_and_check_fullpath(directory);
if ~dirExists
    subDirs = [];
    fullPaths = {};
    return
end

%% Find subdirectories
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

%% Print to standard output
% Count the number of subdirectories
nSubDirs = numel(subDirs);

% Print appropriate message
if nSubDirs == 0
    fprintf('No subdirectories found in %s!!\n', directory);
elseif verbose
    fprintf('%d subdirectories found in %s!\n', nSubDirs, directory);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

directoryDefault = pwd;             % look for subdirectories
                                    %   the present working directory by default

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%