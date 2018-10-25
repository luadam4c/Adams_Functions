function [files, fullPaths] = all_files(varargin)
%% Get all the files (but not subdirectories) in a given directory
% Usage: [files, fullPaths] = all_files(varargin)
%
% Outputs:
%       files       - a files structure for the files
%                   specified as a structure array with fields:
%                       name
%                       folder
%                       date
%                       bytes
%                       isdir
%                       datenum
%       fullPaths   - the full paths to the files
%                   specified as a column cell array of character vectors
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Extension': file extension to limit to
%                   must be a string scalar or a character vector
%                   default == no limits
%
% Requires:
%       cd/construct_and_check_fullpath.m
%       cd/extract_fullpaths.m
%
% Used by: 
%       cd/parse_all_abfs.m
%       cd/plot_all_abfs.m
%       

% File History:
% 2018-10-04 Modified from all_subdirs.m

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
directoryDefault = '';              % construct_and_check_fullpath('') == pwd
extensionDefault = '';              % set later

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
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
extension = iP.Results.Extension;

% Make sure the directory is an existing full path
[directory, dirExists] = construct_and_check_fullpath(directory);
if ~dirExists
    subDirs = [];
    fullPaths = {};
    return
end

%% Find files
% Get the regular expression to match
if ~isempty(extension)
    regExp = sprintf('%s$', extension);
else
    regExp = '';
end

% Get a list of all files and folders in this folder
filesOrDirs = dir(directory);

% Get a logical vector that tells which entries are directories
isDir = transpose([filesOrDirs.isdir]);

% Get a logical vector that tells which entries matches the regular expression
if ~isempty(regExp)
    % Get all file or directory names
    names = transpose({filesOrDirs.name});

    % Test whether each matches the regular expression
    isMatch = cellfun(@any, regexpi(names, regExp));
else
    isMatch = true(size(filesOrDirs));
end

% Keep only those that are not directories and 
%   are matches to the regular expression
files = filesOrDirs(~isDir & isMatch);

% Extract the full paths
fullPaths = extract_fullpaths(files);

%% Print to standard output
% Count the number of files
nFiles = numel(files);

% Print appropriate message
if nFiles == 0
    fprintf('No files with pattern %s found in %s!!\n', regExp, directory);
elseif verbose
    fprintf('%d files with pattern %s found in %s!\n', nFiles, regExp, directory);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%