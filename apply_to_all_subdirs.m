function outputs = apply_to_all_subdirs (myFunction, varargin)
%% Apply the same function to all subdirectories
% Usage: outputs = apply_to_all_subdirs (myFunction, varargin)
% Explanation:
%       TODO
% Example(s):
%       filesAll = apply_to_all_subdirs(@all_files)
%       TODO: contents = apply_to_all_subdirs(dir)
% Outputs:
%       outputs     - outputs from the application to each subdirectory
%                   specified as a TODO
% Arguments:
%       myFunction  - a custom function that takes two equivalent arguments
%                       as normal input
%                       e.g., dir, all_files
%                   must be a function handle
%       varargin    - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/all_subdirs.m
%       cd/create_error_for_nargin.m
%       cd/construct_and_check_fullpath.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-07-22 Created by Adam Lu
% TODO: INCOMPLETE!
% 

%% Hard-coded parameters

%% Default values for optional arguments
directoryDefault = '';          % construct_and_check_fullpath('') == pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, myFunction, varargin{:});
directory = iP.Results.Directory;

%% Preparation
% Make sure the directory is an existing full path
[directory, dirExists] = construct_and_check_fullpath(directory);

%% Do the job
% List all subdirectories
[~, allSubDirPaths] = all_subdirs('Directory', directory);

% Count the number of subdirectories
nSubDirs = numel(allSubDirPaths);

% Apply to all subdirectories
outputs = cell(nSubDirs, 1);
parfor iSubDir = 1:nSubDirs
    outputs{iSubDir} = myFunction('Directory', directory);
end

%% Output results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%