function varargout = apply_to_all_subdirs (myFunction, varargin)
%% Apply the same function (must have 'Directory' as a parameter) to all subdirectories
% Usage: varargout = apply_to_all_subdirs (myFunction, varargin)
% Explanation:
%       TODO
% Example(s):
%       filesAll = apply_to_all_subdirs(@all_files);
%       [~, subDirPaths] = apply_to_all_subdirs(@all_subdirs);
%       TODO: contents = apply_to_all_subdirs(@dir)
% Outputs:
%       varargout   - outputs from the application to each subdirectory
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
directory = construct_and_check_fullpath(directory);

%% Do the job
% List all subdirectories
[~, allSubDirPaths] = all_subdirs('Directory', directory);

% Apply to all subdirectories
if nargout >= 3
    [output1, output2, output3] = cellfun(@(x) myFunction('Directory', x), ...
                                allSubDirPaths, 'UniformOutput', false);
elseif nargout >= 2
    [output1, output2] = cellfun(@(x) myFunction('Directory', x), ...
                                allSubDirPaths, 'UniformOutput', false);
elseif nargout >= 1
    output1 = cellfun(@(x) myFunction('Directory', x), ...
                        allSubDirPaths, 'UniformOutput', false);
else
    cellfun(@(x) myFunction('Directory', x), ...
            allSubDirPaths, 'UniformOutput', false);
end

%% Output results
if nargout >= 1
    varargout{1} = output1;
end
if nargout >= 2
    varargout{2} = output2;
end
if nargout >= 3
    varargout{3} = output3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of subdirectories
nSubDirs = numel(allSubDirPaths);

outputs = cell(nSubDirs, 1);
parfor iSubDir = 1:nSubDirs
    outputs{iSubDir} = myFunction('Directory', directory);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%