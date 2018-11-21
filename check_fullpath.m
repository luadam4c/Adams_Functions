function pathExists = check_fullpath (fullPath, varargin)
%% Checks whether a path or paths exists and prints message if not
% Usage: pathExists = check_fullpath (fullPath, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       pathExists  - whether a path or paths exists
%                   specified as a column logical array
% Arguments:    
%       fullPath    - the full path(s) to file or directory constructed
%                   must be a character vector 
%                       or a cell array or character vectors
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PathType': the type(s) of the path
%                   must be one of:
%                       'folder'
%                       'file'
%                       or a cell array of them
%                   default == based on whether the path has an extension
%
% Used by:    
%       cd/construct_and_check_fullpath.m
%       cd/plot_traces_abf.m

% File History:
% 2018-10-03 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
pathTypeDefault = '';               % set later

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
addRequired(iP, 'fullPath', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['fullPath must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PathType', pathTypeDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['fullPath must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Read from the Input Parser
parse(iP, fullPath, varargin{:});
verbose = iP.Results.Verbose;
pathType = iP.Results.PathType;

%% Do the job for all paths
% If fullPath is a cell array, make sure pathType is also a cell array
if iscell(fullPath) && ~iscell(pathType)
    pathType = repmat({pathType}, size(fullPath));
end

% Check all paths
if iscell(fullPath)
    pathExists = cellfun(@(x, y) check_fullpath_helper(x, y, verbose), ...
                            fullPath, pathType);
else
    pathExists = check_fullpath_helper(fullPath, pathType, verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pathExists = check_fullpath_helper (fullPath, pathType, verbose)

% Detect the path type if not provided
if isempty(pathType)
    % Get the file extension
    [~, ~, fileExt] = fileparts(fullPath);

    % Decide based on the file extension
    if isempty(fileExt)
        pathType = 'folder';
    else
        pathType = 'file';
    end
end

% Return whether the file or directory exists
switch pathType
case 'folder'
    pathExists = isfolder(fullPath);
case 'file'
    pathExists = isfile(fullPath);
otherwise
    error('The path type %s is unrecognized!!', pathType);
end

% Print message if it doesn't exist
if ~pathExists
    fprintf('The %s %s doesn''t exist!!\n', pathType, fullPath);
elseif verbose
    fprintf('The %s %s already exists!!\n', pathType, fullPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%