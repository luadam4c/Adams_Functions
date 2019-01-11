function onPath = is_on_path (folder, varargin)
%% Returns whether folder(s) are on the MATLAB path
% Usage: onPath = is_on_path (folder, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       onPath      - whether a folder is on the MATLAB path
%                   specified as a logical vector
% Arguments:
%       folder      - directory(ies) to check if on MATLAB path
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/is_in_array.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-01-10 Copied from https://www.mathworks.com/matlabcentral/answers/
%               86740-how-can-i-determine-if-a-directory-is-on-the-matlab-path-programmatically
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'folder', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['folder must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, folder, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Get all paths in a cell array
%   Note: path returns paths as a single character array delimited by pathsep
%         pathsep is ':' on Linux
pathCell = regexp(path, pathsep, 'split');

% Test whether the folder is in the path
%   Note: Windows is not case-sensitive, but UNIX systems are

tic
if ispc
    if iscell(folder)
        onPath = cellfun(@(x) any(strcmpi(x, pathCell)), folder);
    elseif isstring(folder)
        onPath = arrayfun(@(x) any(strcmpi(x, pathCell)), folder);
    else
        onPath = any(strcmpi(folder, pathCell));
    end
else
    if iscell(folder)
        onPath = cellfun(@(x) any(strcmp(x, pathCell)), folder);
    elseif isstring(folder)
        onPath = arrayfun(@(x) any(strcmp(x, pathCell)), folder);
    else
        onPath = any(strcmp(folder, pathCell));
    end
end
toc

tic
if ispc
    onPath = is_in_array(folder, pathCell, 'IgnoreCase', false);
else
    onPath = is_in_array(folder, pathCell, 'IgnoreCase', true);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% This is faster, but code is cumbersome
if ispc
    if iscell(folder)
        onPath = cellfun(@(x) any(strcmpi(x, pathCell)), folder);
    elseif isstring(folder)
        onPath = arrayfun(@(x) any(strcmpi(x, pathCell)), folder);
    else
        onPath = any(strcmpi(folder, pathCell));
    end
else
    if iscell(folder)
        onPath = cellfun(@(x) any(strcmp(x, pathCell)), folder);
    elseif isstring(folder)
        onPath = arrayfun(@(x) any(strcmp(x, pathCell)), folder);
    else
        onPath = any(strcmp(folder, pathCell));
    end
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%