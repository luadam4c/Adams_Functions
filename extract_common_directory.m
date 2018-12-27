function parentDir = extract_common_directory (paths, varargin)
%% Extracts the common parent directory of a cell array of file paths
% Usage: parentDir = extract_common_directory (paths, varargin)
% Explanation:
%       TODO
% Example(s):
%       extract_common_directory({'a/b/c', 'a/b/c.m', 'a/b/c/d.m'})
% Outputs:
%       parentDir   - the common parent directory
%                   specified as a character vector
% Arguments:
%       paths       - file paths
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%
% Used by:
%       cd/extract_fileparts.m
%       cd/plot_swd_histogram.m
%       cd/plot_swd_raster.m
%       cd/plot_table.m

% File History:
% 2018-11-27 Created by Adam Lu
% 2018-12-18 Now accepts a character array as the input
% 2018-12-26 Moved code to extract_common_prefix.m


%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'paths', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be a character array or a string array ', ...
        'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paths, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% If empty, just return rempty
if isempty(paths)
    parentDir = '';
    return
end

% Extract only the directory part of each path
directories = extract_fileparts(paths, 'directory');

% Extract the common directory across all directories
if ischar(paths)
    % If a character array, use the directory it is contained in
    parentDir = directories;
else
    % Otherwise, extract the common prefix using filesep as the delimiter
    parentDir = extract_common_prefix(directories, 'Delimiter', filesep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
