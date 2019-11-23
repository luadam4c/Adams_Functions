function [fileDir, fileName, multipleFiles] = ...
                parse_file_or_directory (fileOrDir, varargin)
%% Parses whether the argument is a file or a directory
% Usage: [fileDir, fileName, multipleFiles] = ...
%               parse_file_or_directory (fileOrDir, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fileDir     - directory name
%                   specified as a character vector
%       fileName    - file name
%                   specified as a character vector
%       multipleFiles   - whether there are multiple files
%                       specified as a logical scalar
% Arguments:
%       fileOrDir   - file name or directory name
%                   must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/print_or_show_message.m
%
% Used by:
%       cd/atf2sheet.m
%       cd/convert_sheettype.m

% File History:
% 2018-11-29 Moved from atf2sheet.m
% 

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
addRequired(iP, 'fileOrDir', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, fileOrDir, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if isfile(fileOrDir)                            % it's a file
    % Store the directory containing the file
    [fileDir, fileBase, fileExt] = fileparts(fileOrDir);

    % If the directory is empty, it is the present working directory
    if isempty(fileDir)
        fileDir = pwd;
    end

    % Get the relative path for the file name
    fileName = [fileBase, fileExt];

    % Set flag
    multipleFiles = false;
elseif isfolder(fileOrDir)                      % it's a directory
    % The argument is already the full directory
    fileDir = fileOrDir;

    % Set fileName to empty (will be vary for each file)
    fileName = '';

    % Set flag
    multipleFiles = true;
elseif isfile(fullfile(pwd, fileOrDir))         % it's a file
    % The first argument is just the file name in current directory
    fileDir = pwd;

    % Get the relative path for the file name
    fileName = fileOrDir;

    % Set flag
    multipleFiles = false;
elseif isfolder(fullfile(pwd, fileOrDir))       % it's a directory
    % The first argument is a subdirectory in current directory
    %   Get full path to directory
    fileDir = fullfile(pwd, fileOrDir);

    % Set fileName to empty (will be vary for each file)
    fileName = '';

    % Set flag
    multipleFiles = true;
else
    message = sprintf('The file or directory %s does not exist!', ...
                        fileOrDir);
    mTitle = 'File or Directory Not Found';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                          'MessageMode', 'show', 'Verbose', true, ...
                          'CreateMode', 'replace');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%