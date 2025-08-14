function scriptPath = create_new_mscript (scriptName, varargin)
%% Creates a new MATLAB script starting from a function template
% Usage: scriptPath = create_new_mscript (scriptName, (opt) templateNumber, varargin)
% Explanation:
%       Creates a new MATLAB script in the user's functions directory
%           on fishfish or chalkboard, starting from the function templates
%           in ~/Settings_Matlab
%       By default, the template function_template_minimal.m is copied, but
%           this can be varied with the second argument, which is optional
%           The function template numbers are:
%               1 - function_template_minimal.m
%               2 - function_template_simple.m
%               3 - function_template.m
%
% Example(s):
%       create_new_mscript('test_script1');
%       create_new_mscript('test_script2', 2);
%       create_new_mscript('test_script3', 3);
%
% Outputs:
%       scriptPath      - full path to the new MATLAB script
%                       specified as a character vector
% Arguments:
%       scriptName      - name of new MATLAB script
%                       must be a string scalar or a character vector
%       templateNumber  - template number
%                           1 - function_template_minimal.m
%                           2 - function_template_simple.m
%                           3 - function_template.m
%                       must be a positive integer scalar
%                       default == 1
%       varargin    - 'ScriptDir': Directory to place new script
%                   must be a string scalar or a character vector
%                   default == /home/Matlab/{USER}s_Functions
%                   - 'ScriptParentDir': Directory to create default script
%                                           directory
%                   must be a string scalar or a character vector
%                   default == /home/Matlab if exists, or pwd
%                   - 'TemplateDir': Directory to create default script
%                                           directory
%                   must be a string scalar or a character vector
%                   default == same directory as this file
%
% Requires:
%       cd/check_dir.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/function_template_minimal.m
%       cd/function_template_simple.m
%       cd/function_template.m
%
% Used by:

% File History:
% 2019-08-14 Created by Adam Lu
% 2025-08-13 Made scriptParentDir, templateDir optional arguments and
%               changed default for templateDir
% TODO: Make templateNames an optional argument

%% Hard-coded parameters
homeMatlab = '/home/Matlab';
templateNames = {'function_template_minimal', ...
                'function_template_simple', ...
                'function_template'};
templateStr = 'function_template';

%% Default values for optional arguments
templateNumberDefault = 1;
scriptDirDefault = '';
scriptParentDirDefault = '';
templateDirDefault = '';

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
addRequired(iP, 'scriptName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'templateNumber', templateNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ScriptDir', scriptDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ScriptParentDir', scriptParentDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TemplateDir', templateDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, scriptName, varargin{:});
templateNumber = iP.Results.templateNumber;
scriptDir = iP.Results.ScriptDir;
scriptParentDir = iP.Results.ScriptParentDir;
templateDir = iP.Results.TemplateDir;

%% Preparation
% Extract just the file base for the script name
[scriptDirUser, scriptBase] = fileparts(scriptName);

% Decide on parent directory to create default script directory
if isempty(scriptParentDir)
if exist(homeMatlab, 'dir') == 7
    scriptParentDir = homeMatlab;
else
    scriptParentDir = pwd;
end
end

% Decide on the template directory
if isempty(templateDir)
    % Get full path to current file
    currentPath = mfilename('fullpath');

    % Get directory name
    [templateDir, ~, ~] = fileparts(currentPath);
end

% Construct the template path
templatePath = fullfile(templateDir, [templateNames{templateNumber}, '.m']);

%% Construct script path
% Decide on the script directory
if isempty(scriptDir)
    if ~isempty(scriptDirUser)
        scriptDir = scriptDirUser;
    else
        % Get the user name from the server environment
        userName = getenv('USER');

        % Try again if in windows
        if isempty(userName)
            userName = getenv('USERNAME');
        end

        % Capitalize the first letter
        userNameCapitalized = [upper(userName(1)), userName(2:end)];

        % Construct the appropriate subdirectory name
        scriptDirName = [userNameCapitalized, 's_Functions'];

        % Construct the full path to the user's function subdirectory
        scriptDir = fullfile(scriptParentDir, scriptDirName);
    end
end

% Make sure the directory exists
check_dir(scriptDir);

% Construct the full path to the new function
scriptPath = fullfile(scriptDir, [scriptBase, '.m']);

%% Copy the template
% Copy the template to the new path, not allowing any overwrite
if isfile(scriptPath)
    fprintf('%s already exists! Please remove it first!\n', scriptPath);
    return
else
    copyfile(templatePath, scriptPath);
    fprintf('%s copied to %s!\n', templatePath, scriptPath);
end

% Replace all instances of 'function_template' with the new script name
replaceCommand = sprintf('sed -i s/%s/%s/g %s', ...
                        templateStr, scriptBase, scriptPath);
unix(replaceCommand);

%% Open the new script for edit
edit(scriptPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%