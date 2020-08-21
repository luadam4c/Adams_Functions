function compiledPath = compile_script (mScriptName, varargin)
%% Compiles a MATLAB script/function
% Usage: compiledPath = compile_script (mScriptName, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compiledPath = compile_script('minEASE');
%
% Outputs:
%       compiledPath    - full path to compiled file
%                       specified as a character vector
%
% Arguments:
%       mScriptName - .m script name
%                   must be a string scalar or a character vector
%       varargin    - 'SaveFlag': whether to save spreadsheets
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair 
%                       for all_dependent_files()
%
% Requires:
%       cd/all_dependent_files.m
%       cd/create_error_for_nargin.m
%       cd/force_string_end.m
%
% Used by:
%       cd/minEASE_compile.sh

% File History:
% 2020-08-20 Adapted from minEASE_compile.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
saveFlagDefault = false;

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
addRequired(iP, 'mScriptName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['mScriptName must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, mScriptName, varargin{:});
saveFlag = iP.Results.SaveFlag;

% Keep unmatched arguments for the all_dependent_files() function
otherArguments = iP.Unmatched;

%% Preparation
% If an empty file name is provided, return error
if exist(mScriptName, 'file') ~= 2
    mScriptName = force_string_end(mScriptName, '.m');
    fprintf('The file %s cannot be found!\n', mScriptName);
    return
end

%% Do the job
% Find all dependent files
fileListTable = all_dependent_files(mScriptName, 'SaveFlag', saveFlag, ...
                                    otherArguments);

% Extract the full paths
fullPaths = fileListTable.fullPath;

% Create command to compile script
command = sprintf('mcc -m -v ');
command = [command, char(join(strcat("-a ", fullPaths)))];
command = [command, ' ', mScriptName];

% Compile script
eval(command);

%% Output results
% Output compiled path
compiledPath = fullfile(pwd, mScriptName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
