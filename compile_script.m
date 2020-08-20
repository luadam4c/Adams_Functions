function compiledPath = compile_script (mScriptName, varargin)
%% Compiles a MATLAB script/function
% Usage: compiledPath = compile_script (mScriptName, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compile_script('minEASE');
%
% Outputs:
%       compiledPath     - TODO: Description of compiledPath
%                   specified as a TODO
%
% Arguments:
%       mScriptName - .m script name
%                   must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair 
%                       for all_dependent_files()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_string_end.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2020-08-20 Adapted from minEASE_compile.m
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
addRequired(iP, 'mScriptName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['mScriptName must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, mScriptName, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the all_dependent_files() function
otherArguments = iP.Unmatched;

% Check relationships between arguments
% TODO

%% Preparation
% If an empty file name is provided, return error
if exist(mScriptName, 'file') ~= 2
    mScriptName = force_string_end(mScriptName, '.m');
    fprintf('The file %s cannot be found!\n', mScriptName);
    return
end

%% Do the job
% Find all dependent functions
functionTable = all_dependent_files(mScriptName, otherArguments);

%% Output results
% TODO
compiledPath = mScriptName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
