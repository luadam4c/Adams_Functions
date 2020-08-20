function [output1, output2] = function_template (reqarg1, reqarg2, reqstr3, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: [output1, output2] = function_template (reqarg1, reqarg2, reqstr3, varargin)
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%       output2     - TODO: Description of output2
%                   specified as a TODO
% Arguments:    
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       reqarg2     - TODO: Description of reqarg2
%                   must be a cell array of strings or character arrays
%       reqstr3     - TODO: Description of reqstr3
%                   must be a string scalar or a character vector
%       optarg1     - (opt) TODO: Description of optarg1
%                   must be a TODO
%                   default == TODO
%       optstr2     - (opt) TODO: Description of optstr2
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'str1'        - TODO: Description of str1
%                       'str2'        - TODO: Description of str2
%                   default == 'defaultstr'
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - 'flag2': whether to TODO
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'num3': TODO: Description of num3
%                   must be a positive integer
%                   default == TODO
%                   - 'figtypes4': figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'strs5': string(s) to TODO
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%                   default == 'nostrs5'
%                   - 'folder6': directory to TODO, e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%       COMMON arguments:
%       filebase    - (opt) base of filename (without extension)
%                       e.g. 'A100110_0008_18'
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%
% Requires:
%       TODO: place any custom functions used in this function/script here
%       /TODO:dir/TODO:file
%       cd/isfigtype.m
%
% Used by:    
%       TODO: place any custom functions/scripts that uses this function here
%       /TODO:dir/TODO:file
%
% File History:
% 201X-XX-XX Created by TODO or Adapted from TODO
% 201X-XX-XX Added TODO
% 201X-XX-XX Changed TODO
% 2017-05-09 Now uses isfigtype.m to validate figtypes
% 2017-05-21 Updated for consistent naming conventions
% 2017-05-22 Updated for consistent indentation conventions
% 

%% Parameters used in the body of the function
% (Any number that was hard-coded should be place here)

%% Directory (ies) for placing outputs
directoryArray = {'/dir1/', '/dir2/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with directories
% Find home directory across servers
if exist('/media/adamX/DIRECTORY/', 'dir') == 7
    homeDirectory = '/media/adamX/DIRECTORY/';
elseif exist('/scratch/al4ng/DIRECTORY/', 'dir') == 7
    homeDirectory = '/scratch/al4ng/DIRECTORY/';
else
    error('Valid homeDirectory does not exist!');
end

%% Add directories to search path for required functions across servers
if exist('/home/Matlab/', 'dir') == 7
    functionsDirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsDirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsDirectory does not exist!');
end
addpath_custom(fullfile(functionsDirectory, '/Adams_Functions/')); % for isfigtype.m
addpath_custom(fullfile(functionsDirectory, '/Marks_Functions/')); % for TODO

%% Deal with arguments
% Check number of required arguments
if nargin < 2    % TODO: 2 might need to be changed
    error(['Not enough input arguments, ', ...
            'type ''help function_template'' for usage']);
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'reqarg1', ...                  % TODO: Description of reqarg1
    % TODO: validation function %);
addRequired(iP, 'reqarg2', ...                  % TODO: Description of reqarg2
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
            || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'reqstr3', ...                  % TODO: Description of reqstr3
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after 2016B

% Add optional inputs to the Input Parser
addOptional(iP, 'optarg1', % TODO %, ...        % TODO: Description of optarg1
    % TODO: validation function %);
addOptional(iP, 'optstr2', 'defaultstr', ...    % TODO: Description of optstr2
    @(x) any(validatestring(x, {'str1', 'str2'})));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', % TODO %, ...        % TODO: Description of param1
    % TODO: validation function %);
addParameter(iP, 'flag2', false, ...            % whether to TODO
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'num3', % TODO %, ...          % TODO: Description of num3
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'figtypes4', 'png', ...        % figure type(s) for saving
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'strs5', 'nostrs5', ...        % string(s) to TODO
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) ...
        || min(cellfun(@isstring, x))) || isstring(x) , ...
        ['strs5 must be either a string/character array ', ...
            'or a cell array of strings/character arrays!']));
addParameter(iP, 'folder6', '', ...             % directory to TODO
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, reqarg1, reqarg2, reqstr3, varargin{:});
optarg1 = iP.Results.optarg1;
optstr2 = validatestring(iP.Results.optstr2, {'str1', 'str2'});
param1 = iP.Results.param1;
flag2 = iP.Results.flag2;
num3 = iP.Results.num3;
[~, figtypes4] = isfigtype(iP.Results.figtypes4, 'ValidateMode', true);
strs5 = iP.Results.strs5;
folder6 = iP.Results.folder6;

% Set default arguments
if isempty(folder6)
    folder6 = pwd;
end

% Check relationships between arguments
% TODO

%% Check if needed output directories exist
for k = 1:numel(directoryArray)
    if exist(fullfile(outfolder, directoryArray{k}), 'dir') ~= 7
        newdirectory = fullfile(outfolder, directoryArray{k});
        mkdir(newdirectory);
        fprintf('New directory made: %s\n\n', newdirectory);
    end
end

%% TODO: Describe Big Step 1
% TODO: Describe Small Step 1

% TODO: Describe Small Step 2

%% TODO: Describe Big Step 2
% TODO: Describe Small Step 1

% TODO: Describe Small Step 2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction1 ()
%% A summary of what subfunction1 does


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subfunction2 ()
%% A summary of what subfunction2 does


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
TODO: Place older versions of the code that you want to save here, 
        in case you need it back in the future

%                - 'FigType4': figure type for saving
%                could be anything recognised by the saveas function such as 'png' or 'fig', etc.
%                default == 'png'
addParameter(iP, 'FigType4', 'png', ...            % figure type for saving
    @(x) any(validatestring(x, {'png', 'fig', 'm', 'mfig', ...
                    'jpeg', 'tiff', 'tiffn', 'meta', ...
                    'bmpmono', 'bmp', 'bmp16m', 'bmp256', ...
                    'hdf', 'pbm', 'pbmraw', 'pcxmono', ...
                    'pcx24b', 'pcx256', 'pcx16', 'pgm', ...
                    'pgmraw', 'ppm', 'ppmraw'})));
addParameter(iP, 'FigType4', 'png', ...            % figure type for saving; could be 'png' or 'fig', etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
FigType4 = iP.Results.FigType4;

addParameter(iP, 'figtypes4', {}, ...            % a cell array of figure types for saving
    @(x) iscell(x) && min(isfigtype(x, 'ValidateMode', true)));

%}
