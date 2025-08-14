function [output1, output2] = function_template (reqarg1, reqarg2, reqstr3, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: [output1, output2] = function_template (reqarg1, reqarg2, reqstr3, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%       output2     - TODO: Description of output2
%                   specified as a TODO
%
% Side Effects:
%       Plots TODO
%
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       reqarg2     - TODO: Description of reqarg2
%                   must be a string vector or a cell array of character vectors
%       reqstr3     - TODO: Description of reqstr3
%                   must be a string scalar or a character vector
%       vecs        - vectors to TODO
%                   must be a numeric array or a cell array of numeric arrays
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
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'num3': TODO: Description of num3
%                   must be a positive integer scalar
%                   default == TODO
%                   - 'FigTypes': figure type(s) for saving
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - 'SheetType': sheet type for TODO
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
%                   - 'LineStyle': line style of TODO
%                   must be an unambiguous, case-insensitive match to one of: 
%                       '-'     - solid line
%                       '--'    - dashed line
%                       ':'     - dotted line
%                       '-.'    - dash-dotted line
%                       'none'  - no line
%                   default == '-'
%                   - 'strs5': string(s) to TODO
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == 'nostrs5'
%                   - 'folder6': directory to TODO, e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNo': file number (sweep number)
%                   must be an integer scalar
%                   default == 0
%                   - 'Vectors': vectors
%                   must be a numeric matrix or a cell array of numeric vectors
%                   default == {}
%                   - Any other parameter-value pair for TODO()
%       COMMON arguments:
%       filebase    - (opt) base of filename (without extension)
%                       e.g. 'A100110_0008_18'
%                   must be a string scalar or a character vector
%                   default == 'unnamed'
%
% Requires:
%       TODO: place any custom functions used in this function/script here
%       /TODO:dir/TODO:file
%       /home/Matlab/Adams_Functions/addpath_custom.m
%       /home/Matlab/Adams_Functions/isaninteger.m
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/check_dir.m
%       /home/Matlab/Adams_Functions/check_subdir.m
%       /home/Matlab/Adams_Functions/create_error_for_nargin.m
%       /home/Matlab/Adams_Functions/get_matlab_year.m
%       /home/Matlab/Adams_Functions/struct2arglist.m
%
% Used by:
%       TODO: place any custom functions/scripts that uses this function here
%       /TODO:dir/TODO:file

% File History:
% 20XX-XX-XX Created by TODO or Adapted from TODO
% 20XX-XX-XX Added TODO
% 20XX-XX-XX Changed TODO
% 2017-05-09 Now uses isfigtype.m to validate figTypes
% 2017-05-21 Updated for consistent naming conventions
% 2017-05-22 Updated for consistent indentation conventions
% 2017-05-22 Placed default values in the beginning of file
% 2018-01-24 Added isdeployed
% 2018-02-09 Now uses mfilename as function name
% 2018-05-16 In argument checks, changed min -> all
% 2018-06-10 isdeployed should not be needed for homedirectory
% 2018-10-03 Now uses isfolder() instead of exist()
% 2018-12-17 Updated the response when number of arguments fall short
% 

%% Hard-coded parameters
% TODO: Any number that was hard-coded should be place here
validStrings = {'str1', 'str2'};

%% Default values for optional arguments
optarg1Default  = [];                   % default TODO: Description of optarg1
optstr2Default  = '';                   % default TODO: Description of optstr2
param1Default   = [];                   % default TODO: Description of param1
flag2Default    = false;                % whether to TODO by default
num3Default     = [];                   % default TODO: Description of num3
                                        % TODO: Why this number?
figTypesDefault = 'png';                % default figure type(s) for saving
sheetTypeDefault = 'xlsx';              % default spreadsheet type
lineStyleDefault = '-';                 % default line style
strs5Default = '';                      % default string(s) to TODO
folder6Default = '';                    % default directory to TODO
fileNoDefault = 0;                      % default file number for saving output
vectorsDefault = {};

%% Directory (ies) for placing outputs
directoryArray = {'/dir1/', '/dir2/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with directories
% Find home directory across servers
if isfolder(fullfile(pwd, 'Adams_Functions'))
    homedirectory = pwd;
elseif isfolder('/media/adamX/DIRECTORY')
    homedirectory = '/media/adamX/DIRECTORY/';
elseif isfolder('/scratch/al4ng/DIRECTORY')
    homedirectory = '/scratch/al4ng/DIRECTORY';
else
    homedirectory = pwd;
end

%% If not compiled, add directories to search path for required functions
%   Note: If addpath is used, adding `if ~isdeployed` is important to 
%          avoid errors if the function is used as part of a compiled program
%   Note: addpath takes a long time, so use addpath_custom for checking
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for TODO
    addpath_custom(fullfile(functionsDirectory, 'Adams_Functions'));
    addpath_custom(fullfile(functionsDirectory, 'Downloaded_Functions'));
end

%% Deal with arguments
% Check number of required arguments
if nargin < 2    % TODO: 2 might need to be changed
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'reqarg1', ...                  % TODO: Description of reqarg1
    % TODO: validation function %);
addRequired(iP, 'reqarg2', ...                  % TODO: Description of reqarg2
    @(x) assert(iscellstr(x) || isstring(x), ...
                ['Second input must be a cell array of character arrays ', ...
                'or a string array!']));
addRequired(iP, 'reqstr3', ...                  % TODO: Description of reqstr3
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(iP, 'vecs', ...
    @(x) assert(isnum(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'optarg1', optarg1Default, ...
    % TODO: validation function %);
addOptional(iP, 'optstr2', optstr2Default, ...
    @(x) any(validatestring(x, validStrings)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);
addParameter(iP, 'flag2', flag2Default, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'num3', num3Default, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'num3', num3Default, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['num3 must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'num3', num3Default, ...
    @(x) assert((ischar(x) || isstring(x)) && strcmpi(x, 'all') || ...
                isscalar(x) && isaninteger(x) && x > 0, ...
                ['num3 must be either ''all'', "all" ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'num3', num3Default, ...
    @(x) assert(isnum(x) || iscellnumeric(x), ...
                ['num3 must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'LineStyle', lineStyleDefault, ...
    @(x) all(islinestyle(x, 'ValidateMode', true)));
addParameter(iP, 'strs5', strs5Default, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strs5 must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'folder6', folder6Default, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'FileNo', fileNoDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
addParameter(iP, 'Vectors', vectorsDefault, ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['Vectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Read from the Input Parser
parse(iP, reqarg1, reqarg2, reqstr3, vecs, varargin{:});
optarg1 = iP.Results.optarg1;
optstr2 = validatestring(iP.Results.optstr2, validStrings);
param1 = iP.Results.param1;
flag2 = iP.Results.flag2;
num3 = iP.Results.num3;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
                                    % in /home/Matlab/Adams_Functions
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
                                    % in /home/Matlab/Adams_Functions
[~, lineStyle] = islinestyle(iP.Results.LineStyle, 'ValidateMode', true);
                                    % in /home/Matlab/Adams_Functions
strs5 = iP.Results.strs5;
folder6 = iP.Results.folder6;
fileNo = iP.Results.FileNo;
vectors = iP.Results.Vectors;

% Keep unmatched arguments for the TODO function
otherArguments = iP.Unmatched;
otherArguments = struct2arglist(iP.Unmatched);

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

% Set dependent argument defaults
if isempty(folder6)
    folder6 = pwd;
end

% Check relationships between arguments
% TODO

% Make sure a directory is an existing full path
folder6 = construct_and_check_fullpath(folder6);

% Check if needed output directories exist
check_dir(outFolders);                  % in /home/Matlab/Adams_Functions
check_subdir(outFolder, subdirs);       % in /home/Matlab/Adams_Functions

%% Check for files of a given data type
% Look for files of the given data type in the data directory
%   This returns a row structure array
dataFiles = dir(['*.', dataType]);

% Extract the file names as a row cell array
dataFileNames = {dataFiles.name};

% Convert to full file names
dataFullFileNames = ...
    cellfun(@(x) fullfile(folder, x), dataFileNames, 'UniformOutput', false);

%% Common variables
% Get the current MATLAB release year
matlabYear = get_matlab_year;

% Create log file name
logFileName = fullfile(outFolder, fprintf('%s.log', mfilename));
fid = fopen(logFileName, 'w');
fprintf(fid, 'TODO\n');
fclose(fid);

% Get the time stamp in the format: yymmddThhmm
tempStamp = datestr(clock, 30);     % current time stamp
dateStamp = [tempStamp(1:8)];       % only use date
dateTimeStamp = [tempStamp(1:end-2)];  % take off seconds

%% Preparation
% TODO

%% Do the job
% TODO

%% Output results
% TODO


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

for k = 1:numel(directoryArray)
    if exist(fullfile(outFolder, directoryArray{k}), 'dir') ~= 7
        newdirectory = fullfile(outFolder, directoryArray{k});
        mkdir(newdirectory);
        fprintf('New directory made: %s\n\n', newdirectory);
    end
end

if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory is made: %s\n\n', outFolder);
end

@(x) assert(ischar(x) || iscell(x) && (all(cellfun(@ischar, x)) ...
    || all(cellfun(@isstring, x))) || isstring(x) , ...
    ['strs5 must be either a string/character array ', ...
        'or a cell array of strings/character arrays!']));

if isfolder(fullfile(pwd, 'Adams_Functions'))
    functionsdirectory = pwd;
elseif isfolder('/home/Matlab/')
    functionsDirectory = '/home/Matlab/';
elseif isfolder('/scratch/al4ng/Matlab/')
    functionsDirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsDirectory does not exist!');
end

@(x) assert(iscell(x) && (all(cellfun(@ischar, x)) ...
            || all(cellfun(@isstring, x))), ...
            ['Second input must be a cell array ', ...
            'of strings or character arrays!']));

error(['Not enough input arguments, ', ...
        'type ''help %s'' for usage'], mfilename);
eval(sprintf('help %s', mfilename));

%                   must be a string/character array or 
%                       a cell array of strings/character arrays
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['strs5 must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
