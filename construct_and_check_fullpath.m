function [fullPath, pathExists] = construct_and_check_fullpath (pathName, varargin)
%% Constructs the full path to the file or directory and checks whether it exists
% Usage: [fullPath, pathExists] = construct_and_check_fullpath (pathName, varargin)
% Examples:
%       [abfFullfilename, fileExists] = ...
%           construct_and_check_fullpath(pathName, 'Extension', '.abf');
%       if ~fileExists
%           return 
%       end
% Outputs:
%       fullPath    - the full path(s) to file(s) or directory(s) constructed
%                   specified as a character vector 
%                       or a column cell array or character vectors
%       pathExists  - whether the path(s) exists
%                   specified as a column logical array
% Arguments:
%       pathName    - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': a full directory path, 
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'Suffices': suffix(ces) to add to fileBase
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%                   default == ''
%                   - 'Extension': file extension to use
%                   must be a string scalar or a character vector
%                   default == whatever is provided by the file name
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == {'', NaN}
% Requires:
%       cd/check_fullpath.m
%       cd/construct_fullpath.m
%
% Used by:
%       cd/all_files.m
%       cd/all_subdirs.m
%       cd/load_neuron_outputs.m
%       cd/load_params.m
%       cd/m3ha_parse_mat.m

% File History: 
% 2017-04-11 - Moved from plot_traces_abf.m
% 2018-10-03 - Renamed construct_and_check_fullfilename -> 
%                   construct_and_check_fullfilename; 
% 2018-10-03 - Now uses construct_fullpath.m
% 2018-10-03 - Now uses isfile()
% 2018-10-03 - Renamed construct_and_check_fullfilename -> 
%                   construct_and_check_fullpath
% 2018-11-21 - Updated Used by
% 

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
directoryDefault = '';
sufficesDefault = '';
extensionDefault = '';              % set later
nameValuePairsDefault = {'', NaN};

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

% Add required inputs to an input Parser
addRequired(iP, 'pathName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['pathName must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'Suffices', sufficesDefault, ...
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['Suffices must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'NameValuePairs', nameValuePairsDefault, ...
    @(x) assert(iscell(x) && numel(x) == 2 ...
            && (ischar(x{1}) || iscell(x{1}) || isstring(x{1})) ...
            && isnumeric(x{2}), ...
        ['NameValuePairs must be a 2-element cell array whose ', ...
            'first element is a string/char array or cell array ', ...
            'and whose second element is a numeric array!']));

% Read from the input Parser
parse(iP, pathName, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
suffices = iP.Results.Suffices;
extension = iP.Results.Extension;
namevaluepairs = iP.Results.NameValuePairs;

%% Create full path(s) to file(s) robustly
[fullPath, pathType] = construct_fullpath(pathName, ...
                                'Verbose', verbose, ...
                                'Directory', directory, ...
                                'Suffices', suffices, ...
                                'Extension', extension, ...
                                'NameValuePairs', namevaluepairs);


%% Check whether the full path(s) exist(s)
pathExists = check_fullpath(fullPath, 'Verbose', verbose, ...
                            'PathType', pathType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fprintf('Full path to abf file: %s\n', abfFullFileName);
if exist(abfFullFileName, 'file') ~= 2
    fprintf('This abf file doesn''t exist!');
    return;
end

%% Constructs the full file name to an abf file robustly based on pathName and display message if doesn't exist

if exist(abfFullFileName, 'file') ~= 2
    fileExists = false;
else
    fileExists = true;
end

% Get the file directory and file base, ignoring the file extension
[fileDir, filebase, ~] = fileparts(pathName);

% Append .abf to the file base, then reconstruct the full file name
abffilename = fullfile(fileDir, strcat(filebase, '.abf')];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%