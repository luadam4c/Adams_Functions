function [fullPath, pathType] = construct_fullpath (pathName, varargin)
%% Constructs full path(s) based on file/directory name(s) and optional directory, suffixes or extension
% Usage: [fullPath, pathType] = construct_fullpath (pathName, varargin)
% Explanation:
%       TODO
%
% Examples:
%       fullpaths = construct_fullpath(files, 'Directory', directory);
%       construct_fullpath('funny', 'Directory', '/path/to')
%       construct_fullpath('funny', 'Directory', '/path/to', 'Suffix', 'man')
%       construct_fullpath('funny', 'Directory', '/path/to', 'Suffix', 'man')
%       construct_fullpath('funny', 'Suffixes', {'tall', 'man'})
%
% Outputs:
%       fullPath    - the full path(s) to file(s) or directory(s) constructed
%                   specified as a character vector 
%                       or a column cell array or character vectors
%       pathType    - the type(s) of the path
%                   specified as one of:
%                       'folder'
%                       'file'
%                       or a column cell array of them
%                       
% Arguments:
%       pathName    - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Directory': a full directory path, 
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffixes': suffix(es) to add to fileBase
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%                   default == ''
%                   - 'Extension': file extension to use
%                   must be a string scalar or a character vector
%                   default == whatever is provided by the file name
%   TODO: Change this to 'SuffixNameValuePairs'
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == {'', NaN}
%        
%
% Requires:
%       cd/argfun.m
%       cd/combine_strings.m
%       cd/force_column_cell.m
%       cd/force_string_start.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/atfwrite.m
%       cd/check_dir.m
%       cd/combine_swd_sheets.m
%       cd/locate_dir.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_iox.m
%       cd/plot_grouped_histogram.m
%       cd/plot_histogram.m
%       cd/save_params.m
%       cd/test_var_difference.m
%       ~/RTCl/neuronlaunch.m
%       ~/m3ha/data_dclamp/dclampPassiveFitter.m

% 2017-03-27 Created
% 2017-05-04 Removed ntrials
% 2017-05-04 Added input Parser scheme
% 2017-05-04 Changed paramname, paramvalue -> 'NameValuePairs' 
%                and made it a parameter 
% 2017-05-04 Made 'suffixes' a parameter 
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-10-03 Now accepts file names that are full paths
% 2018-10-03 Added 'Extension' as an optional parameter
% 2018-10-03 Rename construct_fullfilename -> construct_fullpath
% 2018-10-03 Added pathType as an output
% 2018-10-03 Now accepts a cell array of paths as input
% 2018-11-01 Now uses argfun.m, force_column_cell.m, match_format_vector_sets.m
% 2019-10-07 Now uses force_string_start.m on file extension
% TODO: Add 'Prefixes' as an optional argument

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default
directoryDefault = '';              % set later
suffixesDefault = '';
extensionDefault = '';              % set later
nameValuePairsDefault = {'', NaN};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
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
addParameter(iP, 'Suffixes', suffixesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['Suffixes must be either a string/character array ', ...
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
suffixes = iP.Results.Suffixes;
suffixNameValuePairs = iP.Results.NameValuePairs;
extension = iP.Results.Extension;

% TODO: Make optional arguments
prefixes = '';
prefixNameValuePairs = {'', NaN};


%% Preparation
% Match the number of pathNames to directory and extension
pathName = match_format_vector_sets(pathName, directory, 'ForceCellOutputs', false);
pathName = match_format_vector_sets(pathName, extension, 'ForceCellOutputs', false);

%% Do the job for all paths
if iscell(pathName)
    % Force as a column cell array
    pathName = force_column_cell(pathName);
    
    % Match the number of directories and extensions to pathName
    [directory, extension] = ...
        argfun(@(x) match_format_vector_sets(x, pathName, ...
                                            'ForceCellOutputs', true), ...
            directory, extension);

    % Do for each path
    [fullPath, pathType] = ...
        cellfun(@(x, y, z) construct_fullpath_helper(x, verbose, y, ...
                            prefixes, prefixNameValuePairs, ...
                            suffixes, suffixNameValuePairs, z), ...
                pathName, directory, extension, 'UniformOutput', false);
else
    [fullPath, pathType] = ...
        construct_fullpath_helper(pathName, verbose, directory, ...
                                prefixes, prefixNameValuePairs, ...
                                suffixes, suffixNameValuePairs, extension);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fullPath, pathType] = ...
            construct_fullpath_helper (pathName, verbose, directory, ...
                                    prefixes, prefixNameValuePairs, ...
                                    suffixes, suffixNameValuePairs, extension)
%% Create full path to file robustly

% Separate fileDir, fileBase and fileExt from pathName
[fileDir, fileBase, fileExtAuto] = fileparts(pathName);

% Decide on the file extension
if ~isempty(extension)
    % Make sure the file extension, if provided, starts with a dot
    fileExt = force_string_start(extension, '.');
else
    % Use the detected file extension 
    %   Note: this could be empty in case of a directory
    fileExt = fileExtAuto;
end

% Return the path type
if isempty(fileExt)
    pathType = 'folder';
else
    pathType = 'file';
end

% Decide on the directory if not provided
if isempty(directory)
    if ~isempty(fileDir)
        % Use the original file directory
        directory = fileDir;
    else
        % Use the present working directory
        directory = pwd;
    end
end

% Construct final prefix
finalPrefix = combine_strings('EndWithDelimiter', true, ...
                                'Substrings', prefixes, ...
                                'NameValuePairs', prefixNameValuePairs);

% Construct final suffix
finalSuffix = combine_strings('BeginWithDelimiter', true, ...
                                'Substrings', suffixes, ...
                                'NameValuePairs', suffixNameValuePairs);

% Construct full path based on directory, fileBase, final suffix and fileExt
fullPath = fullfile(directory, [finalPrefix, fileBase, finalSuffix, fileExt]);

% Print message
if verbose
    fprintf('Full path to %s: %s\n', pathType, fullPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
            min(cellfun(@isstring, x))) || isstring(x), ...
            ['Suffixes must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
