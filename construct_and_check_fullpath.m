function [fullPath, pathExists] = construct_and_check_fullpath (pathName, varargin)
%% Constructs the full path to the file or directory and checks whether it exists
% Usage: [fullPath, pathExists] = construct_and_check_fullpath (pathName, varargin)
% Explanation:
%       TODO
%
% Examples:
%       [abfFullfilename, fileExists] = ...
%           construct_and_check_fullpath(pathName, 'Extension', '.abf');
%       if ~fileExists
%           return
%       end
%
% Outputs:
%       fullPath    - the full path(s) to file(s) or directory(s) constructed
%                   specified as a character vector 
%                       or a column cell array or character vectors
%       pathExists  - whether the path(s) exists
%                   specified as a column logical array
%
% Arguments:
%       pathName    - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for construct_fullpath()
%
% Requires:
%       cd/check_fullpath.m
%       cd/construct_fullpath.m
%
% Used by:
%       cd/apply_to_all_subdirs.m
%       cd/all_files.m
%       cd/load_neuron_outputs.m
%       cd/load_params.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_parse_mat.m
%       cd/parse_atf_swd.m
%       cd/parse_iox.m
%       cd/parse_gas_trace.m

% File History: 
% 2017-04-11 - Moved from plot_traces_abf.m
% 2018-10-03 - Renamed construct_and_check_fullfilename -> 
%                   construct_and_check_fullfilename; 
% 2018-10-03 - Now uses construct_fullpath.m
% 2018-10-03 - Now uses isfile()
% 2018-10-03 - Renamed construct_and_check_fullfilename -> 
%                   construct_and_check_fullpath
% 2018-11-21 - Updated Used by
% 2019-11-28 - Now passes other arguments to construct_fullpath.m
% 

%% Default values for optional arguments
verboseDefault = false;             % don't print to standard output by default

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

% Add required inputs to an input Parser
addRequired(iP, 'pathName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['pathName must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the input Parser
parse(iP, pathName, varargin{:});
verbose = iP.Results.Verbose;

% Keep unmatched arguments for the construct_fullpath() function
otherArguments = iP.Unmatched;

%% Create full path(s) to file(s) robustly
[fullPath, pathType] = construct_fullpath(pathName, 'Verbose', verbose, ...
                                            otherArguments);

%% Check whether the full path(s) exist(s)
pathExists = check_fullpath(fullPath, 'Verbose', verbose, ...
                            'PathType', pathType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
