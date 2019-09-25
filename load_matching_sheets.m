function [tables1, tables2, distinctParts] = load_matching_sheets (suffix1, suffix2, varargin)
%% Loads spreadsheets with matching strings before given suffixes (incomplete)
% Usage: [tables1, tables2, distinctParts] = load_matching_sheets (suffix1, suffix2, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       tables1     - TODO: Description of tables1
%                   specified as a TODO
%       tables2     - TODO: Description of tables2
%                   specified as a TODO
%       distinctParts    - TODO
%
% Arguments:
%       suffix1     - TODO: Description of suffix1
%                   must be a TODO
%       suffix2     - TODO: Description of suffix1
%                   must be a TODO
%       varargin    - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/all_files.m
%       cd/create_error_for_nargin.m
%       cd/find_matching_files.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-11 Created by Adam Lu
% 

% TODO: Make optional arguments
pathBase = '';
sheetType = 'csv';

%% Hard-coded parameters

%% Default values for optional arguments
directoryDefault = '';          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'suffix1');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, suffix1, varargin{:});
directory = iP.Results.Directory;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Do the job
% Get all files for suffix1
[~, paths1] = ...
    all_files('Directory', directory, 'Keyword', pathBase, ...
                'Suffix', suffix1, 'Extension', sheetType);

% Get all matching files for suffix2
[~, paths2] = ...
    find_matching_files(paths1, 'Directory', directory, ...
                        'Suffix', suffix2, 'Extension', sheetType);

% Read all tables
[tables1, tables2] = ...
    argfun(@(x) cellfun(@readtable, x, 'UniformOutput', false), ...
            paths1, paths2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Extract distinct file parts
distinctParts = extract_distinct_fileparts(paths1);

% Look for corresponding files for suffix2
[~, paths2] = ...
    cellfun(@(x) all_files('Prefix', x, 'Directory', directory, ...
                    'Suffix', suffix2, 'Extension', sheetType), ...
            distinctParts, 'UniformOutput', false);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%