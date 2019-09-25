function [tables1, tables2, distinctParts] = load_matching_sheets (suffix1, suffix2, varargin)
%% Loads spreadsheets with matching strings before given suffixes (incomplete)
% Usage: [tables1, tables2, distinctParts] = load_matching_sheets (suffix1, suffix2, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [tables1, tables2, distinctParts] = load_matching_sheets('pulses', 'SWDs')
%
% Outputs:
%       tables1     - first list of tables
%                   specified as a cell array of tables
%       tables2     - second list of tables
%                   specified as a cell array of tables
%       distinctParts   - distinct parts of file names that identifies tables
%                       specified as a cell array of character vectors
%
% Arguments:
%       suffix1     - suffix for the first list of tables
%                   must be a string scalar or a character vector
%       suffix2     - suffix for the second list of tables
%                   must be a string scalar or a character vector
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
%       cd/plot_relative_events.m

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
addRequired(iP, 'suffix1', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'suffix2', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, suffix1, suffix2, varargin{:});
directory = iP.Results.Directory;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Do the job
% Get all files for suffix1
[~, paths1] = ...
    all_files('Directory', directory, 'Keyword', pathBase, ...
                'Suffix', suffix1, 'Extension', sheetType, ...
                'ForceCellOutput', true);

% Get all matching files for suffix2
[~, paths2, distinctParts] = ...
    find_matching_files(paths1, 'Directory', directory, ...
                        'Suffix', suffix2, 'Extension', sheetType);

% Read all tables
[tables1, tables2] = ...
    argfun(@(x) cellfun(@readtable, x, 'UniformOutput', false), ...
            paths1, paths2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%