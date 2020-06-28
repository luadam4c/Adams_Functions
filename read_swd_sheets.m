function [swdTables, filePaths] = read_swd_sheets (varargin)
%% Loads SWD tables from SWD spreadsheets
% Usage: [swdTables, filePaths] = read_swd_sheets (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       swdTables   - TODO: Description of swdTables
%                   specified as a TODO
%       filePaths   - TODO: Description of filePaths
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'Recursive': whether to search recursively
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FilePaths': names of '_SWDs.csv' sheets to load
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from pwd
%                   - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Suffix': suffix the file name must have
%                   must be a string scalar or a character vector
%                   default == set in all_swd_sheets.m
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   
% Requires:
%       cd/all_swd_sheets.m
%       cd/construct_and_check_fullpath.m
%       cd/issheettype.m
%       cd/read_swd_sheet.m
%
% Used by:
%       cd/plot_swd_histogram.m
%       cd/plot_swd_raster.m

% File History:
% 2018-11-27 Created by Adam Lu
% TODO: Combine with read_swd_sheet.m

%% Hard-coded parameters

%% Default values for optional arguments
recursiveDefault = true;        % search recursively by default
directoryDefault = '';        % set later
filePathsDefault = {};          % detect from pwd by default
verboseDefault = false;         % print to standard output by default
suffixDefault = '';             % set in all_swd_sheets.m
sheetTypeDefault = 'csv';       % default spreadsheet type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Recursive', recursiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FilePaths', filePathsDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
recursive = iP.Results.Recursive;
directory = iP.Results.Directory;
filePaths = iP.Results.FilePaths;
verbose = iP.Results.Verbose;
suffix = iP.Results.Suffix;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Initialize output
swdTables = {};

% Decide on the files to use
if isempty(filePaths)
    % Decide on the directory if not provided
    if isempty(directory)
        % Use the present working directory
        directory = pwd;
    end

    % Find all '_SWDs.csv' sheets in the directory recursively
    [~, filePaths] = all_swd_sheets('Verbose', verbose, ...
                                    'Recursive', recursive, ...
                                    'Directory', directory, ...
                                    'Suffix', suffix, ...
                                    'SheetType', sheetType);

    % Return usage message if no files found
    if isempty(filePaths)
        fprintf('Type ''help %s'' for usage\n', mfilename);
        return
    end
elseif ischar(filePaths)
    % Place in cell array
    filePaths = {filePaths};
end

% Check if each path exists
[filePaths, pathExists] = ...
    construct_and_check_fullpath(filePaths, 'Directory', directory);

% Return if not all paths exist
if ~all(pathExists)
    return
end

%% Do the job
% Read in the tables from the files
if iscell(filePaths)
    swdTables = cellfun(@read_swd_sheet, filePaths, 'UniformOutput', false);
else
    swdTables = read_swd_sheet(filePaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
