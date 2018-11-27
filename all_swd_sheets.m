function [swdSheetFiles, swdSheetPaths] = all_swd_sheets (varargin)
%% Returns all files ending with '_SWDs.csv' under a directory recursively
% Usage: [swdSheetFiles, swdSheetPaths] = all_swd_sheets (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       swdSheetFiles   - file structure(s) for the SWD spreadheet files
%                       specified as a structure array with fields:
%                           name
%                           folder
%                           date
%                           bytes
%                           isdir
%                           datenum
%       swdSheetPaths   - full path(s) to the SWD spreadheet files
%                       specified as a column cell array of character vectors
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Directory': directory to look for SWD table files
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SheetType': sheet type;
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'csv'
%                   
% Requires:
%       cd/all_files.m
%       cd/issheettype.m
%       cd/print_or_show_message.m
%
% Used by:
%       cd/plot_swd_raster.m

% File History:
% 2018-11-27 Moved from plot_swd_raster.m
% 

%% Hard-coded parameters
swdStr = '_SWDs';               % string in file names for SWD spreadsheets

%% Default values for optional arguments
verboseDefault = true;
directoryDefault = '';          % set later
sheetTypeDefault = 'csv';       % default spreadsheet type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
directory = iP.Results.Directory;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Do the job
% Find all SWD spreadsheet files in directory
[swdSheetFiles, swdSheetPaths] = ...
    all_files('Verbose', verbose, 'Recursive', true, ...
                'Directory', directory, ...
                'Suffix', swdStr, 'Extension', ['.', sheetType]);

% Exit function if no spreadsheet files are found
if isempty(swdSheetPaths)
    message = sprintf(['There are no SWD spreadsheets of the', ...
                        ' ending %s.%s in the directory: %s'], ...
                        swdStr, sheetType, directory);
    mTitle = 'No SWD spreadsheets found warning';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon, ...
                            'MessageMode', 'show', 'Verbose', verbose, ...
                            'CreateMode', 'replace');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%