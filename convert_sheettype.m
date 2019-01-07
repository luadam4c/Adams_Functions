function convertedFiles = convert_sheettype (varargin)
%% Converts all spreadsheets to desired sheettype (all .xlsx and .xls files to .csv files by default)
% Usage: convertedFiles = convert_sheettype (fileOrDir (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       convertedFiles  - paths of converted files
%                       specified as a character vector 
%                           or a cell array of character vectors
% Arguments:
%       fileOrDir   - (opt) .atf file name or directory containing .atf files
%                   must be a string scalar or a character vector
%                   default == pwd
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/parse_file_or_directory.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-11-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
fileOrDirDefault = pwd;     % convert all spreadsheet files in the present 
                            %   working directory by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add optional inputs to the Input Parser
addOptional(iP, 'fileOrDir', fileOrDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, varargin{:});
fileOrDir = iP.Results.fileOrDir;

%% Preparation
% Parse first argument
[fileDir, fileName, multipleFiles] = parse_file_or_directory(fileOrDir);

%% Do the job
% TODO

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%