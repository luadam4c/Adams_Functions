function update_file_base_in_matfiles (varargin)
%% Updates slice bases in .mat files to match the file name if changed
% Usage: update_file_base_in_matfiles (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Arguments:
%       varargin    - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/all_files.m
%       cd/extract_fileparts.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-11-30 Created by Adam Lu
% TODO: Generalize to update_var_in_matfile.m

%% Hard-coded parameters

%% Default values for optional arguments
directoryDefault = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;

%% Do the job
% Find all .mat files
[~, matPaths] = all_files('Directory', directory, 'Extension', 'mat');

% Do nothing if no mat file
if numel(matPaths) == 0
    return
end

% Update all files
cellfun(@update_one_file, matPaths)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_one_file (matPath)
%% Updates one file

% Create a matFile object for writing
m = matfile(matPath, 'Writable', true);

% Extract new slice base
newFileBase = extract_fileparts(matPath, 'base');

% Get all
%   Note: this is much faster than who or whos
allVarNames = fieldnames(m);

% Update slice base in the file
if ismember('sliceBase', allVarNames)
    m.sliceBase = newFileBase;
elseif ismember('fileBase', allVarNames)
    m.fileBase = newFileBase;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
