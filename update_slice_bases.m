function update_slice_bases (varargin)
%% Updates slice bases to file bases for all slice data .mat files in a directory
% Usage: update_slice_bases (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Arguments:
%       varargin    - 'Directory': the directory(ies) to search in
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == pwd
%
% Requires:
%
% Used by:
%       cd/m3ha_oscillations_analyze.m

% File History:
% 2020-07-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
directoryDefault = '';          % construct_and_check_fullpath('') == pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;

%% Do the job
% Find all .mat files
[~, allMatPaths] = all_files('Directory', directory, 'Extension', 'mat');

% Find all parsed results .mat files
[~, allParsedMatPaths] = all_files('Directory', directory, ...
                                'Suffix', 'parsed', 'Extension', 'mat');

% Select only the slice data .mat files
allDataMatPaths = setdiff(allMatPaths, allParsedMatPaths);

% Update all slice data .mat files
cellfun(@update_slice_base_helper, allDataMatPaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_slice_base_helper(matPath)

% Read the .mat file with writable permissions
m = matfile(matPath, 'Writable', true);

% Extract the file base
sliceBaseNew = extract_fileparts(matPath, 'base');

% Extract the stored slice base
sliceBaseOld = m.sliceBase;

% Update to file base if not a match
if ~strcmp(sliceBaseOld, sliceBaseNew)
    m.sliceBase = sliceBaseNew;
    fprintf('Slice base updated to %s for %s!\n', sliceBaseNew, matPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%