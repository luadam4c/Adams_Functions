function fullPath = extract_fullpath(files)
%% Extracts full paths from a files structure array
% Usage: fullPath = extract_fullpath(files)
%
% Outputs:
%       fullPath    - full path(s) to the files
%                   specified as a column cell array of character vectors
%                       or a character vector
%
% Arguments
%       files       - files structure (may be array) returned by dir()
%                   must be a structure array
%
% Used by:
%       cd/all_files.m
%       cd/all_subdirs.m

% File History:
% 2018-09-27 Created by Adam Lu
% 2018-10-03 Added the case when files is a single structure
% 2018-10-03 Renamed extract_fullpaths() -> extract_fullpath()

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

% Add required inputs to the Input Parser
addRequired(iP, 'files', @isstruct);

% Read from the Input Parser
parse(iP, files);

%% Return nothing if the files structure is empty 
%   or does not have required fields
if isempty(files) || ~isfield(files, 'folder') || ~isfield(files, 'name')
    fullPath = '';
    return
end

%% Get the full paths
if numel(files) > 1
    % Get the folders in a column cell array
    folders = transpose({files.folder});

    % Get the names in a column cell array
    names = transpose({files.name});

    % Get the full paths
    fullPath = cellfun(@(x, y) fullfile(x, y), folders, names, ...
                        'UniformOutput', false);
else
    % Just put the folder and the name together
    fullPath = fullfile(files.folder, files.name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%