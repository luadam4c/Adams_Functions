function fullPaths = extract_fullpaths (files, varargin)
%% Extracts full paths from a files structure array
% Usage: fullPaths = extract_fullpaths (files, varargin)
%
% Outputs:
%       fullPaths   - full path(s) to the files
%                   specified as a column cell array of character vectors
%                       or a character vector
%
% Arguments
%       files       - files structure (may be array) returned by dir()
%                   must be a structure array
%       varargin    - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/all_files.m
%       cd/all_subdirs.m

% File History:
% 2018-09-27 Created by Adam Lu
% 2018-10-03 Added the case when files is a single structure
% 2018-10-03 Renamed extract_fullpaths() -> extract_fullpath()
% 2018-10-24 Renamed extract_fullpath() -> extract_fullpaths()
% 2018-12-26 Added 'ForceCellOutput' as an optional argument

%% Default values for optional arguments
forceCellOutputDefault = false; % don't force output as a cell array by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'files', @isstruct);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, files, varargin{:});
forceCellOutput = iP.Results.ForceCellOutput;

%% Return nothing if the files structure is empty 
%   or does not have required fields
if isempty(files) || ~isfield(files, 'folder') || ~isfield(files, 'name')
    if forceCellOutput
        fullPaths = {};
    else
        fullPaths = '';
    end
    return
end

%% Get the full paths
if numel(files) > 1 || forceCellOutput
    % Get the folders in a column cell array
    folders = transpose({files.folder});

    % Get the names in a column cell array
    names = transpose({files.name});

    % Get the full paths
    fullPaths = cellfun(@(x, y) fullfile(x, y), folders, names, ...
                        'UniformOutput', false);
else
    % Just put the folder and the name together
    fullPaths = fullfile(files.folder, files.name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%