function directory = locate_dir (candidates, varargin)
%% Locate the first directory that exists out of a list of candidates
% Usage: directory = locate_dir (candidates, varargin)
% Example(s):
%       candidates = {'/tmp/data/m3ha/', '/media/adamX/m3ha/', ...
%                       '/home/adam/m3ha/', '/scratch/al4ng/m3ha/'};
%       homeDirectory = locate_dir(candidates, ...
%                                       'DirectoryType', 'home directory');
% Outputs:
%       directory   - the first directory that exists
%                   specified as a character vector
% Arguments:    
%       candidates  - candidate directories
%                   must be a cell array of character arrays
%       varargin    - 'DirectoryType': the type of directory
%                   must be a string scalar or a character vector
%                   default == 'custom directory'
%                   - 'ContainedSubdir': a subdirectory that must be contained
%                                           under the directory
%                   must be a string scalar or a character vector
%                   default == none
%
% Requires:
%       cd/construct_fullpath.m
%
% Used by:
%       cd/locate_functionsdir.m
%       cd/m3ha_locate_homedir.m

% File History:
% 2018-10-04 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
directoryTypeDefault = 'custom directory';
containedSubdirDefault = '';

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
addRequired(iP, 'candidates', ...               % candidate directories
    @(x) iscellstr(x) || isstring(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DirectoryType', directoryTypeDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'ContainedSubdir', containedSubdirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, candidates, varargin{:});
directoryType = iP.Results.DirectoryType;
containedSubdir = iP.Results.ContainedSubdir;

%% Preparation
% Count the number of candidates
nCandidates = numel(candidates);

% Make sure there are candidates
if nCandidates == 0
    error('There are no candidates provided!!\n\n');
end

% Initialize output
directory = '';

%% Loop through all candidates
for iCandidate = 1:nCandidates
    % Get the current candidate
    candidate = candidates{iCandidate};

    % Construct full path and check if it exists
    fullPath = construct_fullpath(candidate);

    % Check if it is an existing directory
    if isfolder(fullPath) && isfolder(fullfile(fullPath, containedSubdir))
        directory = fullPath;
        return
    end
end

% Return error if no directory located
if isempty(directory)
    error('Valid %s does not exist!', directoryType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%