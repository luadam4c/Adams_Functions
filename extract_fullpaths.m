function fullPaths = extract_fullpaths(filesStructure)
%% Extracts full paths from a files structure
% Usage: fullPaths = extract_fullpaths(filesStructure)
%
% Used by:
%       cd/subdirs.m

% File History:
% 2018-09-27 Created by Adam Lu

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
addRequired(iP, 'filesStructure', @isstruct);

% Read from the Input Parser
parse(iP, filesStructure);

%% Get the full paths
fullPaths = cellfun(@(x, y) fullfile(x, y), ...
                    {filesStructure.folder}, {filesStructure.name}, ...
                    'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%