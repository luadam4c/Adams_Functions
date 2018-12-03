function parentDir = extract_common_parent (filePaths, varargin)
%% Extracts the common parent directory of a cell array of file paths
% Usage: parentDir = extract_common_parent (filePaths, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       parentDir   - the common parent directory
%                   specified as a character vector
% Arguments:
%       filePaths   - file paths
%                   must be a cell array of character vectors or a string array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
%
% Used by:
%       cd/plot_swd_raster.m

% File History:
% 2018-11-27 Created by Adam Lu
% TODO: Fix
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'filePaths', ...
    @(x) assert(isempty(x) || iscellstr(x) || isstring(x), ...
        ['filePaths must be a a a string array ', ...
        'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, filePaths, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Split all paths by filesep to get path parts
pathParts = cellfun(@(x) split(x, filesep), filePaths, 'UniformOutput', false);

% Count the number of parts from all paths
nParts = cellfun(@numel, pathParts);

% Get the minimum number of parts
minNParts = min(nParts);

% Initialize the last common part to minNParts
ctLastCommon = minNParts;

% Run through all path parts until they become different
for iNPart = 1:minNParts
    % Get this part from all paths
    thisPart = cellfun(@(x) x{iNPart}, pathParts, 'UniformOutput', false);

    % Count the number of unique parts
    nUniqueParts = numel(unique(thisPart));

    % If the number of unique parts is not one, make the previous part the last
    %   common part and exit the loop
    if nUniqueParts ~= 1
        ctLastCommon = iNPart - 1;
        break
    end
end

% Construct the common parent directory
tempCell = join(pathParts{1}(1:ctLastCommon), filesep);
parentDir = tempCell{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%