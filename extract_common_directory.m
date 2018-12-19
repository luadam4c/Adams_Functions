function parentDir = extract_common_directory (paths, varargin)
%% Extracts the common parent directory of a cell array of file paths
% Usage: parentDir = extract_common_directory (paths, varargin)
% Explanation:
%       TODO
% Example(s):
%       extract_common_directory({'a/b/c', 'a/b/c.m', 'a/b/c/d.m'})
% Outputs:
%       parentDir   - the common parent directory
%                   specified as a character vector
% Arguments:
%       paths       - file paths
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/extract_fileparts.m
%
% Used by:
%       cd/extract_fileparts.m
%       cd/plot_swd_raster.m
%       cd/plot_table.m

% File History:
% 2018-11-27 Created by Adam Lu
% 2018-12-18 Now accepts a character array as the input
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'paths', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['paths must be a character array or a string array ', ...
        'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paths, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% If empty, just return
if isempty(paths)
    parentDir = paths;
    return
end

% Extract only the directory part of each path
directories = extract_fileparts(paths, 'directory');

% If a character array, use the directory it is contained in
if ischar(paths)
    parentDir = directories;
    return
end

% Split all directories by filesep to get parts
%   Note: split() can take both a cell and a string as the argument
%           and returns a column cell array
partsByFile = arrayfun(@(x) split(x, filesep), directories, ...
                    'UniformOutput', false);

% Extract the first minNParts elements
partsByFileAligned = extract_subvectors(partsByFile, 'AlignMethod', 'leftadjust');

% Place all parts together in a 2-D cell array
%   Each column is a path
%   Each row is a level
partsArray = horzcat(partsByFileAligned{:});

% Separate parts by level
partsByLevel = extract_rows(partsArray);

% Find the number of unique elements in each row
nUniqueEachLevel = count_unique_elements(partsByLevel);

% Find the first row that has more than one unique element
levelFirstDifference = find(nUnique > 1, 1, 'first');

% Use the previous row number
if isempty(levelFirstDifference)
    levelLastCommon = numel(nUniqueEachLevel);
else
    levelLastCommon = levelFirstDifference - 1;
end

% Construct the common parent directory
tempCell = join(partsByFile{1}(1:levelLastCommon), filesep);
parentDir = tempCell{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rowsExtracted = extract_rows(cellArray)
%% Extracts rows from a 2D cell array
% TODO: Pull out to its own function
% TODO: Think about how to make a extract_columns work the same way
%       with an optional argument

% Count the number of levels
nRows = size(cellArray, 1);

% Separate parts by level
rowsExtracted = arrayfun(@(x) cellArray(x, :), transpose(1:nRows), ...
                            'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nUnique = count_unique_elements(vecs)
%% Counts the number of unique elements in each vector
% TODO: Pull out to its own function
% TODO: Use extract_subvectors with a 'Unique' option

% Extract unique elements
uniqueElements = cellfun(@unique, vecs, 'UniformOutput', false);

% Count the numbers
nUnique = count_samples(uniqueElements);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

pathParts = cellfun(@(x) split(x, filesep), directories, 'UniformOutput', false);

nParts = cellfun(@numel, pathParts);

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

tempCell = join(pathParts{1}(1:ctLastCommon), filesep);
parentDir = tempCell{1};

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
