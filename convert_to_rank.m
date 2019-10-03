function ranks = convert_to_rank (array, varargin)
%% Creates a positive integer array from an array showing the ranks of each element
% Usage: ranks = convert_to_rank (array, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       convert_to_rank([20, 3, 3, 10, 10, 20, 3])
%       convert_to_rank({'dog'; 'cat'; 'cat'; 'dog'})
%       convert_to_rank(["dog", "cat", "cat", "dog"])
%       convert_to_rank({'dog'; 'cat'}, 'Ranked', ["dog"; "fly"; "cat"])
%       convert_to_rank(["dog", "fly", "cat"], 'Ranked', {'dog'; 'cat'})
%
% Outputs:
%       ranks        - ranks of each element
%                   specified as a positive integer array (may contain NaN)
%
% Arguments:
%       array       - an array that can be passed into unique()
%       varargin    - 'RankedElements': ranked elements
%                   must be an array of the same type as array
%                   default == unique(array)
%                   - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'  - cand must be identical to the members
%                       'parts'  - cand can be parts of the members
%                       'regexp' - cand is a regular expression
%                   default == 'exact'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_in_list.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_default_grouping.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-15 Now uses find_in_list.m instead
% 

%% Hard-coded parameters
validMatchModes = {'exact', 'parts', 'regexp'};

%% Default values for optional arguments
rankedElementsDefault = [];         % set later
matchModeDefault = 'exact';         % exact matches by default
ignoreCaseDefault = false;          % don't ignore case by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RankedElements', rankedElementsDefault);
addParameter(iP, 'MatchMode', matchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, array, varargin{:});
rankedElements = iP.Results.RankedElements;
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;

% Keep unmatched arguments for the unique() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default ranks of elements
if isempty(rankedElements)
    rankedElements = unique(array, otherArguments{:});
end

%% Do the job
% Count the number of ranks
nRanks = numel(rankedElements);

% Do for each rank
indicesEachRank = cell(nRanks, 1);
parfor iRank = 1:nRanks
    % Get the current ranked element
    if iscell(rankedElements)
        rankedElementThis = rankedElements{iRank};
    else
        rankedElementThis = rankedElements(iRank);
    end

    % Find the indices with this ranked element in the array
    indicesEachRank{iRank} = ...
        find_in_list(rankedElementThis, array, ...
                    'MatchMode', matchMode, 'IgnoreCase', ignoreCase);
end

% Initialize a ranks array with the same size as array
ranks = zeros(size(array));

% Store the ranks of these elements
for iRank = 1:nRanks
    ranks(indicesEachRank{iRank}) = iRank;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% This is too slow:
ranks = find_first_match(array, rankedElements, ...
                        'SearchMode', searchMode, 'IgnoreCase', ignoreCase);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
