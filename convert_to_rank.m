function ranks = convert_to_rank (array, varargin)
%% Creates a positive integer array from an array showing the ranks of each element
% Usage: ranks = convert_to_rank (array, varargin)
% Explanation:
%       TODO
% Example(s):
%       convert_to_rank([20, 3, 3, 10, 10, 20, 3])
%       convert_to_rank({'dog'; 'cat'; 'cat'; 'dog'})
%       convert_to_rank(["dog", "cat", "cat", "dog"])
%       convert_to_rank({'dog'; 'cat'}, 'Ranked', ["dog"; "fly"; "cat"])
%       convert_to_rank(["dog", "fly", "cat"], 'Ranked', {'dog'; 'cat'})
% Outputs:
%       ranks        - ranks of each element
%                   specified as a positive integer array (may contain NaN)
% Arguments:
%       array       - an array that can be passed into unique()
%       varargin    - 'RankedElements': ranked elements
%                   must be an array of the same type as array
%                   default == unique(array)
%                   - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - str must be identical to 
%                                           an element in cellArray
%                       'substrings'    - str can be a substring or 
%                                           a cell array of substrings
%                       'regexp'        - str is considered a regular expression
%                   if searchMode is 'exact' or 'regexp', 
%                       str cannot be a cell array
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_index_in_list.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/create_default_grouping.m

% File History:
% 2019-01-09 Created by Adam Lu
% 

%% Hard-coded parameters
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
rankedElementsDefault = [];         % set later
searchModeDefault = 'exact';        % search exact matches by default
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
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, array, varargin{:});
rankedElements = iP.Results.RankedElements;
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;

% Keep unmatched arguments for the unique() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Set default ranks of elements
if isempty(rankedElements)
    rankedElements = unique(array, otherArguments{:});
end

%% Do the job
% Find the rank for each element in array
ranks = find_index_in_list(array, rankedElements, ...
                            'SearchMode', searchMode, 'IgnoreCase', ignoreCase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
