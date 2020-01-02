function isMatch = is_matching_string (strList, cand, varargin)
%% Returns whether each element in a list of strings matches a candidate
% Usage: isMatch = is_matching_string (strList, cand, varargin)
% Explanation:
%   There are three main search modes (parameter 'SearchMode'):
%       'substrings': allows the candidate to be a substring or substrings 
%                       of a match in strList. (cf. matches())
%       'exact': candidate must be an exact match
%       'regexp': candidate is a regular expression
%   The latter two cases are similar to strcmp()/strcmpi() or regexp()/regexpi()
%   To return indices and matches, use ismatch.m or find_in_strings.m
%
% Example(s):
%       strs1 = {'Mark''s fish', 'Peter''s fish', 'Katie''s sealion'};
%       strs2 = ["Mark's fish", "Peter's fish", "Katie's sealion"];
%       is_matching_string(strs1, 'fish')
%       is_matching_string(strs2, 'Peter')
%       is_matching_string(strs2, {'Katie', 'lion'})
%       is_matching_string(strs1, "fish", 'MaxNum', 1)
%       is_matching_string(strs1, "Fish", 'IgnoreCase', 1)
%       is_matching_string(strs2, 'Fish', 'IgnoreCase', false)
%       is_matching_string(strs1, "sealion", 'SearchMode', 'exact')
%       is_matching_string(strs2, 'sea', 'SearchMode', 'exact', 'ReturnNaN', true)
%       is_matching_string(strs1, "sea.*", 'SearchMode', 'reg')
%       is_matching_string(strs2, 'sea.*', 'SearchMode', 'reg')
%
% Outputs:
%       isMatch     - whether each member of strList contains 
%                       all parts of the candidate
%                   specified as a logical array
%
% Arguments:
%       strList     - a list of strings
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       cand        - candidate string or substring(s)
%                       If cand is a list of substrings, all substrings must 
%                           exist in the string to be matched
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       varargin    - 'MatchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'     - cand must be identical to the members
%                       'parts'     - cand can be parts of the members
%                       'prefix'    - cand is the prefix of the members
%                       'keyword'   - cand is a part of the members
%                       'suffix'    - cand is the suffix of the members
%                       'regexp'    - cand is a regular expression
%                   if search mode is 'exact' or 'regexp', 
%                       cand cannot have more than one elements
%                   default == 'parts'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchOnce': whether to match only once
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/compute_combined_trace.m
%       cd/create_error_for_nargin.m
%       cd/find_custom.m
%
% Used by:
%       cd/ismatch.m

% File History:
% 2019-01-10 Moved from find_in_strings.m
% 2019-01-11 Removed 2nd and 3rd outputs
% 2020-01-01 Added 'prefix', 'keyword', 'suffix' as match modes
% TODO: Implement matchOnce for all conditions

%% Hard-coded constants
validMatchModes = {'exact', 'parts', 'prefix', 'keyword', 'suffix', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'parts';     % match string parts by default
ignoreCaseDefault = false;      % case-sensitive by default
matchOnceDefault = false;       % no restrictions by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'strList', ...          % a list of strings
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['strList must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'cand', ...             % a string/substrings of interest
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['cand must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MatchMode', matchModeDefault, ...
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchOnce', matchOnceDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strList, cand, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;
matchOnce = iP.Results.MatchOnce;

%% Preparation
% Check relationships between arguments
if ~ischar(cand) && numel(cand) > 1 && ~strcmp(matchMode, 'parts')
    error(['Second input cannot have more than one members if ', ...
            '''MatchMode'' is not ''parts''!']);
end

%% Find the indices
switch matchMode
case 'parts'        % cand can be a substring or a list of substrings
    % Find indices that contain cand in strList
    if ischar(cand) || isstring(cand) && numel(cand) == 1
        % Test whether each element of strList contain the substring
        isMatch = contains(strList, cand, 'IgnoreCase', ignoreCase);    
    else                % if cand is a list of substrings
        % Test whether each element contains each substring
        if iscell(cand)
            hasEachCand = cellfun(@(x) contains(strList, x, ...
                                        'IgnoreCase', ignoreCase), ...
                                cand, 'UniformOutput', false);
        elseif isstring(cand)
            hasEachCand = arrayfun(@(x) contains(strList, x, ...
                                        'IgnoreCase', ignoreCase), ...
                                cand, 'UniformOutput', false);
        else
            error('cand is unrecognized!');
        end

        % Test whether each element contains all substrings
        isMatch = compute_combined_trace(hasEachCand, 'all');
    end
case 'exact'        % cand must be an exact match
    % Test whether each string in strList matches the candidate exactly
    if ignoreCase
        isMatch = strcmpi(strList, cand);
    else
        isMatch = strcmp(strList, cand);
    end
case {'prefix', 'keyword', 'suffix', 'regexp'}
    % Construct a regular expression
    switch matchMode
        case 'prefix'
            regExp = ['^', cand, '.*'];
        case 'keyword'
            regExp = ['.*', cand, '.*'];
        case 'suffix'
            regExp = ['.*', cand, '$'];
        case 'regexp'
            regExp = cand;
        otherwise
            error('matchMode unrecognized!')
    end

    % Decide on options for regexpi
    if matchOnce
        matchOption = 'once';
    else
        matchOption = 'all';
    end
    if ignoreCase
        caseOption = 'ignorecase';
    else
        caseOption = 'matchcase';
    end

    % Returns the matched string that matches the regular expression
    %   for each string in strList
    matchedStrings = regexpi(strList, regExp, 'match', caseOption, matchOption);

    % Test whether each string in strList matches the regular expression
    isMatch = ~isemptycell(matchedStrings);
otherwise
    error('matchMode unrecognized!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

validMatchModes = {'exact', 'substrings'};

if strcmpi(matchMode, 'exact')             % String must be exact
elseif strcmpi(matchMode, 'substrings')    % String can be a substring 
                                            % or a cell array of substrings
end

indicesEachStr = contains(strListMod, strMod);

% Find the indices that contain the substring
indicesarray = strfind(strListMod, strMod);
indices = find(~cellfun(@isempty, indicesarray), maxNum);

% Convert substring to lower case if IgnoreCase is set to true
if ignoreCase
    strMod = lower(cand);
else
    strMod = cand;
end

indicesarray = strfind(strListMod, strMod(k));    
indicesEachStr{k} = find(~cellfun(@isempty, indicesarray));

% Convert each substring to lower case if IgnoreCase is set to true
if ignoreCase
    strMod = cellfun(@lower, cand, 'UniformOutput', false);
else
    strMod = cand;
end

nStrs = numel(cand);
indicesEachStr = cell(1, nStrs);
for k = 1:nStrs
    % Test whether each element of strList contain the substring
    isMatch = contains(strListMod, cand(k), 'IgnoreCase', ignoreCase);

    % 
    indicesEachStr{k} = find(isMatch);
end

% Find the indices that contain all substrings by intersection
indices = intersect_over_cells(indicesEachStr);
%       cd/intersect_over_cells.m

% If more than maxNum indices found, 
%   restrict to the first maxNum indices
if length(indices) > maxNum
    indices = indices(1:maxNum);
end

indices = find(~cellfun(@isempty, startIndices), maxNum);

@(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
            min(cellfun(@isstring, x))) || isstring(x), ...
            ['First input must be either a string/character array ', ...
            'or a cell array of strings/character arrays!']));
@(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || ...
            min(cellfun(@isstring, x))), ...
            ['Second input must be a cell array ', ...
            'of strings/character arrays!']));

% Construct a cell array with the same size as strList 
%   but with cand duplicated throughout
str_cell = cell(size(strList));    % a cell array to store copies of cand
for k = 1:numel(strList)
    % Make the kth element the same as cand
    str_cell{k} = cand;
end

% Find indices that corresponds to cand exactly in strList, 
%   case-insensitive if IgnoreCase is set to true
if ignoreCase
    indices = find_custom(cellfun(@strcmpi, strList, str_cell), ...
                            maxNum, 'ReturnNan', returnNan);
else
    indices = find_custom(cellfun(@strcmp, strList, str_cell), ...
                            maxNum, 'ReturnNan', returnNan);
end

strListMod = cellfun(@lower, strList, 'UniformOutput', false);

% Convert each string to lower case if IgnoreCase is set to true
if ignoreCase
    strListMod = lower(strList);
else
    strListMod = strList;
end

if iscell(cand)        % if cand is a list of substrings
    % Test whether each element contains each substring
    if iscell(cand)
        hasEachCand = cellfun(@(x) contains(strList, x, ...
                                    'IgnoreCase', ignoreCase), ...
                            cand, 'UniformOutput', false);
    else
        hasEachCand = arrayfun(@(x) contains(strList, x, ...
                                    'IgnoreCase', ignoreCase), ...
                            cand, 'UniformOutput', false);
    end

    % Test whether each element contains all substrings
    isMatch = compute_combined_trace(hasEachCand, 'all');
else                    % if cand is a single substring
    % Test whether each element of strList contain the substring
    isMatch = contains(strList, cand, 'IgnoreCase', ignoreCase);    
end

if iscell(cand) && ...
    (strcmpi(matchMode, 'exact') || strcmpi(matchMode, 'regexp'))
    error(['Second input cannot be a cell array if ', ...
            'MatchMode'' is ''exact'' or ''rexexp''!']);
end

% Make sure strList is not a character array
if ischar(strList)
    strList = {strList};
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
