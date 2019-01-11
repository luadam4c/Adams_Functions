function varargout = is_matching_string (strList, cand, varargin)
%% Returns whether each element in a list of strings matches a candidate
% Usage: [isMatch, indices, matched] = is_matching_string(strList, cand, varargin)
% Explanation:
%   There are three main search modes (parameter 'SearchMode'):
%       'substrings': allows the candidate to be a substring or substrings 
%                       of a match in strList.
%       'exact': candidate must be an exact match
%       'regexp': candidate is a regular expression
%   The latter two cases are similar to strcmp()/strcmpi() or regexp()/regexpi()
%   However, this function optionally returns the indices as the second output
%       and the matched elements as the third output.
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
%       indices     - indices of strList matching the candidate
%                       could be empty (or NaN if 'ReturnNan' is true)
%       elements    - elements of strList corresponding to those indices
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       strList     - a list of strings
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       cand        - candidate string or substring(s)
%                       If cand is a list of substrings, all substrings must 
%                           exist in the string to be matched
%                   must be a string/character array or 
%                       a cell array of strings/character arrays
%       varargin    - 'SearchMode': the search mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - cand must be identical to 
%                                           an element in strList
%                       'substrings'    - cand can be a substring or 
%                                           a list of substrings
%                       'regexp'        - cand is considered a regular expression
%                   if search mode is 'exact' or 'regexp', 
%                       cand cannot have more than one elements
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNum': maximum number of indices to find
%                   must be empty or a positive integer scalar
%                   default == numel(strList)
%                   - 'ReturnNan': Return NaN instead of empty if nothing found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/compute_combined_trace.m
%       cd/create_error_for_nargin.m
%       cd/find_custom.m
%
% Used by:
%       cd/find_in_strings.m
%       cd/ismember_custom.m

% File History:
% 2019-01-10 Moved from find_in_strings.m

%% Hard-coded constants
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
searchModeDefault = 'substrings';       % default search mode
ignoreCaseDefault = false;              % whether to ignore case by default
maxNumDefault = [];                     % will be changed to numel(strList)
returnNanDefault = false;   % whether to return NaN instead of empty 
                            %   if nothing found by default

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
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, strList, cand, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;
returnNan = iP.Results.ReturnNan;

% Check relationships between arguments
if ~ischar(cand) && numel(cand) > 1 && ...
    (strcmp(searchMode, 'exact') || strcmp(searchMode, 'regexp'))
    error(['Second input cannot have more than one members if ', ...
            '''SearchMode'' is ''exact'' or ''regexp''!']);
end

%% Prepare for the search
% Make sure strList is not a character array
if ischar(strList)
    strList = {strList};
end

% Set the maximum number of indices if not provided
if nargout >= 2 && isempty(maxNum)
    % Count the number of indices in strList
    nIndices = numel(strList);

    % Set maximum number to be the total number of indices
    maxNum = nIndices;
end

%% Find the indices
switch searchMode
case 'substrings'   % cand can be a substring or a list of substrings
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
case 'regexp'       % cand is considered a regular expression
    % Returns the starting index that matches the regular expression
    %   for each string in strList
    if ignoreCase
        startIndices = regexpi(strList, cand);
    else
        startIndices = regexp(strList, cand);
    end

    % Test whether each string in strList matches the regular expression
    isMatch = ~isemptycell(startIndices);
end

%% Find all indices of strings in strList that is a match
if nargout >= 2
    indices = find_custom(isMatch, maxNum, 'ReturnNan', returnNan);
end

%% Return the matched elements too
if nargout >= 3
    if ~isempty(indices) && any(isnan(indices))
        matched = NaN;
    elseif ~isempty(indices) 
        if numel(indices) > 1
            matched = strList(indices);
        else
            matched = strList{indices};
        end
    else
        matched = '';
    end
end

%% Outputs
varargout{1} = isMatch;
if nargout >= 2
    varargout{2} = indices;
end
if nargout >= 3
    varargout{3} = matched;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

validSearchModes = {'exact', 'substrings'};

if strcmpi(searchMode, 'exact')             % String must be exact
elseif strcmpi(searchMode, 'substrings')    % String can be a substring 
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
    (strcmpi(searchMode, 'exact') || strcmpi(searchMode, 'regexp'))
    error(['Second input cannot be a cell array if ', ...
            'SearchMode'' is ''exact'' or ''rexexp''!']);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
