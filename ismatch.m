function varargout = ismatch (list, cand, varargin)
%% Returns whether each element in a list matches a candidate
% Usage: [isMatch, indices, matched] = ismatch(list, cand, varargin)
% Explanation:
%   This function uses is_matching_string.m for strings
%       and isequaln() for everything else
%   This function optionally returns the indices as the second output
%       and the matched elements as the third output.
%
% Example(s):
%       ismatch([3; 4; 4; 3], 4)
%       ismatch(1:10, 4)
%       ismatch([1 NaN 2 1 1 NaN 3 NaN 2 3], 2)
%       ismatch([1 NaN 2 1 1 NaN 3 NaN 2 3], NaN)
%       ismatch({'', '12', 'dg', ''}, '')
%       strs1 = {'Mark''s fish', 'Peter''s fish', 'Katie''s sealion'};
%       strs2 = ["Mark's fish", "Peter's fish", "Katie's sealion"];
%       ismatch(strs1, 'fish')
%       ismatch(strs2, 'Peter')
%       ismatch(strs1, "sealion")
%       ismatch(strs1, "fish", 'MaxNum', 1)
%       ismatch(strs1, "Fish", 'IgnoreCase', 1)
%       ismatch(strs2, 'Fish', 'IgnoreCase', false)
%       ismatch(strs2, {'Katie', 'lion'}, 'MatchMode', 'parts')
%       ismatch(strs1, "sea.*", 'MatchMode', 'reg')
%       ismatch(strs2, 'sea.*', 'MatchMode', 'reg')
%       [~, indices] = ismatch(strs2, 'sea', 'ReturnNaN', true)
%
% Outputs:
%       isMatch     - whether each member of the list contains 
%                       all parts of the candidate
%                   specified as a logical array
%       indices     - indices of the members of the list matching the candidate
%                       could be empty (or NaN if 'ReturnNan' is true)
%       members    - members of the list corresponding to those indices
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       list        - an array
%       cand        - candidate or parts of the candidate
%       varargin    - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'  - cand must be identical to the members
%                       'parts'  - cand can be parts of the members
%                       'regexp' - cand is a regular expression
%                   default == 'exact'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNum': maximum number of indices to find
%                   must be empty or a positive integer scalar
%                   default == numel(list)
%                   - 'ReturnNan': Return NaN instead of empty if nothing found
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_custom.m
%       cd/is_matching_string.m
%       cd/istext.m
%
% Used by:
%       cd/plot_measures.m
%       cd/find_in_strings.m
%       cd/ismember_custom.m
%       cd/unique_groups.m

% File History:
% 2019-01-11 Modified from is_matching_strings.m
% 2019-08-21 Now uses isequaln() instead of isequal() for matching

%% Hard-coded constants
validMatchModes = {'exact', 'parts', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'exact'; % must be exact copy by default
ignoreCaseDefault = false;  % case-sensitive by default
maxNumDefault = [];         % set later
returnNanDefault = false;   % return empty if nothing is found by default

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
addRequired(iP, 'list');
addRequired(iP, 'cand');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MatchMode', matchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...       % maximum number of indices
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));
addParameter(iP, 'ReturnNan', returnNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, list, cand, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;
returnNan = iP.Results.ReturnNan;

% Check relationships between arguments
if strcmp(matchMode, 'regexp') && ~istext(cand)  
    error('Second input must be text if ''MatchMode'' is ''regexp''!');
end

%% Preparation
% Make sure list is not a character array
if ischar(list)
    list = {list};
end

% Set the maximum number of indices if not provided
if nargout >= 2 && isempty(maxNum)
    % Count the number of indices in list
    nIndices = numel(list);

    % Set maximum number to be the total number of indices
    maxNum = nIndices;
end

%% Do the job
if istext(list)
    % Use is_matching_string.m
    isMatch = is_matching_string(list, cand, 'MatchMode', matchMode, ...
                                    'IgnoreCase', ignoreCase);
else
    % Test whether each member is a match to cand
    if iscell(list)
        isMatch = cellfun(@(x) isequaln(x, cand), list);
    else
        isMatch = arrayfun(@(x) isequaln(x, cand), list);
    end
end

%% Find all indices of members in list that is a match
if nargout >= 2
    indices = find_custom(isMatch, maxNum, 'ReturnNan', returnNan);
end

%% Return the matched members too
if nargout >= 3
    if ~isempty(indices) && any(isnan(indices))
        matched = NaN;
    elseif ~isempty(indices) 
        if numel(indices) > 1
            matched = list(indices);
        else
            matched = list{indices};
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
