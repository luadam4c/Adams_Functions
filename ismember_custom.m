function varargout = ismember_custom (cand, list, varargin)
%% Returns whether a particular candidate is a member of a list
% Usage: [isMember, index] = ismember_custom (cand, list, varargin)
% Explanation:
%   This function is similar to ismember() but allows different search modes
%       for strings
%   There are three main search modes (parameter 'MatchMode'):
%       'substrings': allows the candidate to be a substring or substrings 
%                       of a match in list.
%       'exact': candidate must be an exact match
%       'regexp': candidate is a regular expression
%
% Example(s):
%       ismember_custom({'dog'; 'cat'}, ["dog"; "fly"; "cat"])
%       ismember_custom(["dog", "fly", "cat"], {'dog'; 'cat'})
%
%       strs1 = {'Mark''s fish', 'Peter''s fish', 'Katie''s sealion'};
%       strs2 = ["Mark's fish", "Peter's fish", "Katie's sealion"];
%       ismember_custom('fish', strs1)
%       ismember_custom('Peter', strs2)
%       ismember_custom({'Katie', 'lion'}, strs2)
%       ismember_custom({{'Katie', 'lion'}}, strs2)
%       ismember_custom("fish", strs1, 'MaxNum', 1)
%       ismember_custom("Fish", strs1, 'IgnoreCase', 1)
%       ismember_custom('Fish', strs2, 'IgnoreCase', false)
%       ismember_custom("sealion", strs1, 'MatchMode', 'exact')
%       ismember_custom('sea', strs2, 'MatchMode', 'exact', 'ReturnNaN', true)
%       ismember_custom("sea.*", strs1, 'MatchMode', 'reg')
%       ismember_custom('sea.*', strs2, 'MatchMode', 'reg')
%
% Outputs:
%       isMember    - whether cand matches at least a member of list
%                   specified as a logical vector
%       index       - the first matching index in the list
%                   specified as a positive integer (or NaN) vector
%
% Arguments:
%       cand        - candidate(s)
%       list        - a list
%       varargin    - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'     - cand must be identical to the members
%                       'parts'     - cand can be parts of the members
%                       'prefix'    - cand is the prefix of the members
%                       'keyword'   - cand is a part of the members
%                       'suffix'    - cand is the suffix of the members
%                       'regexp'    - cand is a regular expression
%                   default == 'parts'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/ismatch.m
%       cd/istext.m
%
% Used by:
%       cd/extract_vars.m
%       cd/find_first_match.m
%       cd/is_on_path.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_plot_figure07.m

% File History:
% 2019-01-10 Modified from find_in_strings.m
% 2020-01-01 Added 'prefix', 'keyword', 'suffix' as match modes
% 2020-05-13 Now accepts list as any cell array

%% Hard-coded constants
validMatchModes = {'exact', 'parts', 'prefix', 'keyword', 'suffix', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'parts';         % can be parts by default
ignoreCaseDefault = false;          % whether to ignore case by default

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
addRequired(iP, 'cand');
addRequired(iP, 'list');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MatchMode', matchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, cand, list, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;

%% Do the job
% Test whether each member of list matches cand
if istext(cand) && (ignoreCase || ~strcmp(matchMode, 'exact'))
    % Make sure cand is either a string or a cell array
    if ischar(cand)
        cand = {cand};
    end

    % Test whether each member of cand has a match in list
    if nargout <= 1
        % Just test whether each member of cand has a match in list
        isMember = ...
            array_fun(@(x) is_in_strings(x, list, matchMode, ignoreCase), cand);
    else
        % Find the first index of the matching element in list too
        [isMember, index] = ...
            array_fun(@(x) is_in_strings(x, list, matchMode, ignoreCase), cand);
    end
else
    if iscell(list)
        if iscell(cand)
            [isMember, index] = cellfun(@(x) is_an_element(x, list), cand);
        else
            [isMember, index] = is_an_element(cand, list);
        end
    else
    % Test whether each candidate is a member in the list
        if nargout >= 2
            % Find the index too
            %   Note: If not found, zero will be returned by ismember()
            [isMember, index] = ismember(cand, list);

            % Make all zeros NaNs instead
            index(index == 0) = NaN;
        else
            % Just test whether each candidate is a member in the list
            isMember = ismember(cand, list);
        end
    end
end

%% Deal with outputs
varargout{1} = isMember;
if nargout >= 2
    varargout{2} = index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = is_in_strings (cand, list, matchMode, ignoreCase)
% Returns whether list contains a match of cand 
%   and return the index of the first match

% Test whether each member of list matches cand
if nargout == 1
    isMatch = ...
        ismatch(list, cand, 'MaxNum', 1, 'ReturnNan', true, ...
                'MatchMode', matchMode, 'IgnoreCase', ignoreCase);
elseif nargout == 2
    [isMatch, index] = ...
        ismatch(list, cand, 'MaxNum', 1, 'ReturnNan', true, ...
                'MatchMode', matchMode, 'IgnoreCase', ignoreCase);
end

% See if any member is a match
isMember = any(isMatch);

% Outputs
varargout{1} = isMember;
if nargout >= 2
    varargout{2} = index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [isAnElement, index] = is_an_element (cand, cellArray)
%% Returns whether cand is an element of a cell array

% Initialize output
isAnElement = false;
index = NaN;

% Run through each element
for iElement = 1:numel(cellArray)
    if isequal(cand, cellArray{iElement})
        isAnElement = true;
        index = iElement;
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
