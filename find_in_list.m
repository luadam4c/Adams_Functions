function varargout = find_in_list (cand, list, varargin)
%% Returns all indices of a candidate in a list
% Usage: [indices, matched] = find_in_list (cand, list, varargin)
% Explanation:
%   This is the same as find_in_strings.m for text
%   If not text, this is the same as find(list == cand) (for now)
%   Use ismember_custom.m to simply test whether cand is in list
%   Use find_first_match.m to treat each element of cand
%       as a different candidate
%
% Example(s):
%       [i, e] = find_in_list(3, [3, 4, 5, 3])
%       [i, e] = find_in_list([3, 4], [3, 4, 5, 3])
%       [i, e] = find_in_list([3, 4], [3, 4, 5, 3], 'MaxNum', 1)
%       [i, e] = find_in_list([2, 3, 4], [3, 4, 5, 3], 'MaxNum', 1)
%       [i, e] = find_in_list([2, 3, 4], [3, 4, 5, 3], 'MaxNum', 1, 'ReturnNaN', true)
%       [i, e] = find_in_list('dog', {'dog', 'cat', 'dog'})
%
% Outputs:
%       indices     - indices of list matching the candidate
%                       could be empty (or NaN if 'ReturnNan' is true)
%                   specified as a positive integer array
%       elements    - matched elements of list corresponding to those indices
%                   specified as a cell array if more than one indices 
%                       or the element if only one index; or an empty string
% Arguments:
%       cand        - candidate
%                       If cand is a list of substrings, all substrings must 
%                           exist in the string to be matched
%       list        - a list
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
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/find_custom.m
%       cd/find_in_strings.m
%
% Used by:
%       cd/compute_combined_trace.m
%       cd/convert_to_rank.m

% File History:
% 2019-01-13 Modified from find_in_strings.m
% 2020-02-04 Fixed list == cand to find(list == cand)
% 2020-02-04 Now allows cand to be an array of candidates
% 2020-02-04 Fixed varargout and nargout

%% Hard-coded constants
validMatchModes = {'exact', 'parts', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'exact'; % exact matches by default
ignoreCaseDefault = false;  % whether to ignore case by default
maxNumDefault = [];         % will be changed to numel(list)
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
addRequired(iP, 'cand');
addRequired(iP, 'list');

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
parse(iP, cand, list, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;
maxNum = iP.Results.MaxNum;
returnNan = iP.Results.ReturnNan;

%% Preparation
% Check relationships between arguments
if strcmp(matchMode, 'regexp') && ~istext(cand)  
    error('First input must be text if ''MatchMode'' is ''regexp''!');
end

% Translate match mode to search mode
if strcmp(matchMode, 'parts')
    searchMode = 'substrings';
else
    searchMode = matchMode;
end

%% Use ismatch.m
if istext(cand)
    % Use find_in_strings.m
    [varargout{1:nargout}] = ...
        find_in_strings(cand, list, 'searchMode', searchMode, ...
                        'IgnoreCase', ignoreCase, 'MaxNum', maxNum, ...
                        'ReturnNan', returnNan);
else
    % Find all indices
    if numel(cand) == 1
        indices = find_custom(list == cand, maxNum, 'ReturnNaN', returnNan);
    else
        if ~isempty(maxNum) && maxNum == 1 && returnNan
            indices = array_fun(@(x) find_custom(list == x, 1, ...
                                                'ReturnNaN', true), cand);
        else
            try
                indices = array_fun(@(x) find_custom(list == x, maxNum, ...
                                                    'ReturnNaN', returnNan), ...
                                    cand, 'UniformOutput', true);
            catch
                indices = array_fun(@(x) find_custom(list == x, maxNum, ...
                                                    'ReturnNaN', returnNan), ...
                                    cand, 'UniformOutput', false);
            end
        end        
    end

    % Outputs
    varargout{1} = indices;
    if nargout >= 2
        if iscell(indices)
            varargout{2} = array_fun(@(x) list(x), indices, ...
                                    'UniformOutput', false);
        elseif any(isnan(indices))
            varargout{2} = extract_subvectors(list, 'Indices', indices);
        else
            varargout{2} = list(indices);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Restrict to maxNum
if ~isempty(maxNum)
    % Make sure maxNum does not exceed the length of indices
    maxNum = min(maxNum, numel(indices));

    % Restrict to maxNum indices
    indices = indices(1:maxNum);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
