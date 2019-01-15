function varargout = find_in_list (cand, list, varargin)
%% Returns all indices of a candidate in a list
% Usage: [indices, matched] = find_in_list (cand, list, varargin)
% Explanation:
%   This is the same as find_in_strings.m for text
%   If not text, this is the same as testing list == cand (for now)
%   Use ismember_custom.m to simply test whether cand is in list
%   Use find_first_index_in_list.m to treat each element of cand
%       as a different candidate
%
% Example(s):
%       find_in_list(3, [3, 4, 5, 3])
%       find_in_list('dog', {'dog', 'cat', 'dog'})
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
%                   default == 'parts'
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
%       cd/find_in_strings.m
%
% Used by:
%       cd/compute_combined_trace.m

% File History:
% 2019-01-13 Modified from find_in_strings.m

%% Hard-coded constants
validMatchModes = {'exact', 'parts', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'exact'; % can be parts by default
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
    varargout{:} = find_in_strings(cand, list, 'searchMode', searchMode, ...
                                'IgnoreCase', ignoreCase, 'MaxNum', maxNum, ...
                                'ReturnNan', returnNan);
else
    % Find all indices
    indices = list == cand;

    % Restrict to maxNum
    if ~isempty(maxNum)
        % Make sure maxNum does not exceed the length of indices
        maxNum = min(maxNum, numel(indices));

        % Restrict to maxNum indices
        indices = indices(1:maxNum);
    end

    % Outputs
    varargout{1} = indices;
    if nargin > 2
        varargout{2} = list(indices);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
