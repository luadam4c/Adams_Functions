function varargout = find_first_match (candidates, array, varargin)
%% Returns the first matching index and match in an array for candidate(s)
% Usage: [index, matched] = find_first_match (candidates, array, varargin)
% Explanation:
%       TODO
%   Use find_in_strings.m to treat each element of cand
%       as a part of a candidate
% Example(s):
%       find_first_match({'dog'; 'cat'}, ["dog"; "fly"; "cat"])
%       find_first_match(["dog", "fly", "cat"], {'dog'; 'cat'})
% Outputs:
%       index       - the index with matching element in the array
%                   specified as a positive integer array (may contain NaN)
%       matched     - matching element in the array
%                   specified as an array of the same type as elements of array
% Arguments:
%       candidates  - candidates to be matched
%       array       - an array
%       varargin    - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'  - cand must be identical to the members
%                       'parts'  - cand can be parts of the members
%                       'regexp' - cand is a regular expression
%                   default == 'parts'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ismember_custom.m
%
% Used by:
%       cd/convert_to_rank.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Now returns matched
% 

%% Hard-coded parameters
validMatchModes = {'exact', 'parts', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'parts';         % can be parts by default
ignoreCaseDefault = false;          % don't ignore case by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'candidates');
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MatchMode', matchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, candidates, array, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;

%% Do the job
% Find the index in array for each element in candidates
%   Note: If not found, zero will be returned
[~, index] = ismember_custom(candidates, array, ...
                            'MatchMode', matchMode, 'IgnoreCase', ignoreCase);

% Get all matched elements if requested
if nargout >= 2
    matched = array(index);
end

%% Deal with outputs
varargout{1} = index;
if nargout >= 2
    varargout{2} = matched;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
