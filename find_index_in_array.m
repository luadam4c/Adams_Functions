function varargout = find_index_in_array (candidates, array, varargin)
%% Returns the index for candidate(s) in an array, treating strings and character arrays as the same
% Usage: [index, matched] = find_index_in_array (candidates, array, varargin)
% Explanation:
%       TODO
% Example(s):
%       find_index_in_array({'dog'; 'cat'}, ["dog"; "fly"; "cat"])
%       find_index_in_array(["dog", "fly", "cat"], {'dog'; 'cat'})
% Outputs:
%       index       - the index with matching element in the array
%                   specified as a positive integer array (may contain NaN)
%       matched     - matching element in the array
%                   specified as an array of the same type as elements of array
% Arguments:
%       candidates  - candidates to be matched
%       array       - an array
%       varargin    - 'SearchMode': the search mode if candidates are strings
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'         - candidate must be identical to 
%                                           an element in array
%                       'substrings'    - candidate can be a substring or 
%                                           a cell array of substrings
%                       'regexp'        - candidate is considered 
%                                           a regular expression
%                   if searchMode is 'exact' or 'regexp', 
%                       str cannot be a cell array
%                   default == 'substrings'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/find_in_strings.m
%
% Used by:
%       cd/convert_to_rank.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Now returns matched
% 

%% Hard-coded parameters
validSearchModes = {'exact', 'substrings', 'regexp'};

%% Default values for optional arguments
searchModeDefault = 'exact';        % search exact matches by default
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
addParameter(iP, 'SearchMode', searchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validSearchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, candidates, array, varargin{:});
searchMode = validatestring(iP.Results.SearchMode, validSearchModes);
ignoreCase = iP.Results.IgnoreCase;

%% Preparation
% Check type compatibility of candidates and array
% TODO: cellstr and string should be made compatible
% if ~strcmp(class(candidates), class(array))
%     error('candidates and array do not have the same type!!')
% end

%% Do the job
if ischar(candidates) || iscellstr(candidates) || isstring(candidates)
    % Make sure candidates is either a string or a cell array
    if ischar(candidates)
        candidates = {candidates};
    end

    % Use either cellfun or arrayfun
    if iscell(candidates)
        % Find the index in array for each element in candidates
        index = cellfun(@(x) find_in_strings(x, array, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), candidates);
    else
        % Find the index in array for each element in candidates
        index = arrayfun(@(x) find_in_strings(x, array, ...
                                'MaxNum', 1, 'ReturnNan', true, ...
                                'SearchMode', searchMode, ...
                                'IgnoreCase', ignoreCase), candidates);
    end
else
    % Find the index in array for each element in candidates
    %   Note: If not found, zero will be returned
    [~, index] = ismember(candidates, array);

    % Make all zeros NaNs instead
    index(index == 0) = NaN;
end

% Get all matched elements if requested
if nargout > 1
    matched = array(index);
end

%% Deal with outputs
varargout{1} = index;
if nargout > 1
    varargout{2} = matched;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
