function [strValidated, isValid] = validate_string (str, validStrings, varargin)
%% Validate whether a string is an element of a cell array of valid strings
% Usage: [strValidated, isValid] = validate_string (str, validStrings, varargin)
% Explanation:
%       Same as the built-in validatestring() except that an empty string can be 
%       returned if str is not matched in validStrings, and an empty string
%       can be valid if validStrings contain empty
%
% Outputs:
%       strValidated    - validated text, i.e., the shortest match in  
%                           validStrings that contains str as a substring
%                           or an empty string if no match found
% Arguments:
%       str             - text to validate
%                       must be a character vector or a string scalar
%       validStrings    - text to match
%                       must be a cell array of character vectors 
%                           or a string array
%       varargin    - 'ValidateMode': whether to throw error if not a match
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following:
%                       'exact'         - string must be exact
%                       'substring'     - string can be a substring
%                   if 'ValidateMode' is 'true', matching mode is 
%                       automatically 'substring'
%                   default == 'substring'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/find_in_strings.m
%       cd/isemptycell.m
%
% Used by:
%       cd/istype.m
%       cd/minEASE_gui_examine_events.m

% File History:
% 2017-05-25 Moved from isfigtype.m
% 2017-05-25 Made everything more general
% 2017-05-25 Added input parser scheme
% 2017-06-09 Added IgnoreCase and made true the default
% 2018-02-01 Updated description
% 2018-05-15 issheettype.m also uses this now
% 2018-05-16 Now is used by istype.m instead
% 2025-09-02 Now returns isValid and make it true if an empty string
%               is validated in a list of strings including an empty string
%

%% Default values for optional arguments
validateModeDefault = false;            % whether to throw error if not a match
matchModeDefault = 'substring';         % default matching mode
ignoreCaseDefault = true;               % whether to ignore case by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'validate_string';

% Add required inputs to an Input Parser
addRequired(iP, 'str', ...                      % text to validate
    @(x) validateattributes(x, {'char', 'string'}, {}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'validStrings', ...             % text to match
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
                || min(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValidateMode', validateModeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', matchModeDefault, ...     % matching mode
    @(x) any(validatestring(x, {'exact', 'substring'})));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, str, validStrings, varargin{:});
validateMode = iP.Results.ValidateMode;
matchMode = iP.Results.MatchMode;
ignoreCase = iP.Results.IgnoreCase;

% Deal with special case: str is empty
if isempty(str)
    % The validated string is still empty
    strValidated = str;

    % It's valid if the validStrings contain empty
    if any(isemptycell(strValidated))
        isValid = true;
    else
        isValid = false;
    end
    return
end

% Validate string
if validateMode     % throws error if no match found
    % Find the shortest match in validStrings that contains str as a substring 
    %   and throw error if not found
    strValidated = validatestring(str, validStrings);
else                % returns an empty string if no match found
    % Find all possible matches that contains string according to matchMode
    [~, strValidated] = ...
        find_in_strings(str, validStrings, ...
                             'SearchMode', matchMode, ...
                             'IgnoreCase', ignoreCase);

    % If there is more than one match, 
    %   return the one with shortest character length
    if iscell(strValidated)
        % Find the character lengths of all returned figure types
        charlengths = cellfun(@length, strValidated);    

        % Find the index with minimum character length
        [~, indx] = min(charlengths);

        % Return the figure type with minimum character length
        strValidated = strValidated{indx};
    end
end

% Return whether string is validated
if isempty(strValidated)
    isValid = false;
else
    isValid = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
