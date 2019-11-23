function [results, validatedTypes] = istype (strs, validTypes, varargin)
%% Check whether a string or each string in a cell array is a valid type specified by validTypes
% Usage: [results, validatedTypes] = istype (strs, validTypes, varargin)
% Outputs:    
%       results     - indication of whether the string is of a specified type
%                   specified as a logical array
%       validatedTypes  - validated types, if any
%                   specified as a string/char-vec or 
%                       a cell array of strings/char-vecs
%                   returns the shortest match if matchMode == 'substring' 
%                       (sames as validatestring())
% Arguments:
%       strs        - string or strings to check
%                   must be a string/char-vec or 
%                       a cell array of strings/char-vecs
%       varargin    - 'ValidateMode': whether to validate string and 
%                       throw error if string is not a substring of a type
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
%
% Requires:
%       cd/isemptycell.m
%       cd/validate_string.m
%
% Used by:
%       cd/isfigtype.m
%       cd/islinestyle.m
%       cd/issheettype.m
%
% File History:
% 2018-05-16 Modified from issheettype.m
% 

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

% Add required inputs to an Input Parser
addRequired(iP, 'strs', ...                     % string or strings to check
    @(x) assert(ischar(x) || ...
                iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x) , ...
                ['strs must be either a string/character array ', ...
                'or a cell array of strings/character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValidateMode', false, ...     % whether to validate string
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchMode', 'substring', ...  % the matching mode
    @(x) any(validatestring(x, {'exact', 'substring'})));

% Read from the Input Parser
parse(iP, strs, varargin{:});
validateMode = iP.Results.ValidateMode;
matchMode = iP.Results.MatchMode;

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

%% Check strings and validate
if iscell(strs)
    validatedTypes = cellfun(@(x) validate_string(x, validTypes, ...
                                            'ValidateMode', validateMode, ...
                                            'MatchMode', matchMode), ...
                                            strs, 'UniformOutput', false);
    results = ~isemptycell(validatedTypes);
elseif ischar(strs)
    validatedTypes = validate_string(strs, validTypes, ...
                                'ValidateMode', validateMode, ...
                                'MatchMode', matchMode);
    results = ~isempty(validatedTypes);
else
    error(['input argument is in the wrong format! ', ...
            'Type ''help %s'' for usage'], mfilename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%