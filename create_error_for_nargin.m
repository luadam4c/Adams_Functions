function errorText = create_error_for_nargin (functionName, varargin)
%% Creates an error text for not having enough input arguments
% Usage: errorText = create_error_for_nargin (functionName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       errorText   - documentation text
%                   specified as a character array
% Arguments:
%       functionName    - the function name to look for documentation for
%                       must be a string scalar or a character vector
%       varargin    - 'ToPrint': whether to print the error text
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ErrorType': type of error
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'TooFew'    - there are too few input arguments
%                       'TooMany'   - there are too many input arguments
%                   default == 'TooFew'
%
% Requires:
%       cd/print_help.m
%
% Used by:
%       cd/print_help.m
%       ~/Settings_Matlab/function_template.m
%       ~/Settings_Matlab/function_template_simple.m
%       any .m file with a number of arguments check

% File History:
% 2018-12-17 Created by Adam Lu
% 

%% Hard-coded parameters
validStrings = {'TooFew', 'TooMany'};
firstStringTooFew = 'Not enough input arguments!!';
firstStringTooMany = 'Too many input arguments!!';

%% Default values for optional arguments
toPrintDefault = true;          % default: print the string
errorTypeDefault  = 'TooFew';   % default TODO: Description of errorType

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'functionName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ToPrint', toPrintDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(iP, 'ErrorType', errorTypeDefault, ...
    @(x) any(validatestring(x, validStrings)));

% Read from the Input Parser
parse(iP, functionName, varargin{:});
toPrint = iP.Results.ToPrint;
errorType = validatestring(iP.Results.ErrorType, validStrings);

%% Preparation
% Decide on the first string
switch errorType
    case 'TooFew'
        firstString = firstStringTooFew;
    case 'TooMany'
        firstString = firstStringTooMany;
    otherwise
        error('Code logic error!');
end

%% Do the job
% Create the help text
helpText = print_help(functionName, 'ToPrint', false);

% Create the error text
errorText = sprintf('%s\n%s:\n%s', firstString, functionName, helpText);

% Display the error text if requested
if toPrint
    disp(errorText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%