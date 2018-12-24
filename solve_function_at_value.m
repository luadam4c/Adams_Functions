function xValue = solve_function_at_value (functionStr, varargin)
%% Solves x for f(x) at a specific value (default is zero)
% Usage: xValue = solve_function_at_value (functionStr, varargin)
% Explanation:
%       TODO
% Example(s):
%       solve_function_at_value('cos(x)')
%       solve_function_at_value('sin(x)', 1)
%       solve_function_at_value('sin(y)', 'VarStr', 'y')
%       solve_function_at_value('y + 4', 3, 'VarStr', 'y')
% Outputs:
%       xValue      - value of x
%                   specified as a numeric scalar
% Arguments:
%       functionStr     - string for the function expression
%                       must be a string scalar or a character vector
%       functionValue   - (opt) value to set f(x) at
%                       must be a numeric scalar
%                       default == 0
%       varargin    - 'VarStr': string used as a variable
%                   must be a string scalar or a character vector
%                   default == 'x'
%                   - Any other parameter-value pair for the solve() function
%
% Requires:
%       ~/Adams_Functions/create_error_for_nargin.m
%
% Used by:
%       cd/compute_peak_decay.m

% File History:
% 2018-12-24 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
fValueDefault = 0;
varStrDefault = 'x';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'functionStr', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'FunctionValue', fValueDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VarStr', varStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, functionStr, varargin{:});
fValue = iP.Results.FunctionValue;
varStr = iP.Results.VarStr;

% Keep unmatched arguments for the solve() function
otherArguments = iP.Unmatched;

%% Do the job
% Make varStr a symbolic variable
x = sym(varStr);

% Create a symbolic equation
%   Note: str2sym() is introduced in R2017b
eqToSolve = str2sym([functionStr, ' == ', num2str(fValue)]);

% Solve for x in the equation
% TODO: xValueSymbolic = solve(eqToSolve, x, otherArguments);
xValueSymbolic = solve(eqToSolve, x);

% Convert from symbolic to numeric
xValue = double(xValueSymbolic);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%