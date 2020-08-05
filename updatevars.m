function myTable = updatevars (myTable, varValue)
%% Replace a variable in a table or add it if it doesn't exist
% Usage: myTable = updatevars (myTable, varValue)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       myTable     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       myTable     - TODO: Description of reqarg1
%                   must be a TODO
%       varValue    - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Used by:
%       cd/m3ha_simulate_population.m

% File History:
% 2020-08-04 Moved from m3ha_simulate_population.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'myTable');
addRequired(iP, 'varValue');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, myTable, varValue, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
varName = inputname(2);
if is_field(myTable, varName)
    myTable.(varName) = varValue;
else
    myTable = addvars(myTable, varValue, 'NewVariableNames', varName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%