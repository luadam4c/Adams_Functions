function paramsTable = update_param_values (paramsTable, varargin)
%% Updates a parameters table with new values
% Usage: paramsTable = update_param_values (paramsTable, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       paramsTable - updated parameters table
%                   specified as a table
% Arguments:
%       paramsTable - a parameters table with the row names being the 
%                       parameter names and a 'Value' column
%                   must be a table
%       varargin    - name-value pairs that matches row names in the table
%                   must be a list of string-numeric pairs
%
% Used by:    
%       cd/m3ha_create_initial_neuronparams.m

% File History:
% 2018-10-31 Created by Adam Lu
% TODO: Check if parameter exist in rownames
% TODO: Check bounds if 'UpperBound' and 'LowerBound' fields exist
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Allow extraneous options: these will be parameters to be updated
iP.KeepUnmatched = true;

% Add required inputs to the Input Parser
addRequired(iP, 'paramsTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paramsTable, varargin{:});
% param1 = iP.Results.param1;

% Read the name-value pairs unmatched by the input parser
%   These will be name-value pairs for the parameters
allParamValuePairs = iP.Unmatched;

%% Do the job
% Get all parameter names in a cell array
paramNames = fieldnames(allParamValuePairs);

% Get all parameter values in a cell array
paramValues = struct2cell(allParamValuePairs);

% Replace all parameters with new values
paramsTable(paramNames, {'Value'}) = paramValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%