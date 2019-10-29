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
% Requires:
%       cd/is_contained_in.m
%
% Used by:    
%       cd/m3ha_neuron_create_initial_params.m

% File History:
% 2018-10-31 Created by Adam Lu
% 2018-11-14 Now checks if each parameter exists in rownames
% 2018-11-14 Now check bounds if 'UpperBound' and 'LowerBound' fields exist
% 

%% Hard-coded parameters
upperBoundStr = 'UpperBound';
lowerBoundStr = 'LowerBound';

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

%% Preparation
% Get all parameter names from the table
paramsTableRowNames = paramsTable.Properties.RowNames;

% Get all variable names from the table
paramsTableVariableNames = paramsTable.Properties.VariableNames;

%% Do the job
% Get all parameter names in a cell array
paramNames = fieldnames(allParamValuePairs);

% Check if all parameter names exist in the table
if ~is_contained_in(paramNames, paramsTableRowNames)
    fprintf('Original parameters table returned!\n');
    return
end

% Get all parameter values in a cell array
paramValuesCell = struct2cell(allParamValuePairs);

% If upper bounds and lower bounds exist, 
%   check whether the values are within bounds
if is_contained_in({upperBoundStr, lowerBoundStr}, ...
                    paramsTableVariableNames, 'SuppressOutput', true)
    % Convert to a numeric array if poss
    paramValues = cell2mat(paramValuesCell);

    % Extract the upper and lower bounds for each parameter
    upperBounds = paramsTable{paramNames, upperBoundStr};
    lowerBounds = paramsTable{paramNames, lowerBoundStr};

    % Check if all values are within bounds
    if ~check_within_bounds(paramValues, lowerBounds, upperBounds)
        fprintf('Original parameters table returned!\n');
        return
    end
end

% Replace all parameters with new values
paramsTable(paramNames, {'Value'}) = paramValuesCell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
