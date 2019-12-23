function valueStructs = extract_param_values (paramTables, varargin)
%% Extracts parameter values as structure(s)
% Usage: valueStructs = extract_param_values (paramTables, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       valueStructs    - structure(s) where each field is a parameter
%                       specified as a structure array
%
% Arguments:
%       paramTables - parameter table(s) with parameter names as 
%                       row names and at least these variables:
%                           Value
%                   specified as a 2d table or a cell array of 2d tables
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/transpose_table.m
%
% Used by:
%       cd/m3ha_plot_figure03.m

% File History:
% 2019-12-23 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'paramTables', ...
    @(x) validateattributes(x, {'cell', 'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paramTables, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Do for each table
if iscell(paramTables)
    valueStructs = cellfun(@(x) extract_param_values_helper(x), paramTables);
else
    valueStructs = extract_param_values_helper(paramTables);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valueStruct = extract_param_values_helper (paramTable)

% Restrict to the Value column
paramTableRestricted = paramTable(:, 'Value');

% Transpose the table, with no new row names
paramTableTransposed = transpose_table(paramTableRestricted, ...
                                        'RowNames', 'suppress');

% Convert to a structure
valueStruct = table2struct(paramTableTransposed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%