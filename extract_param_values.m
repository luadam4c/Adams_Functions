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
%       varargin    - 'RowsToExtract' - row indices or row names in swpInfo 
%                                       to be used
%                   must be a positive integer vector, a string array 
%                       or a cell array of character vectors
%                   default == []
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ispositiveintegervector.m
%       cd/transpose_table.m
%
% Used by:
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_plot_figure03.m

% File History:
% 2019-12-23 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
rowsToExtractDefault = [];

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
addParameter(iP, 'RowsToExtract', rowsToExtractDefault, ...
    @(x) assert(ispositiveintegervector(x) || iscellstr(x) || isstring(x), ...
                ['strs5 must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));

% Read from the Input Parser
parse(iP, paramTables, varargin{:});
rowsToExtract = iP.Results.RowsToExtract;

%% Do the job
% Do for each table
if iscell(paramTables)
    valueStructs = cellfun(@(x) epv_helper(x, rowsToExtract), paramTables);
else
    valueStructs = epv_helper(paramTables, rowsToExtract);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valueStruct = epv_helper (paramTable, rowsToExtract)
%% Extracts parameter names and values as a structure

% Restrict to the requested rows
tableFiltered = paramTable(rowsToExtract, :);

% Restrict to the Value column
tableSelected = tableFiltered(:, 'Value');

% Transpose the table, with no new row names
tableTransposed = transpose_table(tableSelected, 'RowNames', 'suppress');

% Convert to a structure
valueStruct = table2struct(tableTransposed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%