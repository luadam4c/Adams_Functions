function tableNew = transpose_table (tableOld, varargin)
%% Transposes a table (make row names variable names and vice versa)
% Usage: tableNew = transpose_table (tableOld, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       tableNew = transpose_table(tableOld)
%       tableNew = transpose_table(tableOld, 'RowNames', 'suppress')
%
% Outputs:
%       tableNew    - new table
%                   specified as a 2d table or a cell array of 2d tables
%
% Arguments:    
%       tableOld    - old table
%                   must be a 2d table or a cell array of 2d tables
%       varargin    - 'RowNames': new row names, or 'suppress'
%                   must be a character vector 
%                       or a cell array of character vectors
%                   default == old variable names, or 'suppress'
%                   - 'VariableNames': new variable names
%                   must be a character vector 
%                       or a cell array of character vectors
%                   default == old row names
%
% Requires:
%       cd/create_labels_from_numbers.m
%
% Used by:
%       cd/extract_param_values.m
%       cd/m3ha_neuron_create_simulation_params.m

% File History:
% 2018-10-22 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
rowNamesDefault = {};
variableNamesDefault = {};

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
addRequired(iP, 'tableOld', ...
    @(x) validateattributes(x, {'table', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RowNames', rowNamesDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) ischar(x) || iscellstr(x) || isstring(x));

% Read from the Input Parser
parse(iP, tableOld, varargin{:});
rowNames = iP.Results.RowNames;
variableNames = iP.Results.VariableNames;

%% Do the job
if iscell(tableOld)
    tableNew = ...
        cellfun(@(x) transpose_table_helper(x, rowNames, variableNames), ...
                tableOld, 'UniformOutput', false);
else
    tableNew = transpose_table_helper(tableOld, rowNames, variableNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tableNew = transpose_table_helper (tableOld, rowNames, variableNames)

% Get just the values
tableValuesOld = table2array(tableOld);

% Transpose the table values
tableValuesNew = transpose(tableValuesOld);

% Construct new table
tableNew = array2table(tableValuesNew);

% Decide on the new variable names
if ischar(variableNames) && strcmpi(variableNames, 'suppress')
    % Do nothing
elseif ~isempty(variableNames)
    % Use the user provided variable names
    tableNew.Properties.VariableNames = variableNames;
else
    if ~isempty(tableOld.Properties.RowNames)
        % Use the old row names
        tableNew.Properties.VariableNames = tableOld.Properties.RowNames;
    else
        % Count the number of rows in the old table
        nRowsOld = height(tableOld);

        % Create old row names
        tableNew.Properties.VariableNames = ...
            create_labels_from_numbers(1:nRowsOld, 'Prefix', 'oldRow');
    end
end

% Decide on the new row names
if ischar(rowNames) && strcmpi(rowNames, 'suppress')
    % Do nothing
elseif ~isempty(rowNames)
    % Use the user provided row names
    tableNew.Properties.RowNames = rowNames;
else
    % Use the old variable names
    tableNew.Properties.RowNames = tableOld.Properties.VariableNames;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% The following is slower than cellfun
% Count the number of tables
nTables = numel(tableOld);

% Transpose each table
tableNew = cell(size(tableOld));
parfor iTable = 1:nTables
    tableNew{iTable} = transpose_table_helper(tableOld{iTable});
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
