function [combinedTable] = combine_tables (tablesToCombine, varargin)
%% Robustly combines a cell array of tables into a single table
% Usage: [combinedTable] = combine_tables (tablesToCombine, varargin)
% Explanation:
%       This function takes a cell array of MATLAB tables, which may have
%       different sets of columns, and combines them vertically into a single
%       table. It handles mismatches in columns by creating a union of all
%       column names and adding any missing columns to each table with
%       default `NaN` values before concatenation. This prevents errors that
%       occur when using `vertcat` on tables with different variables.
%
% Example(s):
%       t1 = table([1;2], 'VariableNames', {'A'});
%       t2 = table([3;4], 'VariableNames', {'B'});
%
%       % Example with a struct
%       vars_struct = struct();
%       vars_struct.rep = [1; 2];
%       vars_struct.id = {'x'; 'y'};
%       combined1 = combine_tables({t1, t2}, 'VarsToAdd', vars_struct);
%
%       % Example with a table
%       vars_table = table([10; 20], ["ID1"; "ID2"], 'VariableNames', {'rep', 'id'});
%       combined2 = combine_tables({t1, t2}, 'VarsToAdd', vars_table);
%
% Outputs:
%       combinedTable - The resulting single table with all data.
%                   specified as a MATLAB table
%
% Arguments:
%       tablesToCombine - A cell array where each element is a MATLAB table.
%                       must be a cell array of tables
%       varargin    - 'VarsToAdd': A structure or table used to add identifying columns.
%                   must be a structure or a table
%                   If a struct, each field name is a new variable.
%                   If a table, each column name is a new variable.
%                   If a struct, the value of each field must be a vector
%                   with a length equal to the number of non-empty tables.
%                   If a table, it must have a number of rows equal to
%                   the number of non-empty tables.
%                   For each table `i`, the value `VarsToAdd.fieldName(i)`
%                   (for a struct) or `VarsToAdd{i, 'VarName'}` (for a table)
%                   is repeated to form the new column.
%                   default == struct()
%
% Requires:
%       None
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore_multiple_reps.m

% File History:
% 2025-10-12 Created by Gemini based on logic from virt_moore_multiple_reps.m
% 2025-10-12 Added 'VarsToAdd' optional argument by Gemini
% 2025-10-12 Restored full function template structure.
% 2025-10-12 Modified by Gemini to allow 'VarsToAdd' to also be a table.

%% Hard-coded parameters
% None for this function

%% Default values for optional arguments
varsToAddDefault = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error('combine_tables:NotEnoughInputs', 'This function requires at least one input argument.');
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;

% Add required and optional arguments
addRequired(iP, 'tablesToCombine', ...
    @(x) iscell(x) && all(cellfun(@(x) istable(x) || isempty(x), x)));
addParameter(iP, 'VarsToAdd', varsToAddDefault, @(x) isstruct(x) || istable(x));

% Parse the inputs
parse(iP, tablesToCombine, varargin{:});
tablesToCombine = iP.Results.tablesToCombine;
varsToAdd = iP.Results.VarsToAdd;

%% Preparation
% Remove any empty entries from the cell array to simplify processing
tablesToCombine = tablesToCombine(~cellfun('isempty', tablesToCombine));

% If no non-empty tables are left, return an empty table
if isempty(tablesToCombine)
    combinedTable = table();
    return;
end

%% Do the job
% Determine if 'VarsToAdd' is a non-empty struct or table
isVarsToAddStruct = isstruct(varsToAdd) && ~isempty(fieldnames(varsToAdd));
isVarsToAddTable = istable(varsToAdd) && ~isempty(varsToAdd);

% Add variables from VarsToAdd if provided
if isVarsToAddStruct || isVarsToAddTable
    nTables = numel(tablesToCombine);

    % Get variable names and perform validation based on the type
    if isVarsToAddStruct
        varsToAddNames = fieldnames(varsToAdd);
        % Validate that each vector in varsToAdd has the correct length
        for iVar = 1:numel(varsToAddNames)
            varName = varsToAddNames{iVar};
            if numel(varsToAdd.(varName)) ~= nTables
                error('combine_tables:VarsToAddSizeMismatch', ...
                      'The length of the vector for struct field ''%s'' must match the number of non-empty tables.', varName);
            end
        end
    else % It's a table
        varsToAddNames = varsToAdd.Properties.VariableNames;
        % Validate that the number of rows in the table matches the number of tables
        if height(varsToAdd) ~= nTables
            error('combine_tables:VarsToAddSizeMismatch', ...
                  'The number of rows in the ''VarsToAdd'' table must match the number of non-empty tables.');
        end
    end

    % Add the new columns to each table
    for iTbl = 1:nTables
        currentTable = tablesToCombine{iTbl};
        currentTableVars = currentTable.Properties.VariableNames;
        nEntries = height(currentTable);

        % Check and add to the front in reverse order to maintain desired final order
        for iVar = numel(varsToAddNames):-1:1
            varName = varsToAddNames{iVar};
            if ~ismember(varName, currentTableVars)
                % Get the scalar value to be added, based on the input type
                if isVarsToAddStruct
                    scalarValue = varsToAdd.(varName)(iTbl);
                else % It's a table
                    scalarValue = varsToAdd{iTbl, varName};
                end

                % Create the new column by repeating the scalar value
                if iscell(scalarValue) || isstring(scalarValue)
                    newColumn = repmat(scalarValue, nEntries, 1);
                else
                    newColumn = repelem(scalarValue, nEntries, 1);
                end
                currentTable = addvars(currentTable, newColumn, ...
                                    'Before', 1, 'NewVariableNames', varName);
            end
        end
        tablesToCombine{iTbl} = currentTable;
    end
end


% 1. Get the union of all variable names across all tables, maintaining order
allVarNames = {};
for iTbl = 1:numel(tablesToCombine)
    allVarNames = union(allVarNames, tablesToCombine{iTbl}.Properties.VariableNames, 'stable');
end

% If no variables were found at all, return an empty table
if isempty(allVarNames)
    combinedTable = table();
    return;
end

% 2. Standardize all tables to have the same columns
for iTbl = 1:numel(tablesToCombine)
    currentTable = tablesToCombine{iTbl};

    missingVars = setdiff(allVarNames, currentTable.Properties.VariableNames);
    for iVar = 1:numel(missingVars)
        % Add missing variables as columns of NaN. This is a general approach.
        % A more sophisticated version might try to infer data types.
        currentTable.(missingVars{iVar}) = nan(height(currentTable), 1);
    end

    % Reorder columns of the current table to match the master order
    tablesToCombine{iTbl} = currentTable(:, allVarNames);
end

% 3. Concatenate all the now-standardized tables
combinedTable = vertcat(tablesToCombine{:});

%% Output results
% The result `combinedTable` is returned directly by the function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%