function [outTables, outSheetPaths] = ...
                combine_variables_across_tables (inputs, varargin)
%% Combines measures across different tables
% Usage: [outTables, outSheetPaths] = ...
%               combine_variables_across_tables (inputs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       combine_variables_across_tables(myCellTable, 'Keys', 'Key')
%       combine_variables_across_tables(myCellTable, 'Keys', 'Var')
%       combine_variables_across_tables(myCellTable)
%       combine_variables_across_tables(myCellTable, 'VariableNames', {'Var'; 'Key'})
%       combine_variables_across_tables(myCellTable, 'VariableNames', {'Key', 'Var'})
%       combine_variables_across_tables(myCellTable, 'VariableNames', {'Var', 'Key'})
%       combine_variables_across_tables(myCellTable, 'VariableNames', {'Var', 'Key'; 'Key', 'Var'})
%
% Outputs:
%       outTables   - tables organized with each 
%                   specified as a cell array of tables
%       outSheetPaths   - file names for saved tables
%                   specified as a cell array of character vectors
% Arguments:
%       inputs      - tables organized with each measure as a column
%                       or file names containing those tables
%                   must be a cell array of tables
%                       or a cell array of character vectors
%       varargin    - 'VariableNames': variable (column) names of the table
%                                       to combine; 
%                       Note: each row results in a different table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables as separate tables
%                   - 'InputNames': names of the input tables
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == use distinct parts of the file names
%                               or Input1, Input2, ...
%                   - 'OmitVarName': whether to omit variable name 
%                                       in output table
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Keys': variable(s) used as the joining key(s)
%                   default == row names
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Prefix': prefix for output files
%                   must be a character array
%                   default == ''
%                   - 'SaveFlag': whether to save the output tables
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the outerjoin() function
%
% Requires:
%       cd/apply_over_cells.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_row_labels.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/renamevars.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_measures.m

% File History:
% 2019-03-17 Created by Adam Lu
% 2019-03-24 Now allows variableNames to be 2 dimensional
% TODO: Use parfor instead?

%% Hard-coded parameters

%% Default values for optional arguments
variableNamesDefault = {};  % plot all variables by default
keysDefault = 'Row';        % use the row name by default
inputNamesDefault = {};     % set later
omitVarNameDefault = false; % don't omit variable name by default
outFolderDefault = pwd;
prefixDefault = '';                     % default: prepend nothing
saveFlagDefault = false;    % don't save output tables by default

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
addRequired(iP, 'inputs', ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Keys', keysDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'InputNames', inputNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'OmitVarName', omitVarNameDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveFlag', saveFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, inputs, varargin{:});
varNames = iP.Results.VariableNames;
keys = iP.Results.Keys;
inNames = iP.Results.InputNames;
omitVarName = iP.Results.OmitVarName;
outFolder = iP.Results.OutFolder;
prefix = iP.Results.Prefix;
saveFlag = iP.Results.SaveFlag;

% Keep unmatched arguments for the outerjoin() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Count the number of input tables
nInputs = numel(inputs);

% Read in tables
if iscellstr(inputs)
    % If inputs are file names, use readtable
    inTables = cellfun(@readtable, inputs, 'UniformOutput', false);
else
    % Otherwise the inputs are tables
    inTables = inputs;
end

% Check if row labels are provided
if strcmp(keys, 'Row')
    % Get row labels from each table
    rowNamesAll = cellfun(@(x) x.Properties.RowNames, inTables, ...
                            'UniformOutput', false);

    % Supply row names if all are empty
    if all(isemptycell(rowNamesAll))
        inTables = cellfun(@create_row_labels, inTables, ...
                                'UniformOutput', false);
    end
end

% Decide on input names
if isempty(inNames)
    if iscellstr(inputs)
        % Use the distinct parts of the file names
        inNames = extract_fileparts(inputs, 'distinct');
    else
        % Use Input #
        inNames = create_labels_from_numbers(1:nInputs, 'Prefix', 'Input');
    end
end

% Decide on variable names to combine or restrict table
if ~isempty(varNames)
    % Combine the variables except for 'Row'
    allVars = union_vars_to_extract(keys, varNames);

    % Restrict to provided variable names
    inTables = cellfun(@(x) x(:, allVars), inTables, 'UniformOutput', false);
else
    % Collect all possible variable names
    varNamesAll = cellfun(@(x) x.Properties.VariableNames, inTables, ...
                        'UniformOutput', false);

    % Take the union over all possible variable names
    varNames = apply_over_cells(@union, varNamesAll);

    % Force as a column cell array
    varNames = force_column_vector(varNames);

    % Omit keys from the variable names
    varNames = setdiff(varNames, keys);
end

% Count the number of rows
nRows = size(varNames, 1);

% Create output spreadsheet paths
if saveFlag
    % Combine all variable names that will be used for each file
    %   as a single suffix
    % TODO: Use construct_suffix.m
    varNamesEachFile = arrayfun(@(x) strjoin(varNames(x, :), '_'), ...
                                transpose(1:nRows), 'UniformOutput', false);

    % Create output spreadsheet paths
    outSheetPaths = fullfile(outFolder, ...
                            strcat(prefix, '_', varNamesEachFile, '_all.csv'));
end

%% Do the job
outTables = arrayfun(@(x) combine_variable_across_tables(varNames(x, :), ...
                    inTables, inNames, keys, omitVarName, otherArguments), ...
                    transpose(1:nRows), 'UniformOutput', false);

%% Save results
if saveFlag
    cellfun(@(x, y) writetable(x, y), outTables, outSheetPaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outTable = combine_variable_across_tables(varThis, inTables, ...
                                    inNames, keys, omitVarName, otherArguments)
%% Combine all columns with this variable across tables

% Combine the variables except for 'Row'
allVars = union_vars_to_extract(keys, varThis);

% Count the number of columns to extract
if ischar(allVars)
    nCols = 1;
else
    nCols = numel(allVars);
end

% Extract the variables from each table as a table
tablesToJoin = cellfun(@(x) x(:, allVars), inTables, 'UniformOutput', false);

% Reorder so that the variable to combine is the last column
tablesToJoin = cellfun(@(x) movevars(x, varThis, 'After', nCols), ...
                        tablesToJoin, 'UniformOutput', false);

% Rename the variables to combine
if omitVarName
    % Rename the variable(s) to combine by the input names
    tablesToJoin = cellfun(@(x, y) renamevars(x, varThis, y), ...
                            tablesToJoin, inNames, 'UniformOutput', false);
else
    % Rename the variable(s) to combine by concatenating the input name
    tablesToJoin = cellfun(@(x, y) renamevars(x, varThis, ...
                                            strcat(varThis, '_', y)), ...
                            tablesToJoin, inNames, 'UniformOutput', false);
end

% Join the tables by the key(s)
outTable = apply_over_cells(@outerjoin, tablesToJoin, ...
                            'Keys', keys, 'MergeKeys', true, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function allVars = union_vars_to_extract(keys, vars)

% Force as column cell arrays
[keys, vars] = argfun(@(x) force_column_cell(x, 'ToLinearize', true), keys, vars);

% Combine the keys and variables
allVars = [keys; vars];

% Remove 'Row' from the list of variables
allVars = setdiff(allVars, 'Row');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

tic
% Count the number of variable names
nVars = numel(varNames);
outTables = cell(nVars, 1);
parfor iVar = 1:nVars
    % Get this variable
    varThis = varNames{iVar};

    % Combine all columns with this variable
    outTables{iVar} = ...
        combine_variable_across_tables(varThis, inTables, inNames, keys);
end
toc

% Force as a cell array if not already
[keys, varNames] = argfun(@force_column_cell, keys, varNames);

% Force as a column cell array of column vectors
% [keys, varNames] = argfun(@force_column_cell, keys, varNames);
[keys, varNames] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true, ...
                                    'ForceCellOutput', true), ...
            keys, varNames);


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
