function [outTables, outSheetPaths] = ...
                combine_variables_across_tables (inputs, varargin)
%% Combines measures across different tables
% Usage: [outTables, outSheetPaths] = ...
%               combine_variables_across_tables (inputs, varargin)
% Explanation:
%       TODO
% Example(s):
%       load_examples;
%       combine_variables_across_tables(myCellTable, 'KeyVariable', 'Key')
%       combine_variables_across_tables(myCellTable)
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
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - 'KeyVariable': variable used as the joining key
%                   default == row names
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/apply_over_cell.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/clc2_plot_measures.m

% File History:
% 2019-03-15 Created by Adam Lu
% 

%% Hard-coded parameters
keyVar = 'setNumber';

%% Default values for optional arguments
variableNamesDefault = {};  % plot all variables by default
outFolderDefault = pwd;

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
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, inputs, varargin{:});
varNames = iP.Results.VariableNames;
outFolder = iP.Results.OutFolder;

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

% Decide on table names
if iscellstr(inputs)
else
    % Use nInputs TODO
    create_labels_from_numbers();
end

% Decide on variable names or restrict table
if ~isempty(varNames)
    % Restrict to provided variable names
    inTables = cellfun(@(x) x(:, varNames), inTables, 'UniformOutput', false);
else
    % Use all variable names
    varNames = cellfun(@(x) x.Properties.VariableNames, inTables, ...
                        'UniformOutput', false);
end

% Count the number of variable names
nVars = numel(varNames);

% Create output spreadsheet paths
outSheetPaths = fullfile(outFolder, strcat(varNames, '_all.csv'));

%% Do the job
% TODO
outTables = cell(nVars, 1);
parfor iVar = 1:nVars
    % Get this variable
    varThis = varNames{iVar};

    % Extract the key column and the readout column from each inTable
    tablesToJoin = x(:, {keyVar, varThis}), inTables);

    % Rename the variable by the original table's name


    % Join the tables by the key column
    if 
    outTables{iVar} = apply_over_cell(@outerjoin, tablesToJoin, ...
                                        'MergeKeys', true);


end

%% Save results
cellfun(@(x, y) writetable(x, y), outTables, outSheetPaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%