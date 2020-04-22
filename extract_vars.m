function varExtracted = extract_vars (tableOrPath, varNames, varargin)
%% Extracts variable(s) (column(s)) from a table
% Usage: varExtracted = extract_vars (tableOrPath, varNames, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       extract_vars(myTable1, 'Var')
%       extract_vars(myTable1, 'Var', 'RowsToExtract', [1, 3, 4])
%       extract_vars(myTable1, 'Var', 'RowConditions', {'Key', 'b'})
%
% Outputs:
%       varExtracted    - variable contents extracted
%                       specified as a column vector
%
% Arguments:
%       tableOrPath - table or spreadsheet path
%                   must be a table 
%                       or a string scalar or a character vector
%       varNames    - variable name(s) to extract
%                   must be a string scalar or a character vector
%                       or a cell array of character vectors
%       varargin    - 'RowsToExtract': rows to extract
%                   must be a numeric array,
%                       a string scalar or a character vector, 
%                       or a cell array of character vectors
%                   default == 'all' (no restrictions)
%                   - 'RowConditions': conditions (name-value pairs) 
%                                       for row selection
%                   must be a scalar struct or a cell array
%                   default == [] (no restrictions)
%                   - 'ForceCellOutput': whether to force output as a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/arglist2struct.m
%       cd/compute_combined_trace.m
%       cd/create_error_for_nargin.m
%       cd/ismatch.m
%       cd/ispositiveintegervector.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_simulate_population.m
%       cd/plot_table_parallel.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting75.m

% File History:
% 2019-12-03 Moved from singleneuronfitting.m
% 2019-12-23 Now outputs a cell array if extracting more than one variable
% 2020-04-19 Now allows the first argument to be a cell array
% TODO: Allow other strings for rowsToExtract, such as 'first', 'last', etc.
% TODO: Add 'OutputMode' to optionally extract as many outputs

%% Hard-coded parameters

%% Default values for optional arguments
rowsToExtractDefault = 'all';
rowConditionsDefault = [];
forceCellOutputDefault = false; % don't force output as a cell array by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tableOrPath', ...
    @(x) assert(iscell(x) || ischar(x) || isstring(x) || istable(x), ...
        ['tableOrPath must be a character array or a string array ', ...
            'or a table!']));
addRequired(iP, 'varNames', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['varNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
% TODO: Validation
addParameter(iP, 'RowsToExtract', rowsToExtractDefault, ...
    @(x) assert(ispositiveintegervector(x) || iscellstr(x) || isstring(x), ...
                ['RowsToExtract must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'RowConditions', rowConditionsDefault);
addParameter(iP, 'ForceCellOutput', forceCellOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, tableOrPath, varNames, varargin{:});
rowsToExtract = iP.Results.RowsToExtract;
rowConditions = iP.Results.RowConditions;
forceCellOutput = iP.Results.ForceCellOutput;

%% Do the job
% Parse the first argument
if iscell(tableOrPath)
    % Recursively apply to each table
    varExtracted = ...
        cellfun(@(x) extract_vars(x, varNames, varargin{:}), ...
                tableOrPath, 'UniformOutput', false);
    return;
elseif istable(tableOrPath)
    % First argument already a table
    tableOfVars = tableOrPath;
elseif isfile(tableOrPath)
    % Read in the table
    tableOfVars = readtable(tableOrPath);
elseif ischar(tableOrPath) || isstring(tableOrPath)
    fprintf('%s does not exist!', tableOrPath);
else
    error('Wrong first argument!');
end      

% Make sure the table is not empty
if isempty(tableOfVars)
    varExtracted = [];
    return
end

% Restrict to the rows of interest
if isnumeric(rowsToExtract)
    tableOfVars = tableOfVars(rowsToExtract, :);
elseif ischar(rowsToExtract) && ~strcmpi(rowsToExtract, 'all')
    error('Not implemented yet!');
elseif ~isempty(rowConditions)
    % Find rows to restrict to based on conditions if provided
    % Force rowConditions as a scalar structure
    if iscell(rowConditions)
        rowConditions = arglist2struct(rowConditions);
    end
    
    % List the conditional variable names
    conditionVarNames = fieldnames(rowConditions);
    
    % Detemine whether each row matches the condition, for each condition
    isMatchingRowEachCondition = ...
        cellfun(@(x) ismatch(tableOfVars.(x), rowConditions.(x)), ...
                            conditionVarNames, 'UniformOutput', false);
                        
    % Determine whether each row is to be used
    rowsToExtract = compute_combined_trace(isMatchingRowEachCondition, 'all');
    
    % Extract those rows
    tableOfVars = tableOfVars(rowsToExtract, :);
end

% Make sure the table is not empty
if isempty(tableOfVars)
    varExtracted = [];
    return
end

% Extract variable
if (iscellstr(varNames) || isstring(varNames)) && ...
        (numel(varNames) > 1 || forceCellOutput)
    if iscell(varNames)
        varExtracted = cellfun(@(x) tableOfVars{:, x}, varNames, ...
                                'UniformOutput', false);
    else
        varExtracted = arrayfun(@(x) tableOfVars{:, x}, varNames, ...
                                'UniformOutput', false);
    end
else
    varExtracted = tableOfVars{:, varNames};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
