function varExtracted = extract_vars (tableOrPath, varName, varargin)
%% Extracts variable(s) (column(s)) from a table
% Usage: varExtracted = extract_vars (tableOrPath, varName, varargin)
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
%       tableOrPath     - table or spreadsheet path
%                       must be a table 
%                           or a string scalar or a character vector
%       varName         - variable name to extract
%                       must be a string scalar or a character vector
%       varargin        - 'RowsToExtract': rows to extract
%                       must be a numeric array,
%                           a string scalar or a character vector, 
%                           or a cell array of character vectors
%                       default == 'all' (no restrictions)
%                       - 'RowConditions': conditions (name-value pairs) 
%                                           for row selection
%                       must be a scalar struct or a cell array
%                       default == [] (no restrictions)
%
% Requires:
%       cd/arglist2struct.m
%       cd/compute_combined_trace.m
%       cd/create_error_for_nargin.m
%       cd/ismatch.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting75.m

% File History:
% 2019-12-03 Moved from singleneuronfitting.m
% TODO: Output a cell array if extracting more than one variable
% TODO: Add 'OutputMode' to optionally extract as many outputs

%% Hard-coded parameters

%% Default values for optional arguments
rowsToExtractDefault = 'all';
rowConditionsDefault = [];

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
% TODO: Validation
addRequired(iP, 'tableOrPath');
addRequired(iP, 'varName');

% Add parameter-value pairs to the Input Parser
% TODO: Validation
addParameter(iP, 'RowsToExtract', rowsToExtractDefault);
addParameter(iP, 'RowConditions', rowConditionsDefault);

% Read from the Input Parser
parse(iP, tableOrPath, varName, varargin{:});
rowsToExtract = iP.Results.RowsToExtract;
rowConditions = iP.Results.RowConditions;

%% Do the job
% Parse the first argument
if istable(tableOrPath)
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

% Restrict to the rows of interest
if isnumeric(rowsToExtract)
    tableOfVars = tableOfVars(rowsToExtract, :);
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

% Extract variable
varExtracted = tableOfVars{:, varName};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%