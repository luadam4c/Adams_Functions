function table = modify_table (table, varargin)
%% Modify (possibly only parts of) a table by applying a specific function
% Usage: table = modify_table (table, myFunction (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples; 
%       modify_table(myTable1, @movingaveragefilter, 'VariableNames', 'Var')
%
% Outputs:
%       table    - output table
%                   specified as a table
%
% Arguments:
%       table     - input table
%                   must be a table
%       myFunction  - (opt) a custom function
%                   must be a function handle
%                   default == none
%       varargin    - 'VariableNames': variable (column) names of the table
%                                       to modify;
%                       Note: all other variables are returned as is
%                               doesn't have to be exact match by default
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == all variables
%                   - 'UseVarFun': whether to apply the function to each column
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SheetName' - spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == '' (don't save)
%                   - Any other parameter-value pair for myFunction()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/plot_measures.m

% File History:
% 2019-08-28 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
myFunctionDefault = function_handle.empty;
variableNamesDefault = {};  % plot all variables by default
useVarFunDefault = [];
sheetNameDefault = '';

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
addRequired(iP, 'table', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'myFunction', myFunctionDefault, ...
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VariableNames', variableNamesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VariableNames must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'UseVarFun', useVarFunDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, table, varargin{:});
myFunction = iP.Results.myFunction;
varNames = iP.Results.VariableNames;
useVarFun = iP.Results.UseVarFun;
sheetName = iP.Results.SheetName;

% Keep unmatched arguments for the myFunction() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% If no function provided, just return table as is
if isempty(myFunction)
    return
end

% Get all the variable names
columnNames = table.Properties.VariableNames;

% Decide on variables to restrict to
if isempty(varNames)
    varNames = columnNames;
else
    % Restrict to certain variables if requested
    varNames = columnNames(contains(columnNames, varNames));
end

%% Do the job
if useVarFun
    % Display warning if necessary
    if ~isempty(otherArguments)
        disp('otherArguments ignored when useVarFun is used!');
    end

    % Extract the table with just the variables
    tableOfInterest = table(:, varNames);

    % Apply function to each column separately
    tableModified = varfun(myFunction, tableOfInterest);

    % Place back in the table
    table(:, varNames) = tableModified;
else
    % Extract the data from just the variables
    dataExtracted = table{:, varNames};

    % Apply the myFunction to dataExtracted
    dataModified = myFunction(dataExtracted, otherArguments{:});

    % Place it back in the table
    table{:, varNames} = dataModified;
end

%% Save table
if ~isempty(sheetName)
    writetable(table, sheetName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%