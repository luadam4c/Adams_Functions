function [popTable] = compute_population_average (inTable, varargin)
%% Computes the population mean and confidence intervals from a table or time table
% Usage: [popTable] = compute_population_average (inTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       popTable    - population table
%                   specified as a table
%
% Arguments:
%       inTable     - input table
%                   must be a table
%       varargin    - 'VarStr': string for variable of interest
%                   must be a string scalar or a character vector
%                   default == '' (extract all columns)
%                   - 'SheetName' - spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == '' (don't save)
%
% Requires:
%       cd/combine_phase_numbers.m
%       cd/compute_stats.m
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/write_timetable.m
%
% Used by:
%       cd/plot_measures.m

% File History:
% 2019-08-27 Moved from plot_measures.m
% 2019-08-27 Now accepts tables or timetables
% 2019-08-27 Made 'VarStr' and 'SheetName' optional arguments
% 

%% Hard-coded parameters
popVarNames = {'mean', 'lower95', 'upper95'};
phaseNumStr = 'phaseNumber';

%% Default parameters for optional arguments
varStrDefault = '';
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

% Add required inputs to the Input Parser
addRequired(iP, 'inTable');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'VarStr', varStrDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, inTable, varargin{:});
varStr = iP.Results.varStr;
sheetName = iP.Results.SheetName;

%% Preparation
% Decide whether input is a time table
if istimetable(inTable)
    isTimeTable = true;
else
    isTimeTable = false;
end

% Prepend the variable name to the prefix
if ~isempty(varStr)
    popVarNamesToUse = strcat(varStr, '_', popVarNames);
else
    popVarNamesToUse = popVarNames;
end

%% Do the job
% Get all the variable names
columnNames = inTable.Properties.VariableNames;

% Extract the data
if ~isempty(varStr)
    % Select the columns that contains a given variable string
    columnsToAverage = contains(columnNames, varStr);

    % Extract the data from these columns
    popData = inTable{:, columnsToAverage};
else
    % Extract everything
    popData = inTable{:, :};
end

% Select the columns that contains a given variable string
phaseNumberColumns = contains(columnNames, phaseNumStr);

% Extract the data from these columns
phaseNumbersAll = inTable{:, phaseNumberColumns};

% Try combining the phase numbers
phaseNumbers = combine_phase_numbers(phaseNumbersAll);

% Only output a phaseNumber column if there are no conflicts
if ~isempty(phaseNumbers)

end

% Compute the means and bounds of 95% confidence intervals
[means, lower95s, upper95s] = ...
    argfun(@(x) compute_stats(popData, x, 2, 'IgnoreNan', true), ...
            'mean', 'lower95', 'upper95');

if isTimeTable
    % Save the row times from inTable
    rowTimes = inTable.Properties.RowTimes;

    % Create output table
    popTable = timetable(means, lower95s, upper95s, ...
                        'VariableNames', popVarNamesToUse, ...
                        'RowTimes', rowTimes);
else
    % Create output table
    popTable = table(means, lower95s, upper95s, ...
                        'VariableNames', popVarNamesToUse);

    % Use the same row names as before
    popTable.Properties.RowNames = inTable.Properties.RowNames;
end

%% Save the table
if ~isempty(sheetPath)
    if isTimeTable
        write_timetable(popTable, sheetPath);
    else
        writetable(popTable, sheetPath);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%