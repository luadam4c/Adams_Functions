function [tableToPrint] = write_table(T, varargin)
%% A wrapper for writetable that converts cell array columns to character vectors
% Usage: [tableToPrint] = write_table(T, varargin)
% Explanation:
%       This function is a wrapper for the built-in writetable(). Before writing,
%       it identifies all columns in the input table T that are cell arrays and
%       converts their contents into character vectors using print_cellstr.m.
%       This is useful for creating CSV files where complex cell contents
%       need to be represented as a single string. The modified table is then
%       passed to writetable() with any additional arguments provided.
%
% Example(s):
%       % Create a sample table with a cell array column
%       testData = table({'S01'; 'S02'}, [101; 102], {{'data1', 1}; {'data2', 2}}, ...
%                        'VariableNames', {'SubjectID', 'Trial', 'Results'});
%       % Convert and write the table to a file
%       convertedTable = write_table(testData, 'test_output.csv');
%       % Display the converted table
%       disp(convertedTable);
%
% Outputs:
%       tableToPrint - The modified table where cell array columns have been
%                      converted to columns of character vectors.
%                   specified as a table
%
% Arguments:
%       T           - The input table to process and write.
%                   must be a table
%       varargin    - 'Delimiter': used to delimit separate entries
%                   must be a character array or string
%                   default == ''
%                   - 'ColDelimiter': used to delimit separate columns
%                   must be a character array or string
%                   default == ', '
%                   - 'RowDelimiter': used to delimit separate rows
%                   must be a character array or string
%                   default == '; '
%                   - 'Prefix': prefix to prepend to each element
%                   must be a character array or string
%                   default == ''
%                   - 'Suffix': suffix to append to each element
%                   must be a character array or string
%                   default == ''
%                   - Any other parameter-value pairs for writetable(),
%                     including the mandatory filename as the first argument.
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_parameter_value_pairs.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/virt_analyze_sniff_whisk.m

% File History:
% 2025-09-12 Created by Gemini, based on a user request and example code
%            from virt_analyze_sniff_whisk.m.

%% Hard-coded parameters

%% Default values for optional arguments
delimiterDefault = ' ';                 % default: separate by ' '
rowDelimiterDefault = '';               % default: set in print_cellstr.m
colDelimiterDefault = '';               % default: set in print_cellstr.m
prefixDefault = '';                     % default: set in print_cellstr.m
suffixDefault = '';                     % default: set in print_cellstr.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Use an inputParser to validate that the first argument is a table.
iP1 = inputParser;
iP1.FunctionName = mfilename;
addRequired(iP1, 'T', @istable);
parse(iP1, T);

% Remove any parameter-value pairs
[paramsList, inputList] = extract_parameter_value_pairs(varargin, 'OutputMode', 'cell');

% Use a second inputParser to validate parameters for this function
iP2 = inputParser;
iP2.FunctionName = mfilename;
iP2.KeepUnmatched = true;                      % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP2, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP2, 'RowDelimiter', rowDelimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP2, 'ColDelimiter', colDelimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP2, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP2, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP2, paramsList{:});
delimiter = iP2.Results.Delimiter;               % Examples: ',' '/'
rowDelimiter = iP2.Results.RowDelimiter;
colDelimiter = iP2.Results.ColDelimiter;
prefix = iP2.Results.Prefix;
suffix = iP2.Results.Suffix;

% Keep unmatched arguments for the writetable() function
otherArguments = struct2arglist(iP2.Unmatched);

%% Preparation
% Create a copy of the table to modify, so the original table is not changed.
tableToPrint = T;

% Get the names of all columns in the table
colNames = tableToPrint.Properties.VariableNames;

%% Do the job
% Iterate over each column by name
for i = 1:length(colNames)
    colName = colNames{i};
    
    % Access the data within the current column
    colData = tableToPrint.(colName);
    
    % Check if the column's data type is a cell array
    if iscell(colData)
        % If it is a cell array, apply print_cellstr to each element.
        % cellfun applies the provided function to each cell in the column.
        % The result is a new cell array where each cell contains a
        % character vector, which is suitable for writing to a text file.
        tableToPrint.(colName) = ...
            cellfun(@(cellElement) print_cellstr(cellElement, ...
                        'Delimiter', delimiter, 'RowDelimiter', rowDelimiter, ...
                        'ColDelimiter', colDelimiter, 'Prefix', prefix, ...
                        'Suffix', suffix, 'ToPrint', false, ...
                        'OmitQuotes', true, 'OmitBraces', true, ...
                        'OmitNewline', true), ...
                    colData, 'UniformOutput', false);
    end
end

%% Output results
% If there are arguments for writetable (e.g., a filename),
% call the original function with the newly modified table.
if ~isempty(varargin)
    writetable(tableToPrint, inputList{:}, otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
