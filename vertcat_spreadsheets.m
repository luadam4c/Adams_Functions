function varargout = vertcat_spreadsheets (inputTablesOrFileNames, varargin)
%% Combine spreadsheets using readtable, vertcat, then writetable
% Usage: outputTable = vertcat_spreadsheets (inputTablesOrFileNames, varargin)
% Explanation:
%       TODO
%
% Example:
%       TODO
%
% Side Effects:
%       TODO
%
% Outputs:
%       outputTable - combined table
%                   specified as a 2-D table
%
% Arguments:    
%       inputTablesOrFileNames  - input tables or spreadsheet file names
%                       must be a cell array of tables
%                           or a cell array of character arrays
%       varargin    - 'OutputFileName': output spreadsheet file name
%                   must be a string scalar or a character vector
%                   default == ''
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%
% Used by:
%       cd/combine_swd_sheets.m
%       cd/combine_swd_resp_data.m
%       cd/m3ha_simulate_population.m
%       cd/parse_all_swds.m
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%
% File History:
% 2018-05-15 Created by Adam Lu
% 2018-12-26 Turned 'OutputFileName' into an optional argument
% 2020-07-23 Now uses varargout
% 2020-08-19 Now allows the argument to be a cell array of tables
% TODO: Use apply_over_cells.m

%% Hard-coded parameters

%% Default values for optional arguments
outputFileNameDefault = '';

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
addRequired(iP, 'inputTablesOrFileNames', ...
    @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x) || istable(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputFileName', outputFileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, inputTablesOrFileNames, varargin{:});
outputFileName = iP.Results.OutputFileName;

%% Prepare
% Force file names as a cell array
inputTablesOrFileNames = force_column_cell(inputTablesOrFileNames);

% Get the number of spreadsheets to combine
nTables = numel(inputTablesOrFileNames);

%% Read and concatenate
outputTable = table.empty;
for iTable = 1:nTables
    % Read in the table for this input file
    if istable(inputTablesOrFileNames{iTable})
        inputTable = inputTablesOrFileNames{iTable};
    else
        inputTable = readtable(inputTablesOrFileNames{iTable});
    end

    % Vertically concatenate the new input table to the existing outputTable
    outputTable = vertcat(outputTable, inputTable);
end

%% Save output
if ~isempty(outputFileName)
    writetable(outputTable, outputFileName);
end

%% Outputs
varargout{1} = outputTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% The following will reorder the rows, so is not ideal
outputTable = outerjoin(outputTable, inputTable, 'MergeKeys', true);

outputTable = [];
% Concatenate the input table with the output table
if isempty(outputTable)
    % The inputTable is the first one; assign it as the outputTable
    outputTable = inputTable;
else
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
