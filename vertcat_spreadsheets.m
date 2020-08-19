function varargout = vertcat_spreadsheets (inputFileNames, varargin)
%% Combine spreadsheets using readtable, vertcat, then writetable
% Usage: outputTable = vertcat_spreadsheets (inputFileNames, varargin)
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
%       inputFileNames  - input spreadsheet file names
%                       must be a cell array of strings or character arrays
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
%       cd/combine_swd_pleth_data.m
%       cd/m3ha_simulate_population.m
%       cd/parse_all_swds.m
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%
% File History:
% 2018-05-15 Created by Adam Lu
% 2018-12-26 Turned 'OutputFileName' into an optional argument
% 2020-07-23 Now uses varargout
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
addRequired(iP, 'inputFileNames', ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputFileName', outputFileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, inputFileNames, varargin{:});
outputFileName = iP.Results.OutputFileName;

%% Prepare
% Force file names as a cell array
inputFileNames = force_column_cell(inputFileNames);

% Get the number of spreadsheets to combine
nTables = numel(inputFileNames);

%% Read and concatenate
outputTable = table.empty;
for iTable = 1:nTables
    % Read in the table for this input file
    inputTable = readtable(inputFileNames{iTable});

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
