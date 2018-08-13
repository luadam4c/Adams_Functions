function vertcat_spreadsheets (inputFileNames, outputFileName, varargin)
%% Combine spreadsheets using readtable, vertcat, then writetable
% Usage: vertcat_spreadsheets (inputFileNames, outputFileName, varargin)
% Explanation:
%       TODO
% Example:
%       TODO
% Side Effects:
%       TODO
% Arguments:    
%       inputFileNames  - input spreadsheet file names
%                       must be a cell array of strings or character arrays
%       outputFileName  - output spreadsheet file name
%                       must be a string scalar or a character vector
%
% Used by:    
%       /home/Matlab/EEG_gui/combine_EEG_gui_outputs.m
%
% File History:
% 2018-05-15 Created by Adam Lu
% 

%% Hard-coded parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'inputFileNames', ...
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
                || min(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));
%    @(x) iscellstr(x) || isstring(x));
addRequired(iP, 'outputFileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, inputFileNames, outputFileName, varargin{:});

%% Prepare
% Get the number of spreadsheets to combine
nTables = numel(inputFileNames);

%% Read and concatenate
outputTable = [];
for iTable = 1:nTables
    % Read in the table for this input file
    inputTable = readtable(inputFileNames{iTable});

    % Concatenate the input table with the output table
    if isempty(outputTable)
        % The inputTable is the first one; assign it as the outputTable
        outputTable = inputTable;
    else
        % Vertically concatenate the new input table to the existing outputTable
        outputTable = vertcat(outputTable, inputTable);
    end
end

%% Save output
writetable(outputTable, outputFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% The following will reorder the rows, so is not ideal
outputTable = outerjoin(outputTable, inputTable, 'MergeKeys', true);

%}
