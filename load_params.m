function paramTables = load_params (fileNames, varargin)
%% Loads parameters from file(s) into a table
% Usage: paramTables = load_params (fileNames, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       paramTables - parameter table(s) with parameter names as 
%                       row names and at least these variables:
%                           Value
%                           LowerBound
%                           UpperBound
%                   specified as a 2d table or a cell array of 2d tables
%
% Arguments:    
%       fileNames   - parameter file name(s)
%                       Note: If it's a spreadsheet file, 
%                               the first column must be parameter names
%                   must be a string array or a character vector
%                       or a cell array of character vectors
%
% Requires:
%       cd/construct_and_check_fullpath.m
%       cd/force_logical.m
%       cd/issheettype.m
%       cd/renamevars_custom.m
%
% Used by:
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_pfiles2csv.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_xolotl_create_neuron.m

% File History:
% 2018-10-16 Created by Adam Lu
% 2018-10-21 Now returns a table
% 2018-10-22 Now reads in .p files too
% 2018-10-22 Now standardizes the variable names
% 2019-12-04 Added default extension

%% Hard-coded parameters
defaultExtension = 'csv';

% TODO: Make optional argument
directory = '';

%% Default values for optional arguments

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
addRequired(iP, 'fileNames', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, fileNames, varargin{:});

%% Preparation
% Construct full paths and check whether the files exist
%   TODO: Expand to accept optional 'Directory', 'Suffix', etc.
[fullPaths, pathExists] = ...
    construct_and_check_fullpath(fileNames, 'Directory', directory);

% If no path exists, try adding the default extension to the end of all paths
if ~all(pathExists)
    [fullPaths, pathExists] = ...
        construct_and_check_fullpath(fileNames, 'Directory', directory, ...
                                        'Extension', defaultExtension);
end

% Return if not all paths exist
if ~all(pathExists)
    paramTables = [];
    return
end

%% Do the job
if iscell(fullPaths)
    paramTables = cellfun(@load_params_helper, fullPaths, ...
                            'UniformOutput', false);
else
    paramTables = load_params_helper(fullPaths);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramTable = load_params_helper (fullPath)
%% Loads a file into a table

% Get the file extension
[~, ~, fileExt] = fileparts(fullPath);

% Read the file into a table
if issheettype(fileExt)
    % Read the spreadsheet file into a table with the first column
    %   being the row names and the first column being the variable names 
    %   (Name, Value, LowerBound, UpperBound, etc.)
    paramTable = readtable(fullPath, ...
                            'ReadRowNames', true, ...
                            'ReadVariableNames', true);

    % Get all variable names
    variableNames = paramTable.Properties.VariableNames;

    % Replace the variable names if needed
    paramTable = renamevars_custom(paramTable, {'val', 'min', 'max'}, ...
                            {'Value', 'LowerBound', 'UpperBound'}, ...
                            'SearchMode', 'substrings', 'IgnoreCase', true, ...
                            'MaxNum', 1);
else
    % Import the data as a cell array
    paramsCell = importdata(fullPath);

    % Set the first column to be row names
    rowNames = paramsCell(:, 1);

    % Get the rest of the cell array
    paramsCellRest = paramsCell(:, 2:end);

    % Count the number of variables
    nVariables = size(paramsCellRest, 2);

    % Construct the header
    header = cell(nVariables, 1);
    for iVariable = 1:nVariables
        if iVariable == 1
            header{iVariable} = 'Value';
        elseif iVariable == 2
            header{iVariable} = 'LowerBound';
        elseif iVariable == 3
            header{iVariable} = 'UpperBound';
        else
            header{iVariable} = strcat('Column', num2str(iVariable));
        end
    end

    % Convert the data to a table
    paramTable = cell2table(paramsCellRest, ...
                            'RowNames', rowNames, 'VariableNames', header);
end

% Convert any binary numeric column to a logical column
paramTableTemp = varfun(@force_logical, paramTable);
paramTableTemp.Properties.VariableNames = paramTable.Properties.VariableNames;
paramTableTemp.Properties.RowNames = paramTable.Properties.RowNames;
paramTable = paramTableTemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       paramNames  - parameter names
%                   specified as a column cell vector of character vectors
%       paramValues - parameter values
%                   specified as a column numeric vector
%       paramLBs    - parameter lower bounds
%                   specified as a column numeric vector
%       paramUBs    - parameter upper bounds
%                   specified as a column numeric vector

% Read from the table
paramNames = neuronParamsTable.Name;
paramValues = neuronParamsTable.Value;
paramLBs = neuronParamsTable.LowerBound;
paramUBs = neuronParamsTable.UpperBound;

% Assume the data is a MATLAB data file
mFile = matfile(fileName);

% Get all the variable names
varNames = who(mFile);

% Read in the first variable
paramsCell = mFile.(varNames{1});

paramTable = cell(size(fullPath));
parfor iFile = 1:numel(fullPath)
    paramTable{iFile} = load_params_helper(fullPath{iFile});
end

paramTable = renamevars_custom(paramTable, 'val', 'Value', ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);
paramTable = renamevars_custom(paramTable, 'min', 'LowerBound', ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);
paramTable = renamevars_custom(paramTable, 'max', 'UpperBound', ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

idxVal = find_in_strings('val', variableNames, ...
                              'SearchMode', 'substrings', ...
                              'IgnoreCase', true, 'MaxNum', 1);
if ~isempty(idxVal)
    paramTable.Properties.VariableNames{idxVal} = 'Value';
end
idxMin = find_in_strings('min', variableNames, ...
                              'SearchMode', 'substrings', ...
                              'IgnoreCase', true, 'MaxNum', 1);
if ~isempty(idxMin)
    paramTable.Properties.VariableNames{idxMin} = 'LowerBound';
end
idxMax = find_in_strings('max', variableNames, ...
                              'SearchMode', 'substrings', ...
                              'IgnoreCase', true, 'MaxNum', 1);
if ~isempty(idxMax)
    paramTable.Properties.VariableNames{idxMax} = 'UpperBound';
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
