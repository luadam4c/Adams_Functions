function [cellIdsSelected, cellInfo, swpInfo] = m3ha_select_cells (varargin)
%% Selects cells with sweeps to use for all pharm-gIncr pairs
% Usage: [cellIdsSelected, cellInfo, swpInfo] = m3ha_select_cells (varargin)
% Explanation: 
%   Finds all cell IDs and cell names (from 'CellInfo') for the rows 
%       restricted by 'ToUse' in 'SwpInfo'
%
% Outputs:
%       TODO
%
% Arguments:
%       varargin    - 'DataMode': data mode
%                   must be a one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                       3 - all data 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 2
%                   - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                       toUse       - whether the sweep is to be used (optional)
%                   default == m3ha_load_sweep_info
%                   - 'CellInfo': cell name info
%                   must a 2D table with row indices being cell IDs 
%                       and with fields:
%                       cellName - cell names
%                   default == detected from swpInfo
%                   - 'RowsToUse' - row indices or row names in swpInfo 
%                                       of sweeps to fit
%                   must be a positive integer vector, a string array 
%                       or a cell array of character vectors
%                   default == []
%                   - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == ~/m3ha/data_dclamp/take4/special_cases
%
% Requires:
%       cd/isemptycell.m
%       cd/ispositiveintegervector.m
%       cd/m3ha_create_cell_info_table.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_organize_sweep_indices.m
%       cd/m3ha_select_sweeps.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/m3ha_rank_neurons.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m and beyond

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2017-05-22 Improved comments
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Now uses print_cellstr.m
% 2018-11-15 Improved documentation
% 2018-11-15 CellIds is now a numeric array
% 2018-11-19 Moved code out to m3ha_organize_sweep_indices.m
% 2018-12-05 Now allows rowsToUse to be a cellstr
% 2018-12-05 Made everything optional arguments
% 2018-12-05 Now selects cells that have sweeps for all conditions

%% Hard-coded parameters
cellNameStr = 'cellName';
toUseStr = 'toUse';

%% Default values for optional arguments
dataModeDefault = 2;
swpInfoDefault = table.empty;
cellInfoDefault = table.empty;
rowsToUseDefault = [];
casesDirDefault = '~/m3ha/data_dclamp/take4/special_cases';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'CellInfo', cellInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'RowsToUse', rowsToUseDefault, ...
    @(x) assert(ispositiveintegervector(x) || iscellstr(x) || isstring(x), ...
                ['strs5 must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'CasesDir', casesDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, varargin{:});
dataMode = iP.Results.DataMode;
swpInfo = iP.Results.SwpInfo;
cellInfo = iP.Results.CellInfo;
rowsToUse = iP.Results.RowsToUse;
casesDir = iP.Results.CasesDir;

%% Preparation
% Load sweep information if not provided
%   Note: the file names are read in as row names
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

% Generate a table of cell names from swpInfo if not provided
%   TODO? and sweep indices organized by pharm-gIncr-vHold-sweepNumber
if isempty(cellInfo)
    cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);
end

%% Select cell IDs to fit
% Update swpInfo so that there is a toUse column
%   using m3ha_select_sweeps if rowsToUse not provided
if ~isempty(rowsToUse)
    if ismember(toUseStr, swpInfo.Properties.VariableNames)
        % Print message
        fprintf('Warning: toUse column already exists in swpInfo.\n');
        fprintf('RowsToUse will be ignored!\n');
    else
        % Initialize all rows to not be fitted
        toUse = false(nRows, 1);

        % Add this to swpInfo
        swpInfo = addvars(swpInfo, toUse);

        % If rowsToUse is not empty, turn toUse on for those to fit
        swpInfo{rowsToUse, toUseStr} = true;
    end
else
    % Select the sweep indices that will be fitted
    %   Remove gIncr == 15%, 50%, 800% and problematic traces
    swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'DataMode', dataMode, ...
                                    'CasesDir', casesDir);
end

% Organize sweep indices by g incr, pharm conditions for each cell
sweepIndicesByCondition = m3ha_organize_sweep_indices('SwpInfo', swpInfo);

% Print message
fprintf('Selecting the cells to fit ... \n');

% Decide whether each cell is to be selected (has sweeps for all conditions)
isSelected = cellfun(@(x) ~any(isemptycell(x(:))), sweepIndicesByCondition);

% Store in cellInfo
cellInfo = addvars(cellInfo, sweepIndicesByCondition, isSelected);

% Find the cell IDs to be selected
cellIdsSelected = find(isSelected);

%% Print results to standard output
% Count the number of selected cells
nSelected = length(cellIdsSelected);

% Get the corresponding cell names
cellNamesSelected = cellInfo{cellIdsSelected, cellNameStr};

% Print the number of cells to fit
fprintf('There are %d cells that will be fitted: \n', nSelected);

% Print cell names to standard output
print_cellstr(cellNamesSelected, 'Delimiter', '\n', ...
            'OmitQuotes', true, 'OmitBraces', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for k = 1:nCellsToUse
    fprintf('%s\n', cellNamesSelected{k});
end
fprintf('\n');

% If no indices found, don't use this cell
if isempty(swpIndThisCell{iGIncr, iPCond})
    toUse = false;
    break;
end

% Count the total number of cells to fit
nCellsToUse = ctSelected;

% Initialize a variable to count selected cells
ctSelected = 0;

% Loop through all possible cells
for iCell = 1:nCells   
    % First assume the cell will be used
    toUse = true;
    
    % Look for all traces that can be fitted for this cell
    swpIndThisCell = cell(nGIncr, nPCond);
    for iPCond = 1:nPCond       % for each pharmacological condition
        for iGIncr = 1:nGIncr   % for each g incr in 100%, 200%, 400%
            % Get the indices for all traces of this g incr x pharm condition
            indGincrPcond = swpIdxSCPGV(:, iCell, iPCond, iGIncr + 2, :);

            % Find indices of traces for this pharm-g incr pair
            swpIndThisCell{iGIncr, iPCond} = ...
                intersect(indGincrPcond, rowsToUse, 'sorted');
        end
    end

    % If any pharm-g incr pair has no trace found, don't use this cell
    if any(any(isemptycell(swpIndThisCell)))
        toUse = false;
    end

    % If the cell is still to be used
    if toUse
        % Update cell count
        ctSelected = ctSelected + 1;

        % Store cell ID #
        cellIdsSelected(ctSelected) = iCell;             

        % Store cell name
        %   Note: This takes the first 7 characters of the first sweep's name
        %           e.g., A100810
        cellNamesSelected{ctSelected} = fileNames{swpIndThisCell{1, 1}(1)}(1:7);

        % Store sweep indices for all pharm-g incr pairs of this cell
        sweepIndicesByCondition{ctSelected} = swpIndThisCell; 
    end
end

nSelected = ctSelected;

@(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));

if ~isempty(rowsToUse)
    % Extract just the cell IDs for the sweeps to fit
    cellIdAllRows = swpInfoToUse{rowsToUse, cellidrowStr};

    % Get the sorted unique cell IDs
    cellIdsSelected = sort(unique(cellIdAllRows));
else
    % Select all cells by default
    cellIdsSelected = transpose(1:nCells);
end

% Count the number of gIncr conditions to fit
nGrowToUse = length(growToUse);

% Count the number of pharm conditions to fit
nProwToUse = length(prowToUse);

% Possible conductance amplitude scaling percentages for fitting
growToUse = [100; 200; 400]; 

% Possible pharm conditions for fitting
%   1 - Control
%   2 - GAT1 Block
%   3 - GAT3 Block
%   4 - Dual Block
prowToUse = [1; 2; 3; 4];

% Restrict table to rows to fit
if ~isempty(rowsToUse)
    swpInfoToUse = swpInfo{rowsToUse, :};
else
    swpInfoToUse = swpInfo;
end

% Count the number of cells
nCells = height(cellInfo);


[cellIdsSelected, cellNamesSelected] = m3ha_select_cells (varargin)

% Compute the number of rows
nRows = height(swpInfo);

% Generate a logical vector
toUse = false(nRows, 1);
toUse(rowsToUse) = true;

% Add the logical vector as a column to the table
swpInfo = addvars(swpInfo, toUse);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
