function [cellIdsSelected, cellInfo, swpInfo] = m3ha_select_cells (varargin)
%% Selects cells with sweeps to use for all pharm-gIncr pairs
% Usage: [cellIdsSelected, cellInfo, swpInfo] = m3ha_select_cells (varargin)
% Explanation: 
%   Finds all cell IDs and cell names (from 'CellInfo') for the rows 
%       restricted by 'ToUse' in 'SwpInfo'
%
% Examples: 
%       TODO
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
%                                       to be used
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
if isempty(cellInfo)
    cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);
end

%% Select cell IDs to fit
% Update swpInfo so that there is a toUse column
swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'RowsToUse', rowsToUse, ...
                                'DataMode', dataMode, 'CasesDir', casesDir);

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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
