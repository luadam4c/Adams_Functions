function cellInfo = m3ha_generate_cell_info (varargin)
%% Generates a table of cell information from the sweep information table
% Usage: cellInfo = m3ha_generate_cell_info (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       cellInfo    - a table of cell info
%                   specified as a 2D table with row indices being cell IDs 
%                       and with fields:
%                       cellName   - cell names
%                       TODO? swpIndices - sweep indices in original swpInfo
%                                       organized by pharm-gincr-vhold-sweep
%                                       conditions
% Arguments:
%       varargin    - 'SwpInfo': a table of sweep info, with each row named by 
%                       the matfile name containing the raw data
%                   must be a 2D table with row names being file names
%                       and with the fields:
%                       cellidrow   - cell ID
%                       TODO?
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                       vrow        - holding voltage level (mV)
%                       swpnrow     - sweep number
%                   default == loaded from 
%                       ~/m3ha/data_dclamp/take4/dclampdatalog_take4.csv
%                   - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/m3ha_load_sweep_info.m
%
% Used by:
%       cd/m3ha_select_cells.m

% File History:
% 2018-11-19 Created by Adam Lu
% 

%% Hard-coded parameters
cellNameStr = 'cellName';

%% Default values for optional arguments
swpInfoDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
swpInfo = iP.Results.SwpInfo;

%% Preparation
% Read in swpInfo if not provided
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

%% Generate cell name info
% Print message
fprintf('Generating cell name info ... \n');

% Extract all cell IDs
cellIdAllRows = swpInfo{:, 'cellidrow'};

% Extract all file names
fileNames = swpInfo.Properties.RowNames;

% Get the maximum cell ID
maxCellId = max(cellIdAllRows);

% For each possible cell ID, look for the first row with the cell ID in question
indFirstRow = arrayfun(@(x) find(cellIdAllRows == x, 1), ...
                        transpose(1:maxCellId), 'UniformOutput', false);

% Find the appropriate cell names for each row
cellNames = cellfun(@(x) find_cell_name(x, fileNames), indFirstRow, ...
                    'UniformOutput', false);

% Put in a table
cellInfo = table(cellNames, 'VariableName', {cellNameStr});

%% Organize sweep indices by pharm, g incr, vHold, sweep # for each cell
% TODO?
fprintf('Organizing sweep indices ... \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cellName = find_cell_name (idxFirstRow, fileNames)
%% Find the appropriate cell name

% Decide on the appropriate cell name
if isempty(idxFirstRow)
    % If no row is found, give the cell name 'none'. 
    cellName = 'none';
else
    %   Otherwise, use the first 7 characters of the file name
    cellName = fileNames{idxFirstRow}(1:7);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the sorted unique cell IDs
uniqueCellIds = sort(unique(cellIdAllRows));

% Count the number of selected cells
nCells = length(uniqueCellIds);

% Get the corresponding cell names
cellNamesSelected = cellInfo{uniqueCellIds, 'cellName'};

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%