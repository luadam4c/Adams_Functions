function cellInfo = m3ha_generate_cell_info (swpInfo, varargin)
%% TODO: A summary of what the function does (must be a single unbreaked line)
% Usage: cellInfo = m3ha_generate_cell_info (swpInfo, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       cellInfo    - a table of cell info
%                   specified as a 2D table with row indices being cell IDs 
%                       and with fields:
%                       cellName  - cell names
% Arguments:
%       swpInfo     - a table of sweep info, with each row named by 
%                       the matfile name containing the raw data
%                   must be a 2D table with row names being file names
%                       and with the fields:
%                       cellidrow - cell ID
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/m3ha_select_cells.m

% File History:
% 2018-11-19 Created by Adam Lu
% 

%% Hard-coded parameters
cellNameStr = 'cellName';

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'swpInfo', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Read from the Input Parser
parse(iP, swpInfo, varargin{:});

%% Do the job
% Print message
fprintf('Generating cell info ... \n');

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