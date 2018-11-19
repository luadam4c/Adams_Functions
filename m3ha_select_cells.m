function [cellIdsSelected, cellNamesSelected] = ...
                m3ha_select_cells (rowsToFit, swpInfo, varargin)
%% Select the cells that will be used for fitting from the sweeps to fit
% Usage: [cellIdsSelected, cellNamesSelected] = ...
%               m3ha_select_cells (rowsToFit, swpInfo, varargin)
% Explanation: 
%   Finds all cell ID and cell names for the rows rowsToFit in swpInfo
%
% Outputs:
%       TODO
%
% Arguments:
%       rowsToFit - row indices in swpInfo of sweeps to fit
%                   must be a positive integer vector
%       swpInfo     - a table of sweep info, with each row named by 
%                       the matfile name containing the raw data
%                   must be a 2D table with row names being file names
%                       and with the fields:
%                       cellidrow - cell ID
%       varargin    - 'CellInfo': data mode
%                   must a 2D table with row indices being cell IDs 
%                       and with fields:
%                       cellName - cell names
%                   default == []
%
% Requires:
%       cd/m3ha_generate_cell_info.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/singleneuronfitting42.m and later versions

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2017-05-22 Improved comments
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Now uses print_cellstr.m
% 2018-11-15 Improved documentation
% 2018-11-15 CellIds is now a numeric array
% 2018-11-19 Moved code out to m3ha_organize_sweep_indices.m

%% Default values for optional arguments
cellInfoDefault = [];

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
addRequired(iP, 'rowsToFit', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));
addRequired(iP, 'swpInfo', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'CellInfo', cellInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Read from the Input Parser
parse(iP, rowsToFit, swpInfo, varargin{:});
cellInfo = iP.Results.CellInfo;

%% Preparation
% Generate cellInfo from swpInfo if not provided
if isempty(cellInfo)
    cellInfo = m3ha_generate_cell_info(swpInfo);
end

%% Do the job
% Print message
fprintf('Selecting the cells to fit ... \n');

% Extract just the cell IDs for the sweeps to fit
cellIdAllRows = swpInfo{rowsToFit, 'cellidrow'};

% Get the sorted unique cell IDs
cellIdsSelected = sort(unique(cellIdAllRows));

% Count the number of selected cells
nSelected = length(cellIdsSelected);

% Get the corresponding cell names
cellNamesSelected = cellInfo{cellIdsSelected, 'cellName'};

%% Print results to standard output
% Print the number of cells to fit
fprintf('There are %d cells that will be fitted: \n', nSelected);

% Print cell names to standard output
print_cellstr(cellNamesSelected, 'Delimiter', '\n', ...
            'OmitQuotes', true, 'OmitBraces', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for k = 1:nCellsToFit
    fprintf('%s\n', cellNamesSelected{k});
end
fprintf('\n');

% If no indices found, don't use this cell
if isempty(swpIndThisCell{iGIncr, iPCond})
    toUse = false;
    break;
end

% Count the total number of cells to fit
nCellsToFit = ctSelected;

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
                intersect(indGincrPcond, rowsToFit, 'sorted');
        end
    end

    % If any pharm-g incr pair has no trace found, don't use this cell
    if any(any(cellfun(@isempty, swpIndThisCell)))
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
        swpIndGincrPcondAllCells{ctSelected} = swpIndThisCell; 
    end
end

nSelected = ctSelected;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%