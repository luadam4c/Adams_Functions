function [swpIndByCondAllCells, conditions] = m3ha_organize_sweep_indices (varargin)
%% Organize sweep indices by g incr, pharm conditions for each cell
% Usage: [swpIndByCondAllCells, conditions] = m3ha_organize_sweep_indices (varargin)
% Explanation: 
%       TODO
%
% Outputs:
%       TODO
%
% Arguments:
%       varargin    - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                       toUse       - whether the sweep is to be used (optional)
%                   default == m3ha_load_sweep_info
%                   - 'ToUse' - whether to fit each sweep
%                   must be a binary vector
%                   default == true(nRows, 1)
%                   - 'RowsToUse' - row indices or row names in swpInfo 
%                                       of sweeps to use
%                   must be a positive integer vector, a string array 
%                       or a cell array of character vectors
%                   default == []
%
% Requires:
%       cd/m3ha_load_sweep_info.m
%       cd/outer_product.m
%       cd/print_cellstr.m
%
% Used by:
%       cd/m3ha_select_cells.m

% File History:
% 2018-11-19 Moved from m3ha_select_cells.m
% 2018-12-05 Now returns indices grouped for all possible cells
% 2018-12-05 Now uses SwpInfo and RowsToUse
% 2018-12-06 Now uses a toUse column in SwpInfo

%% Hard-coded parameters
cellidrowStr = 'cellidrow';
prowStr = 'prow';
growStr = 'grow';
toUseStr = 'toUse';

%% Default values for optional arguments
swpInfoDefault = [];
toUseDefault = [];
rowsToUseDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'ToUse', toUseDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'binary', 'vector'}));
addParameter(iP, 'RowsToUse', rowsToUseDefault, ...
    @(x) assert(ispositiveintegervector(x) || iscellstr(x) || isstring(x), ...
                ['RowsToUse must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));

% Read from the Input Parser
parse(iP, varargin{:});
swpInfo = iP.Results.SwpInfo;
toUse = iP.Results.ToUse;
rowsToUse = iP.Results.RowsToUse;

%% Preparation
% Read in swpInfo if not provided
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

% Get the number of rows
nRows = height(swpInfo);

% Decide whether to fit each sweep if not provided
if isempty(toUse)
    % Decide based on whether swpInfo contains a toUse column
    %   then whether rowsToUse is provided
    if ismember(toUseStr, swpInfo.Properties.VariableNames)
        toUse = swpInfo{:, toUseStr};
    elseif ~isempty(rowsToUse)
        % Initialize all rows to not be fitted
        toUse = false(nRows, 1);

        % Add this to swpInfo
        swpInfo = addvars(swpInfo, toUse);

        % If rowsToUse is not empty, turn toUse on for those to fit
        swpInfo{rowsToUse, toUseStr} = true;

        % Extract from swpInfo
        toUse = swpInfo{:, toUseStr};
    else
        % All rows will fitted
        toUse = true(nRows, 1);
    end
end

%% Do the job
% Extract all cell numbers from swpInfo
cellIdrow = swpInfo{:, cellidrowStr};

% Get the maximum cell IDs
maxCellId = max(cellIdrow);

% Extract grow and prow from swpInfo
grow = swpInfo{:, growStr};
prow = swpInfo{:, prowStr};

% Get all unique conductance amplitude scaling percentages to fit
gCondToUse = unique(swpInfo{toUse, growStr});

% Get all unique pharmacological conditions to fit
pCondToUse = unique(swpInfo{toUse, prowStr});

% Count the number of distinct conductance amplitude scaling percentages
nGCondToUse = length(gCondToUse);

% Count the number of distinct pharmacological conditions 
nPCondToUse = length(pCondToUse);

% Initialize a cell array to store sweep indices for each cell
swpIndByCondAllCells = cell(maxCellId, 1);

%% Do the job
% Print message
fprintf('Organizing sweep indices for all cells ... \n');

% Generate strings for the conditions
gCondToUseStrs = insertBefore(cellstr(num2str(gCondToUse)), 1, ...
                                'Conductance amplitude % = ');
pCondToUseStrs = insertBefore(cellstr(num2str(pCondToUse)), 1, ...
                                'Pharm condition = ');

% Store the conditions used
conditions = outer_product(gCondToUseStrs, pCondToUseStrs);

% Loop through all possible cells
for cellIdThis = 1:maxCellId
    % Initialize a cell array to store sweep indices for each condition
    swpIndThisCell = cell(nGCondToUse, nPCondToUse);

    % Loop through each pharm-gIncr pair
    for iG = 1:nGCondToUse
        for iP = 1:nPCondToUse
            % Get the current pharm and gIncr conditions
            gThis = gCondToUse(iG);
            pThis = pCondToUse(iP);

            % Find the sweep indices that match 
            %   this cell, pharm and gIncr condition
            swpIndThisCell{iG, iP} = ...
                find(toUse & cellIdrow == cellIdThis & ...
                        grow == gThis & prow == pThis);
        end
    end

    % Store sweep indices for all pharm-g incr pairs of this cell
    swpIndByCondAllCells{cellIdThis} = swpIndThisCell; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
OLD CODE:

swpIndByCondAllCells = ...
    m3ha_organize_sweep_indices (swpIdxBySCPGV, swpIndToUse, cellIdsToUse)

% Count the number of possible cells
nCellsToUse = length(cellIdsToUse);

% Initialize a cell array to store sweep indices for each cell
swpIndByCondAllCells = cell(nCellsToUse, 1);

for iCellToUse = 1:nCellsToUse
    % Store sweep indices for all pharm-g incr pairs of this cell
    swpIndByCondAllCells{iCellToUse} = swpIndThisCell; 
end

% Count the number of cells
nCells = size(swpIdxBySCPGV, 2);

% Count the number of possible conductance amplitude scaling percentages
nGcond = length(gCondAll);

% Count the number of possible pharm conditions 
nPCond = length(pCondAll);

% Look for all traces that can be fitted for this cell
for iPCond = 1:nPCond       % for each pharmacological condition
    for iGIncr = 1:nGcond   % for each g incr in 100%, 200%, 400%
        % Get the indices for all traces of this g incr x pharm condition
        indGincrPcond = swpIdxBySCPGV(:, iCell, iPCond, iGIncr + 2, :);

        % Find indices of traces for this pharm-g incr pair
        swpIndThisCell{iGIncr, iPCond} = ...
            intersect(indGincrPcond, swpIndToUse, 'sorted');
    end
end

if iscellstr(rowsToUse) || isstring(rowsToUse)
    swpInfo.Properties.RowNames
else
    toUse(rowsToUse) = true;
end

%       cd/singleneuronfitting42.m and later versions

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%