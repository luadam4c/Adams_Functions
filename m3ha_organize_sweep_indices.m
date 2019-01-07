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
%                               the matfile name containing the raw data
%                   must a 2D table with row names being file names
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                   default == loaded from 
%                       ~/m3ha/data_dclamp/take4/dclampdatalog_take4.csv
%                   - 'ToFit' - whether to fit each sweep
%                   must be a binary vector
%                   default == true(nRows, 1)
%                   - 'RowsToFit' - row indices or row names in swpInfo 
%                                       of sweeps to fit
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
% 2018-12-05 Now uses SwpInfo and RowsToFit
% 2018-12-06 Now uses a toFit column in SwpInfo

%% Hard-coded parameters
cellidrowStr = 'cellidrow';
prowStr = 'prow';
growStr = 'grow';
toFitStr = 'toFit';

%% Default values for optional arguments
swpInfoDefault = [];
toFitDefault = [];
rowsToFitDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'ToFit', toFitDefault, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'binary', 'vector'}));
addParameter(iP, 'RowsToFit', rowsToFitDefault, ...
    @(x) assert(ispositiveintegervector(x) || iscellstr(x) || isstring(x), ...
                ['RowsToFit must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));

% Read from the Input Parser
parse(iP, varargin{:});
swpInfo = iP.Results.SwpInfo;
toFit = iP.Results.ToFit;
rowsToFit = iP.Results.RowsToFit;

%% Preparation
% Read in swpInfo if not provided
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

% Get the number of rows
nRows = height(swpInfo);

% Decide whether to fit each sweep if not provided
if isempty(toFit)
    % Decide based on whether swpInfo contains a toFit column
    %   then whether rowsToFit is provided
    if ismember(toFitStr, swpInfo.Properties.VariableNames)
        toFit = swpInfo{:, toFitStr};
    elseif ~isempty(rowsToFit)
        % Initialize all rows to not be fitted
        toFit = false(nRows, 1);

        % Add this to swpInfo
        swpInfo = addvars(swpInfo, toFit);

        % If rowsToFit is not empty, turn toFit on for those to fit
        swpInfo{rowsToFit, toFitStr} = true;

        % Extract from swpInfo
        toFit = swpInfo{:, toFitStr};
    else
        % All rows will fitted
        toFit = true(nRows, 1);
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
gCondToFit = unique(swpInfo{toFit, growStr});

% Get all unique pharmacological conditions to fit
pCondToFit = unique(swpInfo{toFit, prowStr});

% Count the number of distinct conductance amplitude scaling percentages
nGCondToFit = length(gCondToFit);

% Count the number of distinct pharmacological conditions 
nPCondToFit = length(pCondToFit);

% Initialize a cell array to store sweep indices for each cell
swpIndByCondAllCells = cell(maxCellId, 1);

%% Do the job
% Print message
fprintf('Organizing sweep indices for all cells ... \n');

% Generate strings for the conditions
gCondToFitStrs = insertBefore(cellstr(num2str(gCondToFit)), 1, ...
                                'Conductance amplitude % = ');
pCondToFitStrs = insertBefore(cellstr(num2str(pCondToFit)), 1, ...
                                'Pharm condition = ');

% Store the conditions used
conditions = outer_product(gCondToFitStrs, pCondToFitStrs);

% Loop through all possible cells
for cellIdThis = 1:maxCellId
    % Initialize a cell array to store sweep indices for each condition
    swpIndThisCell = cell(nGCondToFit, nPCondToFit);

    % Loop through each pharm-gIncr pair
    for iG = 1:nGCondToFit
        for iP = 1:nPCondToFit
            % Get the current pharm and gIncr conditions
            gThis = gCondToFit(iG);
            pThis = pCondToFit(iP);

            % Find the sweep indices that match 
            %   this cell, pharm and gIncr condition
            swpIndThisCell{iG, iP} = ...
                find(toFit & cellIdrow == cellIdThis & ...
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
    m3ha_organize_sweep_indices (swpIdxBySCPGV, swpIndToFit, cellIdsToFit)

% Count the number of possible cells
nCellsToFit = length(cellIdsToFit);

% Initialize a cell array to store sweep indices for each cell
swpIndByCondAllCells = cell(nCellsToFit, 1);

for iCellToFit = 1:nCellsToFit
    % Store sweep indices for all pharm-g incr pairs of this cell
    swpIndByCondAllCells{iCellToFit} = swpIndThisCell; 
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
            intersect(indGincrPcond, swpIndToFit, 'sorted');
    end
end

if iscellstr(rowsToFit) || isstring(rowsToFit)
    swpInfo.Properties.RowNames
else
    toFit(rowsToFit) = true;
end

%       cd/singleneuronfitting42.m and later versions

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%