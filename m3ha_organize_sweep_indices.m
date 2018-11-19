function swpIndGincrPcondAllCells = m3ha_organize_sweep_indices (swpIdxBySCPGV, swpIndToFit, cellIdsToFit)
%% Organize sweep indices by g incr, pharm conditions for each cell
% Usage: swpIndGincrPcondAllCells = m3ha_organize_sweep_indices (swpIdxBySCPGV, swpIndToFit, cellIdsToFit)
% Explanation: 
%       TODO
%
% Outputs:
%       TODO
%
% Arguments:
%       TODO
%
% Requires:
%       cd/print_cellstr.m
%
% Used by:
%       cd/singleneuronfitting42.m and later versions

% File History:
% 2018-11-19 Moved from m3ha_select_cells.m
% TODO: Input Parser

%% Hard-coded parameters
% The following must be consistent with dclampDataExtractor.m
gIncrAll = [100; 200; 400]; % possible conductance amplitude scaling percentages
pCondAll = [1; 2; 3; 4];    % possible pharm conditions 
                            %   1 - Control
                            %   2 - GAT1 Block
                            %   3 - GAT3 Block
                            %   4 - Dual Block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of possible cells
nCellsToFit = length(cellIdsToFit);

% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gIncrAll);

% Count the number of possible pharm conditions 
nPCond = length(pCondAll);

% Initialize a cell array to store sweep indices for each cell
swpIndGincrPcondAllCells = cell(nCellsToFit, 1);

%% Do the job
% Print message
fprintf('Organizing sweep indices ... \n');

% Loop through all possible cells
for iCellToFit = 1:nCellsToFit
    % Get the current cell ID
    cellId = cellIdsToFit(iCellToFit);

    % Initialize a cell array to store sweep indices for each condition
    swpIndThisCell = cell(nGIncr, nPCond);

    % Look for all traces that can be fitted for this cell
    for iPCond = 1:nPCond       % for each pharmacological condition
        for iGIncr = 1:nGIncr   % for each g incr in 100%, 200%, 400%
            % Get the indices for all traces of this g incr x pharm condition
            indGincrPcond = swpIdxBySCPGV(:, cellId, iPCond, iGIncr + 2, :);

            % Find indices of traces for this pharm-g incr pair
            swpIndThisCell{iGIncr, iPCond} = ...
                intersect(indGincrPcond, swpIndToFit, 'sorted');
        end
    end

    % Store sweep indices for all pharm-g incr pairs of this cell
    swpIndGincrPcondAllCells{iCellToFit} = swpIndThisCell; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%