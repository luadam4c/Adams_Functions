function [cellIds, cellNames, swpIndGincrPcondAllCells] = ...
            m3ha_select_cells (swpIdxSCPGV, swpIndToFit, fileNames)
%% Select the cells that will be used for fitting
% Usage: 
% Explanation: 
%   Finds cells with traces present in swpIndToFit
%       for all g_incr - pharm condition pairs
%
% Example(s):
%       TODO
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
%
% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2017-05-22 Improved comments
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Now uses print_cellstr.m
% 2018-11-15 Improved documentation

%% Hard-coded parameters
% The following must be consistent with dclampDataExtractor.m
gIncrAll = [100; 200; 400]; % possible conductance amplitude scaling percentages
pCondAll = [1; 2; 3; 4];    % possible pharm conditions 
                            %   1 - Control
                            %   2 - GAT1 Block
                            %   3 - GAT3 Block
                            %   4 - Dual Block
cellIdAll = 1:1:49;         % possible cell ID #s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of possible cells
nCells = length(cellIdAll);

% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gIncrAll);

% Count the number of possible pharm conditions 
nPCond = length(pCondAll);

%% Do the job
% Print message
fprintf('Selecting what cells to fit... \n');

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
                intersect(indGincrPcond, swpIndToFit, 'sorted');
        end
    end

    % If any pharm-g incr pair has no trace found, don't use this cell
    if any(cellfun(@isempty, swpIndThisCell))
        toUse = false;
    end

    % If the cell is still to be used
    if toUse
        % Update cell count
        ctSelected = ctSelected + 1;

        % Store cell ID #
        cellIds{ctSelected} = iCell;             

        % Store cell name
        %   Note: This takes the first 7 characters of the first sweep's name
        %           e.g., A100810
        cellNames{ctSelected} = fileNames{swpIndThisCell{1, 1}(1)}(1:7);

        % Store sweep indices for all pharm-g incr pairs of this cell
        swpIndGincrPcondAllCells{ctSelected} = swpIndThisCell; 
    end
end

%% Print results to standard output
% Print the number of cells to fit
fprintf('There are %d cells that will be fitted: \n', ctSelected);

% Print cell names to standard output
print_cellstr(cellNames, 'Delimiter', '\n', ...
            'OmitQuotes', true, 'OmitBraces', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

for k = 1:nCellsToFit
    fprintf('%s\n', cellNames{k});
end
fprintf('\n');

% If no indices found, don't use this cell
if isempty(swpIndThisCell{iGIncr, iPCond})
    toUse = false;
    break;
end

% Count the total number of cells to fit
nCellsToFit = ctSelected;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%