function [fileNames, rowConditions, figurePositions] = ...
                m3ha_select_raw_traces (rowMode, columnMode, ...
                    attemptNumber, cellIdToFit, swpInfo, cellInfo)
%% Select raw traces to import
% Usage: [fileNames, rowConditions, figurePositions] = ...
%               m3ha_select_raw_traces (rowMode, columnMode, ...
%                   attemptNumber, cellIdToFit, swpInfo, cellInfo)
%
% Requires:
%       cd/m3ha_determine_row_conditions.m
%
% Used by:
%       cd/singleneuronfitting42.m and later versions
%
% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2017-08-11 Added Attempt #5 for fitting across cells
% 2017-08-21 Modified Attempt #4 for fitting across trials so that
%               the "most representative trial" is selected
% 2017-08-21 Added maxNoise as an argument
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation and code clarity
% 2018-12-06 Now uses cellIdThis, fileNamesToFit, swpInfo, cellInfo

%% Hard-coded parameters
prowStr = 'prow';
growStr = 'grow';
burstTimeStr = 'bursttime';
ltsPeakTimeStr = 'ltspeaktime';
maxNoiseStr = 'maxnoise';
toFitStr = 'toFit';
cellNameStr = 'cellName';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Extract from swpInfo
fnrow = swpInfo.Properties.RowNames;
grow = swpInfo{:, growStr};
prow = swpInfo{:, prowStr};
burstTime = swpInfo{:, burstTimeStr};
ltsPeakTime = swpInfo{:, ltsPeakTimeStr};
maxNoise = swpInfo{:, maxNoiseStr};
toFit = swpInfo{:, toFitStr};

% Get all unique conductance amplitude scaling percentages to fit
gCondToFit = unique(swpInfo{toFit, growStr});

% Get all unique pharmacological conditions to fit
pCondToFit = unique(swpInfo{toFit, prowStr});

% Count the number of distinct conductance amplitude scaling percentages
nGCondToFit = length(gCondToFit);

% Count the number of distinct pharmacological conditions 
nPCondToFit = length(pCondToFit);

% Get the current cell name(s) from cellInfo
cellName = cellInfo{cellIdToFit, cellNameStr};
if numel(cellName) == 1
    cellName = cellName{1};
end

%% Determine row conditions
% Determine what each row would be
rowConditions = ...
    m3ha_determine_row_conditions(rowMode, columnMode, attemptNumber, ...
                                    gCondToFit, pCondToFit);

% Count the number of rows
nRows = size(rowConditions, 1);

%% Determine what each sweep satisfies
% Determine whether each sweep has bursts
hasBurst = burstTime > 0;

% Determine whether each sweep has LTSs
hasLts = ltsPeakTime > 0;

% Determine whether each sweep is from the cell(s) to fit
%   Note: cellName may have more than one elements
fromCell = contains(fnrow, cellName);

% Determine whether each sweep is of each pharm condition
%   Note: isPCond is a cell array of logical vectors
isPCond = arrayfun(@(x) prow == x, pCondToFit, 'UniformOutput', false);

% Determine whether each sweep is of each gIncr condition
%   Note: isGCond is a cell array of logical vectors
isGCond = arrayfun(@(x) grow == x, gCondToFit, 'UniformOutput', false);

% Determine whether each sweep satisfies the conditions for each row
%   Note: isRowCond is a cell array of logical vectors
if size(rowConditions, 2) == 1
    % Each row is the iPCond
    iPEachRow = rowConditions;

    % Use isPCond
    %   Note: isRowCond might be a reordering of isPCond
    isRowCond = arrayfun(@(x) isPCond{x}, iPEachRow, 'UniformOutput', false);
elseif size(rowConditions, 2) == 2
    % Extract the iPCond for each row
    iPEachRow = rowConditions(:, 1);

    % Extract the iGincr for each row
    iGEachRow = rowConditions(:, 2);

    % Use isPCond and isGCond
    isRowCond = arrayfun(@(x, y) isPCond{x} & isGCond{y}, ...
                        iPEachRow, iGEachRow, 'UniformOutput', false);
else
    error('This rowConditions is not recognized!');
end

%% Select sweeps
% Print message
fprintf('Selecting the traces for this cell to fit ... \n');

% Initialize nColumns
nColumns = 0;

% Initialize various schemes
% TODO
swpIndTemp1 = cell(1, nPCondToFit);
swpIndTemp2 = cell(1, nPCondToFit);
swpIndTemp3 = cell(1, nPCondToFit);
swpIndTemp4 = cell(1, nPCondToFit);
swpIndTemp5 = cell(nGCondToFit, nPCondToFit);

% Decide on the sweeps to use for each condition
switch columnMode
case 1
    % Print message
    fprintf('Fitting across trials of cell %s ... \n', cellName);

    % Decide on the indices based on attempt number
    switch attemptNumber
    case {1, 2}
        % Print message
        if attemptNumber == 1
            fprintf('Attempt #1: Using 4 traces of %s @ 200% gIncr ... \n', ...
                    cellName);
        elseif attemptNumber == 2
            fprintf('Attempt #2: Using all traces of %s @ 200% gIncr ... \n', ...
                    cellName);
        end

        % Find the sweep indices for each pharmacological condition 
        %   @ 200% g_incr from this cell to be fitted
        swpIndRow = cellfun(@(x) find(x & isGCond{2} & fromCell & toFit), ...
                                isPCond, 'UniformOutput', false);

        % Check if there is enough data or not
        if any(cellfun(@isempty, swpIndRow))
            error('There is not enough data for %s @ 200% gIncr!', cellName);
        end

        % Decide on the number of columns per row
        if attemptNumber == 1
            % Use only the first column
            nColumns = 1;
        elseif attemptNumber == 2
            % Make sure each row has the same number of traces
            nColumns = min(cellfun(@length, swpIndRow));
        end
    case 3
        % Print message
        fprintf('Attempt #3: Using all traces of %s ... \n', cellName);

        % Find the sweep indices for each row condition
        %   from this cell to be fitted
        swpIndRow = cellfun(@(x) find(x & fromCell & toFit), ...
                            isRowCond, 'UniformOutput', false);

        % Check if there is enough data or not
        if any(cellfun(@isempty, swpIndRow))
            error('There is not enough data for %s!', cellName);
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, swpIndRow));
    case 4
        % Print message
        fprintf('Attempt #4: Using 12 traces (one "best trial" for each ');
        fprintf('condition) of %s ... \n', cellName);

        % Take only the most representative trace from each condition
        if size(rowConditions, 2) == 1
            for iRow = 1:nRows
                % Get the corresponding iPCond
                iPThis = iPEachRow(iRow);

                % Initialize
                swpIndThisRow = [];

                % Get three traces
                for iG = 1:nGCondToFit
                    % Determine whether each sweep is this condition
                    isThisCond = isGCond{iG} & isPCond{iPThis} & ...
                                    fromCell & toFit;

                    % Select the "most representative trace"
                    swpIdxSelected = ...
                        select_trace(isThisCond, hasLts, hasBurst, maxNoise);

                    if isempty(swpIdxSelected)
                        error(['Most representative trace not found', ...
                                ' for iG == %d, iP == %d!'], iG, iPThis);
                    else
                        swpIndThisRow = [swpIndThisRow; swpIdxSelected];
                    end
                end

                % Store
                swpIndRow{iRow} = swpIndThisRow;
            end
        elseif size(rowConditions, 2) == 2
            swpIndRow = arrayfun(@(x) select_trace(x & fromCell & toFit, ...
                                            hasLts, hasBurst, maxNoise), ...
                                isRowCond, 'UniformOutput', false);
        end

        % Check if there is enough data or not
        if any(cellfun(@isempty, swpIndRow))
            error('There is not enough data for %s!', cellName);
        end

        % Find the number of traces per row desired
        if size(rowConditions, 2) == 1
            nColumns = 3;
        elseif size(rowConditions, 2) == 2
            nColumns = 1;
        end
    end
case 2
% TODO TODO TODO: FIX
    fprintf('Fitting across cells ... \n');
    switch attemptNumber
    case 1
        % Attempt #1: Find cells with LTSs present for all pharm conditions @ 200% g_incr
        fprintf('Attempt #1: Find cells with LTSs present ');
        fprintf('for all pharm conditions @ 200p g_incr ... \n');
        ct = 0;                 % counts cells with an lts for all pharm conditions
        for iC = 1:length(cellIdAll)
            toUse = true;
            for iP = 1:nPCondToFit       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasLts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCondToFit   % for each pharmacological condition @ 200% g_incr
                    swpIndRow{iP}(ct) = swpIndTemp1{iP}(1);
                end
            end
        end
        nColumns = ct;                     % In our data set, nColumns == 13
    case 2
        % Attempt #2: Find cells with bursts present for all pharm conditions @ 200% g_incr
        fprintf('Attempt #2: Find cells with bursts present ');
        fprintf('for all pharm conditions @ 200p g_incr ... \n');
        ct = 0;                % counts cells with an lts for all pharm conditions
        for iC = 1:length(cellIdAll)
            toUse = true;
            for iP = 1:nPCondToFit       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasBursts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCondToFit   % for each pharmacological condition @ 200% g_incr
                    swpIndRow{iP}(ct) = swpIndTemp1{iP}(1);
                end
            end
        end
        nColumns = ct;                     % In our data set, nColumns == 8
    case 3
        % Attempt #3: Find cells with bursts present for 3 out of 4 pharm conditions @ 200% g_incr, 
        %        and choose the "best" sweep for each cell
        fprintf('Attempt #2: Find cells with bursts present ');
        fprintf('for 3 out of 4 pharm conditions @ 200p g_incr, \n'); 
        fprintf('        and choose the "best" sweep for each cell ... \n');
        ct = 0;             % counts cells with a burst for 3 out of 4 pharm conditions
        for iC = 1:length(cellIdAll)
            toUse = true;
            ngoodb = 0;
            for iP = 1:nPCondToFit       % for each pharmacological condition @ 200% g_incr
                swpIndTemp2{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), ...
                                swpIndHasBursts, 'sorted');
                swpIndTemp3{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), ...
                                swpIndHasLts, 'sorted');
                swpIndTemp4{iP} = ...
                    setdiff(swpIdxSCPGV(:, iC, iP, 4, :), 0);
                if ~isempty(swpIndTemp2{iP})
                    ngoodb = ngoodb + 1;
                end
                if isempty(swpIndTemp4{iP})
                    toUse = false;
                end
            end
            if ngoodb < 3
                toUse = false;
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCondToFit   % for each pharmacological condition @ 200% g_incr
                    if ~isempty(swpIndTemp2{iP})
                        swpIndRow{iP}(ct) = swpIndTemp2{iP}(1);
                    elseif ~isempty(swpIndTemp3{iP})
                        swpIndRow{iP}(ct) = swpIndTemp3{iP}(1);
                    else
                        swpIndRow{iP}(ct) = swpIndTemp4{iP}(1);
                    end
                end
            end
        end
        nColumns = ct;                     % In our data set, nColumns == 22
    case 4
        % Attempt #4: Find cells within indices to fit present for all pharm conditions and all g_incr
        fprintf('Attempt #4: Find cells within indices to fit present \n');
        fprintf('        for all pharm conditions and all g_incr ... \n');
        ct = 0;                % counts cells
        for iC = 1:length(cellIdAll)
            toUse = true;
            for iP = 1:nPCondToFit       % for each pharmacological condition
                for iG = 1:nGCondToFit    % for each g_incr
                    swpIndTemp5{iG, iP} = ...
                        intersect(swpIdxSCPGV(:, iC, iP, iG + 2, :), ...
                                    swpIndToFit, 'sorted');
                    if isempty(swpIndTemp5{iG, iP})
                        toUse = false;
                        break;
                    end
                end
            end
            if toUse 
                ct = ct + 1;
                iRow = 0;
                for iP = 1:nPCondToFit   % for each pharmacological condition
                    for iG = 1:nGCondToFit    % for each g incr
                        iRow = iRow + 1;
                        swpIndRow{iRow}(ct) = swpIndTemp5{iG, iP}(1);
                    end
                end
            end
        end
        nColumns = ct;                % In our data set, nColumns == 36
    case 5
        % Attempt #5: Choose the "most representative" sweep from each cell in iCellToFit
        fprintf(['Attempt #5: Choose the "most representative" sweep for each cell ', ...
                    'in cells #%d to #%d ... \n'], iCellToFit(1), iCellToFit(end));
        % Get number of cells to fit
        ncellstofit = length(iCellToFit);

        % Take a trace from each cell for each condition
        %   with priority given to a trace with bursts,
        %   then to a trace with LTSs
        if strcmpi(rowType, 'pharm')
            % Number of traces per row is 3 times number of cells to fit
            nColumns = 3*ncellstofit;

            for iRow = 1:nRows          % for each pharm condition
                % Initialize swpIndRow
                swpIndRow{iRow} = [];     

                % Build swpIndRow
                for iG = 1:nGCondToFit   % for each g incr
                    for iFit = iCellToFit             % for each cell to fit

                        indThisCond = swpIndByCondition{iFit}{iG, iPEachRow(iRow)};

                        % Select "most representative trace"
                        swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                                    swpIndHasBursts, maxNoise);
                        swpIndRow{iRow} = [swpIndRow{iRow}; swpIdxSelected];
                    end
                end

                % Verify length of swpIndRow{iRow}
                if length(swpIndRow{iRow}) ~= nColumns
                    error(['Not enough traces selected for row #%d! ', ...
                            'Something is wrong!\n\n'], iRow);
                end
            end            
        elseif strcmpi(rowType, 'pharm-gincr')
            % Number of traces per row equals the number of cells to fit
            nColumns = ncellstofit;

            % Take only the first trace from each cell for each condition
            for iRow = 1:nRows          % for each pharm condition
                % Initialize swpIndRow
                swpIndRow{iRow} = [];     

                % Build swpIndRow
                for iFit = iCellToFit             % for each cell to fit
                    % TODO
                    indThisCond = swpIndByCondition{iFit}{iGEachRow(iRow), ...
                                                            iPEachRow(iRow)};

                    % Select "most representative trace"
                    swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                                    swpIndHasBursts, maxNoise);
                    swpIndRow{iRow} = [swpIndRow{iRow}; swpIdxSelected];
                end

                % Verify length of swpIndRow{iRow}
                if length(swpIndRow{iRow}) ~= nColumns
                    error(['Not enough traces selected for row #%d! ', ...
                            'Something is wrong!\n\n'], iRow);
                end
            end
        end
    end

otherwise
    error('Column mode undefined!');
end

% Make sure traces are found
if nColumns == 0
    error('No traces found!')
end

% Count the number of traces
nSweeps = nRows * nColumns;

% Get the figure positions and sweep indices for each trace
swpIndices = zeros(nSweeps, 1);
ct = 0;
for iRow = 1:nRows
    for iCol = 1:nColumns
        % Increment count
        ct = ct + 1;

        % Get the current sweep index
        swpIndices(ct) = swpIndRow{iRow}(iCol);
    end
end

% Get the file names
fileNames = fnrow(swpIndices);

% Get the figure positions for each trace
% TODO: Do we need figurePositions here?
figurePositions = cell(nSweeps, 1);
ct = 0;
for iRow = 1:nRows
    for iCol = 1:nColumns
        % Increment count
        ct = ct + 1;

        % Get the current figure position
        figurePositions{ct} = [iRow, iCol];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function swpIdxSelected = select_trace(isThisCond, hasLts, hasBursts, maxNoise)

% Return error if there are no traces for this condition
if ~any(isThisCond)
    error('There must be traces for this cell & condition!\n');
end

% Decide whether each trace in this condition has an LTS or not
isThisCondHasNoLts = isThisCond & ~hasLts;
isThisCondHasLts = isThisCond & hasLts;

% Decide whether each trace in this condition with LTS has a burst or not
isThisCondHasLtsOnly = isThisCondHasLts & ~hasBursts;
isThisCondHasBurst = isThisCond & hasBursts;

% Count the number of traces in this condition
nIndThisCond = sum(isThisCond);

% Count the number of traces in this condition without LTS
nIndThisCondHasNoLts = sum(isThisCondHasNoLts);

% Count the number of traces that have an LTS but not a burst
nIndThisCondHasLtsOnly = sum(isThisCondHasLtsOnly);

% Count the number of traces that have a burst
nIndThisCondHasBurst = sum(isThisCondHasBurst);

% Find the corresponding indices
indThisCondHasNoLts = find(isThisCondHasNoLts);
indThisCondHasLtsOnly = find(isThisCondHasLtsOnly);
indThisCondHasBurst = find(isThisCondHasBurst);

% Select one "most representative" trace
if nIndThisCondHasNoLts > nIndThisCond / 2
    % Look for the sweep with minimum noise in all traces with no LTS
    [~, iTemp] = min(maxNoise(indThisCondHasNoLts));
    swpIdxSelected = indThisCondHasNoLts(iTemp);
elseif nIndThisCondHasLtsOnly > nIndThisCondHasBurst
    % Look for the sweep with minimum noise in all traces with LTS but no burst
    [~, iTemp] = min(maxNoise(indThisCondHasLtsOnly));
    swpIdxSelected = indThisCondHasLtsOnly(iTemp);
elseif nIndThisCondHasBurst >= nIndThisCondHasLtsOnly
    % Look for the sweep with minimum noise in all traces with burst
    [~, iTemp] = min(maxNoise(indThisCondHasBurst));
    swpIdxSelected = indThisCondHasBurst(iTemp);
else
    % Something is wrong 
    error('Error in LTS or burst data??\n');   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%        ng200p_ind(iP) = ct;
%        nColumns = min(min(ng200p_ind), 20);        
% number of cells to take per pharmacological condition
%                                
% TO FIX: doesn't work if too many cells are chosen

nColumns = min(nColumns, 5);

swpIndRow{iRow} = [swpIndRow{iRow}; ...
            swpIndGincrPcond{iG, iPEachRow(iRow)}(1)];

% Select one trace with the priority given to 
%   bursts, then to LTSs
if ~isempty(indThisCondHasBurst)
    swpIdxSelected = indThisCondHasBurst(1);
elseif ~isempty(indThisCondHasLts)
    swpIdxSelected = indThisCondHasLts(1);
else
    swpIdxSelected = indThisCond(1);
end

% Find all indices for this cell
indThisCell = swpIndByCondition{iFit};

% Find all indices for this condition for this cell
indThisCond = indThisCell{iGEachRow(iRow), iPEachRow(iRow)};

% Find all indices of traces that have an LTS
indThisCondHasLts = ...
    intersect(indThisCond, swpIndHasLts, 'sorted');

% Find all indices of traces that have a burst
indThisCondHasBurst = ...
    intersect(indThisCond, swpIndHasBursts, 'sorted');

% Select one trace with the priority given to 
%   bursts, then LTSs
if ~isempty(indThisCondHasBurst)
    swpIdxSelected = indThisCondHasBurst(1);
elseif ~isempty(indThisCondHasLts)
    swpIdxSelected = indThisCondHasLts(1);
else
    swpIdxSelected = indThisCond(1);
end

for iRow = 1:nRows
    swpIndRow{iRow} = [];
    for iG = 1:nGCondToFit
        swpIndRow{iRow} = [swpIndRow{iRow}; ...
                    swpIndGincrPcond{iG, iPEachRow(iRow)}];
    end
end

for iRow = 1:nRows
    swpIndRow{iRow} = transpose([swpIndGincrPcond{:, iRow}]);
end

for iRow = 1:nRows
    swpIndRow{iRow} = ...
        swpIndGincrPcond{iGEachRow(iRow), iPEachRow(iRow)};
end


% Get the file names
fileNames = cell(nRows, nColumns);
for iRow = 1:nRows
    for iCol = 1:nColumns
        % Get the current sweep index
        if ~isempty(swpIndG200P{1})
            swpIdx = swpIndG200P{iRow}(iCol);
        elseif ~isempty(swpIndPCond{1})
            swpIdx = swpIndPCond{iRow}(iCol);
        elseif ~isempty(swpIndRow{1})
            swpIdx = swpIndRow{iRow}(iCol);
        end

        % TODO
        fileNames{iRow, iCol} = fnrow{swpIdx};
    end
end

swpIndE091710 = find(~cellfun(@isempty, strfind(fnrow, 'E091710')));

% Get the sweep indices for each g incr-pharm condition
swpIndGincrPcond = swpIndByCondition{iCellToFit};

[fileNames, swpIndices, figurePositions] = ...
                m3ha_select_raw_traces (columnMode, rowConditions, ...
                    attemptNumber, iCellToFit, swpIdxSCPGV, swpIndToFit, ...
                    swpIndByConditionAllCells, cellNamesToFit, swpInfo)

cellIdAll = 1:1:49;         % possible cell ID #s

fnrow = swpInfo.fnrow;

cellName = cellNamesToFit{iCellToFit};

% Get the sweep indices for each g incr-pharm condition
swpIndGincrPcond = swpIndByConditionAllCells{iCellToFit};

swpIndG200P{iP} = ...
    intersect(swpIdxSCPGV(:, :, iP, 4, :), swpIndE091710, 'sorted');

for iP = 1:nPCondToFit        % for each pharmacological condition
    swpIndPCond{iP} = ...
        intersect(swpIdxSCPGV(:, :, iP, :, :), swpIndE091710, 'sorted');
    if isempty(swpIndPCond{iP})
        error('Does not have enough data for E091710!')
    end
end

% Determine whether each sweep is from the cell E091710
isE091710 = contains(fnrow, 'E091710');

if attemptNumber < 3
    fprintf('Fitting across trials of cell E091710 ... \n');
else
    fprintf('Fitting across trials of cell %s ... \n', cellName);
end

% Find all sweeps for the cell E091710
swpIndE091710 = find(isE091710);

swpIndG200P = cell(1, nPCondToFit);
swpIndPCond = cell(1, nPCondToFit);

if ~isempty(swpIndG200P{1})
    swpIndices(ct) = swpIndG200P{iRow}(iCol);
elseif ~isempty(swpIndPCond{1})
    swpIndices(ct) = swpIndPCond{iRow}(iCol);
elseif ~isempty(swpIndRow{1})
    swpIndices(ct) = swpIndRow{iRow}(iCol);
end

swpIndRow = cell(1, nRows);

% The following must be consistent with dclampDataExtractor.m
gCondToFit = [100; 200; 400]; % possible conductance amplitude scaling percentages
pCondToFit = [1; 2; 3; 4];    % possible pharm conditions 
                            %   1 - Control
                            %   2 - GAT1 Block
                            %   3 - GAT3 Block
                            %   4 - Dual Block

% Find all sweep indices with bursts
swpIndHasBursts = find(hasBurst);

% Find all sweep indices with LTSs
swpIndHasLts = find(hasLts);

fprintf('Attempt #2: Using all traces of %s ... \n', cellName);

% Attempt #3: Use all traces of each cell to fit @ 100%, 200%, 400% g_incr
fprintf(['Attempt #3: Using all traces of ', ...
            '%s @ 100p, 200p, 400p g_incr ... \n'], cellName);

% Use all data already filtered
if strcmpi(rowType, 'pharm')
    swpIndRow = arrayfun(@(x) reshape([swpIndGincrPcond{:, x}], [], 1), ...
                        transpose(1:nRows), 'UniformOutput', false);

elseif strcmpi(rowType, 'pharm-gincr') 
    swpIndRow = arrayfun(@(x, y) swpIndGincrPcond{x, y}, ...
                        iGEachRow, iPEachRow, 'UniformOutput', false);
end

indThisCond = swpIndGincrPcond{iGEachRow(iRow), iPEachRow(iRow)};

if strcmpi(rowType, 'pharm')
    % Each row is a pharmacological condition
    swpIndRow = arrayfun(@(x) find(isPCond{x} & fromCell & toFit), ...
                        iPEachRow, 'UniformOutput', false);
elseif strcmpi(rowType, 'pharm-gincr') 
    swpIndRow = arrayfun(@(x, y) find(isGCond{x} & isPCond{y} & ...
                                        fromCell & toFit), ...
                        iGEachRow, iPEachRow, 'UniformOutput', false);
end

% Display warning if row type is incorrect
if ~strcmpi(rowType, 'pharm')
    error('%s expected to be ''pharm'' for attempt number %g', ...
            rowType, attemptNumber);
end

% Each row is one pharm condition
rowType = 'pharm';
% Each row is a pharm condition paired with a g incr
rowType = 'pharm-gincr';

% Find the sweep indices for each row
if strcmpi(rowType, 'pharm')
    % Take only the first trace from each condition
    for iRow = 1:nRows
        swpIndRow{iRow} = [];
        for iG = 1:nGCondToFit
            % 
            indThisCond = swpIndGincrPcond{iG, iPEachRow(iRow)};

            % Select "most representative trace"
            swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                            swpIndHasBursts, maxNoise);
            if isempty(swpIdxSelected)
                error(['Most representative trace not found', ...
                        ' for iG == %d, iP == %d!'], ...
                        iG, iPEachRow(iRow));
            else
                swpIndRow{iRow} = [swpIndRow{iRow}; swpIdxSelected];
            end
        end
    end

    % Number of traces per row is three
    nColumns = 3;
elseif strcmpi(rowType, 'pharm-gincr')
    % Take only the first trace from each condition
    for iRow = 1:nRows
        % Get the current iG and iP
        iGThis = iGEachRow(iRow);
        iPThis = iPEachRow(iRow);

        % Find the indices for this condition
        indThisCond = find(isGCond{iGThis} & isPCond{iPThis} & ...
                            fromCell & toFit);

        % Select "most representative trace"
        swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                        swpIndHasBursts, maxNoise);
        if isempty(swpIdxSelected)
            error(['Most representative trace not found', ...
                    ' for iG == %d, iP == %d!'], ...
                    iGEachRow(iRow), iPEachRow(iRow));
        else
            swpIndRow{iRow} = swpIdxSelected;
        end
    end

    % Number of traces per row is one
    nColumns = 1;
end


for iRow = 1:nRows
    % Get the current iG and iP
    iGThis = iGEachRow(iRow);
    iPThis = iPEachRow(iRow);

    % Find the indices for this condition
    indThisCond = isGCond{iGThis} & isPCond{iPThis} & fromCell & toFit;

    % Select "most representative trace"
    swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                    swpIndHasBursts, maxNoise);
    if isempty(swpIdxSelected)
        error(['Most representative trace not found', ...
                ' for iG == %d, iP == %d!'], ...
                iGEachRow(iRow), iPEachRow(iRow));
    else
        swpIndRow{iRow} = swpIdxSelected;
    end
end

indThisCond = swpIndGincrPcond{iG, iPEachRow(iRow)};

% Select "most representative trace"
swpIdxSelected = select_trace(indThisCond, swpIndHasLts, ...
                                swpIndHasBursts, maxNoise);

% Find all indices for this condition for this cell
nIndThisCond = length(indThisCond);
if nIndThisCond <= 0
    error('There must be traces for this cell & condition!\n');
end

% Find all indices of traces that do not have an LTS
indThisCondHasNoLts = setdiff(indThisCond, swpIndHasLts);
nIndThisCondHasNoLts = length(indThisCondHasNoLts);

% Find all indices of traces that have an LTS
indThisCondHasLts = ...
    intersect(indThisCond, swpIndHasLts, 'sorted');

% Find all indices of traces that have an LTS but not a burst
indThisCondHasLtsOnly = setdiff(indThisCondHasLts, swpIndHasBursts);
nIndThisCondHasLtsOnly = length(indThisCondHasLtsOnly);

% Find all indices of traces that have a burst
indThisCondHasBurst = ...
    intersect(indThisCond, swpIndHasBursts, 'sorted');
nIndThisCondHasBurst = length(indThisCondHasBurst);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%