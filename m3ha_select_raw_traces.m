function [figurePositions, swpIndices, fileNames] = ...
                m3ha_select_raw_traces (columnMode, rowConditions, ...
                    attemptNumber, iCellToFit, swpIdxSCPGV, swpIndToFit, ...
                    swpIndByCondition, cellNamesToFit, swpInfo)
%% Select raw traces to import
% Usage: [figurePositions, swpIndices, fileNames] = ...
%               m3ha_select_raw_traces (columnMode, rowConditions, ...
%                   attemptNumber, iCellToFit, swpIdxSCPGV, swpIndToFit, ...
%                   swpIndByCondition, cellNamesToFit, swpInfo)
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
% 2017-08-21 Added maxnoise as an argument
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation and code clarity

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
% Count the number of possible conductance amplitude scaling percentages
nGIncr = length(gIncrAll);

% Count the number of possible pharm conditions 
nPCond = length(pCondAll);

% Count the number of rows
nRows = size(rowConditions, 1);

% Determine the row type
%   Note: Must be consistent with m3ha_determine_row_conditions.m
if size(rowConditions, 2) == 1
    % Each row is one pharm condition
    rowType = 'pharm';

    % Each row is the iPCond
    iPEachRow = rowConditions;
elseif size(rowConditions, 2) == 2
    % Each row is a pharm condition paired with a g incr
    rowType = 'pharm-gincr';

    % Extract the iPCond for each row
    iPEachRow = rowConditions(:, 1);

    % Extract the iGincr for each row
    iGEachRow = rowConditions(:, 2);
else
    error('This rowConditions is not recognized!');
end

% Initialize nColumns
nColumns = 0;

% Extract from swpInfo
fnrow = swpInfo.fnrow;
bursttime = swpInfo.bursttime;
ltspeaktime = swpInfo.ltspeaktime;
maxnoise = swpInfo.maxnoise;

%% Do the job
% Print message
fprintf('Selecting the traces for this cell to fit ... \n');

% Find all sweep indices with bursts
swpIndHasBursts = find(bursttime > 0);

% Find all sweep indices with LTSs
swpIndHasLts = find(ltspeaktime > 0);

% Find all sweeps for the cell E091710
swpIndE091710 = find(~cellfun(@isempty, strfind(fnrow, 'E091710')));

% Initialize various schemes
swpIndG200P = cell(1, nPCond);
swpIndPCond = cell(1, nPCond);
swpIndRow = cell(1, nRows);
swpIndTemp1 = cell(1, nPCond);
swpIndTemp2 = cell(1, nPCond);
swpIndTemp3 = cell(1, nPCond);
swpIndTemp4 = cell(1, nPCond);
swpIndTemp5 = cell(nGIncr, nPCond);

% Decide on the sweeps to use for each condition
switch columnMode
case 1
    % Get the current cell name
    cellName = cellNamesToFit{iCellToFit};

    % Get the sweep indices for each g incr-pharm condition
    swpIndGincrPcond = swpIndByCondition{iCellToFit};

    if attemptNumber < 3
        fprintf('Fitting across trials of cell E091710 ... \n');
    else
        fprintf('Fitting across trials of cell %s ... \n', ...
                cellName);
    end
    switch attemptNumber
    case 0
        % Attempt #0: Use 4 traces of E091710 @ 200% g_incr
        fprintf('Attempt #0: Use 4 traces of E091710 @ 200p g_incr ... \n');
        for iP = 1:nPCond        % for each pharmacological condition @ 200% g_incr
            swpIndG200P{iP} = ...
                intersect(swpIdxSCPGV(:, :, iP, 4, :), swpIndE091710, 'sorted');
            if isempty(swpIndG200P{iP})
                error('Does not have enough data for E091710!')
            end
        end
        nColumns = 1;
    case 1
        % Attempt #1: Use all traces of E091710 @ 200% g_incr
        fprintf('Attempt #1: Use all traces of E091710 @ 200p g_incr ... \n');
        for iP = 1:nPCond        % for each pharmacological condition @ 200% g_incr
            swpIndG200P{iP} = ...
                intersect(swpIdxSCPGV(:, :, iP, 4, :), swpIndE091710, 'sorted');
            if isempty(swpIndG200P{iP})
                error('Does not have enough data for E091710!')
            end
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, swpIndG200P));
    case 2
        % Attempt #2: Use all traces of E091710
        fprintf('Attempt #2: Use all traces of E091710 ... \n');
        for iP = 1:nPCond        % for each pharmacological condition
            swpIndPCond{iP} = ...
                intersect(swpIdxSCPGV(:, :, iP, :, :), swpIndE091710, 'sorted');
            if isempty(swpIndPCond{iP})
                error('Does not have enough data for E091710!')
            end
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, swpIndPCond));
    case 3
        % Attempt #3: Use all traces of each cell to fit @ 100%, 200%, 400% g_incr
        fprintf(['Attempt #3: Use all traces of ', ...
                    '%s @ 100p, 200p, 400p g_incr ... \n'], cellName);

        % Use all data already filtered
        if strcmpi(rowType, 'pharm')
            swpIndRow = arrayfun(@(x) reshape([swpIndGincrPcond{:, x}], [], 1), ...
                                transpose(1:nRows), 'UniformOutput', false);

        elseif strcmpi(rowType, 'pharm-gincr') 
            swpIndRow = arrayfun(@(x, y) swpIndGincrPcond{x, y}, ...
                                iGEachRow, iPEachRow, 'UniformOutput', false);
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, swpIndRow));
    case 4
        % Attempt #4: Use 12 traces of each cell to fit 
        %   (one "best trial" for each of 4 pharm conditions X 100%, 200%, 400% g_incr)
        fprintf('Attempt #4: Use 12 traces (one "best trial" for each condition)');
        fprintf(' of %s @ 100p, 200p, 400p g_incr ... \n', cellName);
        if strcmpi(rowType, 'pharm')
            % Take only the first trace from each condition
            for iRow = 1:nRows
                swpIndRow{iRow} = [];
                for iG = 1:nGIncr
                    % Select "most representative trace"
                    swpIdxSelected = select_trace(swpIndGincrPcond, ...
                                                iG, iPEachRow(iRow), ...
                                                swpIndHasLts, swpIndHasBursts, ...
                                                maxnoise);
                    if isempty(swpIdxSelected)
                        error(['Most representative trace not found', ...
                                ' for iG == %d, iP == %d!'], ...
                                iG, iPEachRow(iRow));
                    else
                        swpIndRow{iRow} = [swpIndRow{iRow}; swpIdxSelected];
                    end
                end
            end
            % Number of traces per row
            nColumns = 3;
        elseif strcmpi(rowType, 'pharm-gincr')
            % Take only the first trace from each condition
            for iRow = 1:nRows
                % Select "most representative trace"
                swpIdxSelected = select_trace(swpIndGincrPcond, ...
                                            iGEachRow(iRow), ...
                                            iPEachRow(iRow), ...
                                            swpIndHasLts, swpIndHasBursts, maxnoise);
                if isempty(swpIdxSelected)
                    error(['Most representative trace not found', ...
                            ' for iG == %d, iP == %d!'], ...
                            iGEachRow(iRow), iPEachRow(iRow));
                else
                    swpIndRow{iRow} = swpIdxSelected;
                end
            end
            % Number of traces per row
            nColumns = 1;
        end
    end
case 2        
    fprintf('Fitting across cells ... \n');
    switch attemptNumber
    case 1
        % Attempt #1: Find cells with LTSs present for all pharm conditions @ 200% g_incr
        fprintf('Attempt #1: Find cells with LTSs present ');
        fprintf('for all pharm conditions @ 200p g_incr ... \n');
        ct = 0;                 % counts cells with an lts for all pharm conditions
        for iC = 1:length(cellIdAll)
            toUse = true;
            for iP = 1:nPCond       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasLts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCond   % for each pharmacological condition @ 200% g_incr
                    swpIndG200P{iP}(ct) = swpIndTemp1{iP}(1);
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
            for iP = 1:nPCond       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasBursts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCond   % for each pharmacological condition @ 200% g_incr
                    swpIndG200P{iP}(ct) = swpIndTemp1{iP}(1);
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
            for iP = 1:nPCond       % for each pharmacological condition @ 200% g_incr
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
                for iP = 1:nPCond   % for each pharmacological condition @ 200% g_incr
                    if ~isempty(swpIndTemp2{iP})
                        swpIndG200P{iP}(ct) = swpIndTemp2{iP}(1);
                    elseif ~isempty(swpIndTemp3{iP})
                        swpIndG200P{iP}(ct) = swpIndTemp3{iP}(1);
                    else
                        swpIndG200P{iP}(ct) = swpIndTemp4{iP}(1);
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
            for iP = 1:nPCond       % for each pharmacological condition
                for iG = 1:nGIncr    % for each g_incr
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
                for iP = 1:nPCond   % for each pharmacological condition
                    for iG = 1:nGIncr    % for each g incr
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
                for iG = 1:nGIncr   % for each g incr
                    for iFit = iCellToFit             % for each cell to fit
                        % Select "most representative trace"
                        swpIdxSelected = select_trace(swpIndByCondition{iFit}, ...
                                                    iG, ...
                                                    iPEachRow(iRow), ...
                                                    swpIndHasLts, swpIndHasBursts, ...
                                                    maxnoise);
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
                    % Select "most representative trace"
                    swpIdxSelected = select_trace(swpIndByCondition{iFit}, ...
                                                iGEachRow(iRow), ...
                                                iPEachRow(iRow), ...
                                                swpIndHasLts, swpIndHasBursts, ...
                                                maxnoise);
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
if nColumns == 0
    error('No traces found!')
end

% Count the number of traces
nSweeps = nRows * nColumns;

% Get the figure positions and sweep indices for each trace
figurePositions = cell(nSweeps, 1);
swpIndices = zeros(nSweeps, 1);
ct = 0;
for iRow = 1:nRows
    for iCol = 1:nColumns
        % Get the current figure position
        figurePositions{ct} = [iRow, iCol];

        % Get the current sweep index
        if ~isempty(swpIndG200P{1})
            swpIndices(ct) = swpIndG200P{iRow}(iCol);
        elseif ~isempty(swpIndPCond{1})
            swpIndices(ct) = swpIndPCond{iRow}(iCol);
        elseif ~isempty(swpIndRow{1})
            swpIndices(ct) = swpIndRow{iRow}(iCol);
        end
    end
end

% Get the file names
fileNames = fnrow(swpIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function swpIdxSelected = select_trace(allIndices, iG, iP, swpIndHasLts, swpIndHasBursts, maxnoise)

% Find all indices for this condition for this cell
indThisCond = allIndices{iG, iP};
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

% Select one "most representative" trace
if nIndThisCondHasNoLts > nIndThisCond / 2
    % Look for the trial with minimum noise in all trials with no LTS
    [~, iTemp] = min(maxnoise(indThisCondHasNoLts));
    swpIdxSelected = indThisCondHasNoLts(iTemp);
elseif nIndThisCondHasLtsOnly > nIndThisCondHasBurst
    % Look for the trial with minimum noise in all trials with LTS but no burst
    [~, iTemp] = min(maxnoise(indThisCondHasLtsOnly));
    swpIdxSelected = indThisCondHasLtsOnly(iTemp);
elseif nIndThisCondHasBurst >= nIndThisCondHasLtsOnly
    % Look for the trial with minimum noise in all trials with burst
    [~, iTemp] = min(maxnoise(indThisCondHasBurst));
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
    for iG = 1:nGIncr
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

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%