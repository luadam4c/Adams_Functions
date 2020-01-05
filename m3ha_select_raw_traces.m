function [fileNames, rowConditions, figurePositions] = ...
                m3ha_select_raw_traces (cellNameOrId, varargin)
%% Select raw traces to import for specific cells
% Usage: [fileNames, rowConditions, figurePositions] = ...
%               m3ha_select_raw_traces (cellNameOrId, varargin)
%
% Explanation: 
%       TODO
%
% Examples: 
%       [fileNames, rowCond, figPos] = m3ha_select_raw_traces('D101310');
%
% Outputs:
%       TODO
%
% Arguments:
%       cellNameOrId- cell name or cell ID in cellInfo table
%                   must be a character vector or a string scalar
%                       or a numeric scalar
%       varargin    - 'ColumnMode': column mode
%                   must be empty or one of:
%                       1 - across trials
%                       2 - across cells TODO
%                   default == 1
%                   - 'RowMode': row mode
%                   must be empty or one of:
%                       1 - each row is a pharm condition
%                       2 - each row is a pharm, g incr pair
%                   default == 1 if columnMode == 1
%                              2 if columnMode == 2
%                   - 'AttemptNumber': attempt number
%                   must be empty or one of:
%                       FOR columnMode == 1
%                           1 - Use 4 traces @ 200% gIncr for this data mode
%                           2 - Use all traces @ 200% gIncr for this data mode
%                           3 - Use all traces for this data mode
%                           4 - Use 1 trace for each pharm x gIncr 
%                                   for this data mode
%                           5 - Use 4 traces @ 400% gIncr for this data mode
%                       FOR columnMode == 2
%                           1 - Find cells with LTSs present 
%                                   for all pharm conditions @ 200% g_incr
%                           2 - Find cells with bursts present 
%                                   for all pharm conditions @ 200% g_incr
%                           3 - Find cells with bursts present 
%                                   for 3 out of 4 pharm conditions @ 200% g_incr, 
%                                   and choose the "best" sweep for each cell
%                           4 - Find cells within indices to fit 
%                                   present for all pharm conditions and all g_incr
%                           5 - Choose the "best" sweep from each cell in iCellToFit
%                   default == 4 if columnMode == 1
%                              5 if columnMode == 2
%                   - 'DataMode': data mode
%                   must be empty or one of:
%                       0 - all data (no restrictions)
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                       3 - all data 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 0
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
%                                       of sweeps to fit
%                   must be a positive integer vector, a string array 
%                       or a cell array of character vectors
%                   default == set in m3ha_select_sweeps.m
%                   - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == set in m3ha_select_sweeps.m
%
% Requires:
%       cd/istext.m
%       cd/m3ha_create_cell_info_table.m
%       cd/m3ha_determine_row_conditions.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_select_sweeps.m
%
% Used by:
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_rank_neurons.m
%       cd/singleneuronfitting42.m and later versions

% File History:
% 2017-05-20 Moved from singleneuronfitting2.m
% 2017-05-22 Changed line width and indentation
% 2017-08-11 Added Attempt #5 for fitting across cells
% 2017-08-21 Modified Attempt #4 for fitting across trials so that
%               the "most representative trial" is selected
% 2017-08-21 Added maxNoise as an argument
% 2018-11-15 Moved to Adams_Functions
% 2018-11-15 Improved documentation and code clarity
% 2018-12-06 Now uses cellIdThis, fileNamesToUse, swpInfo, cellInfo
% 2019-12-17 Now selects the "best trace" for Attempt #1 across trials
% 2019-12-21 Now uses input parser
% 2019-12-23 Fixed default dataMode to be 0 so that the toUse column
%               already set in swpInfo would be respected

%% Hard-coded parameters
pharmStr = 'prow';
gIncrStr = 'grow';
burstDelayStr = 'bursttime';
ltsDelayStr = 'ltspeaktime';
maxNoiseStr = 'maxnoise';
toUseStr = 'toUse';
cellNameStr = 'cellName';

%% Default values for optional arguments
columnModeDefault = [];
rowModeDefault = [];
attemptNumberDefault = [];
dataModeDefault = [];
swpInfoDefault = table.empty;
cellInfoDefault = table.empty;
rowsToUseDefault = [];
casesDirDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'cellNameOrId', ...
    @(x) validateattributes(x, {'char', 'string', 'cell', 'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ColumnMode', columnModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RowMode', rowModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'AttemptNumber', attemptNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
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
parse(iP, cellNameOrId, varargin{:});
columnMode = iP.Results.ColumnMode;
rowMode = iP.Results.RowMode;
attemptNumber = iP.Results.AttemptNumber;
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

% Set default data mode
if isempty(dataMode)
    dataMode = 0;
end

% Update swpInfo so that there is a toUse column
swpInfo = m3ha_select_sweeps('SwpInfo', swpInfo, 'RowsToUse', rowsToUse, ...
                                'DataMode', dataMode, 'CasesDir', casesDir);

% Decide on the cell name(s)
if istext(cellNameOrId)
    cellName = cellNameOrId;
elseif isnumeric(cellNameOrId)
    % Generate a table of cell names from swpInfo if not provided
    if isempty(cellInfo)
        cellInfo = m3ha_create_cell_info_table('SwpInfo', swpInfo);
    end

    % Extract correct cell name(s)
    cellName = cellInfo{cellNameOrId, cellNameStr};
else
    error('cellNameOrId unrecognized!');
end

% Extract from swpInfo
fnrow = swpInfo.Properties.RowNames;
grow = swpInfo{:, gIncrStr};
prow = swpInfo{:, pharmStr};
burstTime = swpInfo{:, burstDelayStr};
ltsPeakTime = swpInfo{:, ltsDelayStr};
maxNoise = swpInfo{:, maxNoiseStr};
toUse = swpInfo{:, toUseStr};

% Get all unique conductance amplitude scaling percentages to fit
gCondToUse = unique(swpInfo{toUse, gIncrStr});

% Get all unique pharmacological conditions to fit
pCondToUse = unique(swpInfo{toUse, pharmStr});

% Count the number of distinct conductance amplitude scaling percentages
nGCondToUse = length(gCondToUse);

% Count the number of distinct pharmacological conditions 
nPCondToUse = length(pCondToUse);

% Set default column mode
if isempty(columnMode)
    columnMode = 1;
end

% Set default row mode
if isempty(rowMode)
    if columnMode == 1
        rowMode = 1;
    else
        rowMode = 2;
    end
end

% Set default attempt number
if isempty(attemptNumber)
    if columnMode == 1
        attemptNumber = 4;
    else
        attemptNumber = 5;
    end
end

%% Determine row conditions
% Update rowMode based on colMode and attemptNumber
if rowMode ~= 1 && ...
    ((colMode == 1 && attemptNumber <= 2) || ...
    (colMode == 2 && attemptNumber <= 3))
    fprintf(['Row Mode changed to 1 for ', ...
            'column mode %d and attempt number %d!\n'], ...
            colMode, attemptNumber);
    rowMode = 1;
end

% Determine what each row would be
rowConditions = m3ha_determine_row_conditions(rowMode, pCondToUse, gCondToUse);

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
isPCond = arrayfun(@(x) prow == x, pCondToUse, 'UniformOutput', false);

% Determine whether each sweep is of each gIncr condition
%   Note: isGCond is a cell array of logical vectors
isGCond = arrayfun(@(x) grow == x, gCondToUse, 'UniformOutput', false);

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
swpIndTemp1 = cell(1, nPCondToUse);
swpIndTemp2 = cell(1, nPCondToUse);
swpIndTemp3 = cell(1, nPCondToUse);
swpIndTemp4 = cell(1, nPCondToUse);
swpIndTemp5 = cell(nGCondToUse, nPCondToUse);

% Decide on the sweeps to use for each condition
switch columnMode
case 1
    if iscell(cellName)
        cellName = cellName{1};
    end

    % Print message
    fprintf('Fitting across trials of cell %s ... \n', cellName);

    % Decide on the indices based on attempt number
    switch attemptNumber
    case {1, 2, 5}
        % Print message
        if attemptNumber == 1
            fprintf('Attempt #1: Using 4 traces of %s @ 200%% gIncr ... \n', ...
                    cellName);
        elseif attemptNumber == 2
            fprintf('Attempt #2: Using all traces of %s @ 200%% gIncr ... \n', ...
                    cellName);
        elseif attemptNumber == 5
            fprintf('Attempt #5: Using 4 traces of %s @ 400%% gIncr ... \n', ...
                    cellName);
        end

        % Return whether gCondToUse is 200
        isGIncr200 = gCondToUse == 200;

        % Return whether gCondToUse is 400
        isGIncr400 = gCondToUse == 400;

        % Find the sweep indices for each pharmacological condition 
        %   @ 200% g_incr from this cell to be fitted
        if attemptNumber == 1
            swpIndRow = cellfun(@(x) select_trace(x & isGCond{isGIncr200} & ...
                                                    fromCell & toUse, ...
                                                hasLts, hasBurst, maxNoise), ...
                                    isPCond, 'UniformOutput', false);
        elseif attemptNumber == 2
            swpIndRow = cellfun(@(x) find(x & isGCond{isGIncr200} & ...
                                            fromCell & toUse), ...
                                    isPCond, 'UniformOutput', false);
        elseif attemptNumber == 5
            swpIndRow = cellfun(@(x) select_trace(x & isGCond{isGIncr400} & ...
                                                    fromCell & toUse, ...
                                                hasLts, hasBurst, maxNoise), ...
                                    isPCond, 'UniformOutput', false);
        end

        % Check if there is enough data or not
        if any(cellfun(@isempty, swpIndRow))
            error('There is not enough data for %s @ 200%% gIncr!', cellName);
        end

        % Decide on the number of columns per row
        if attemptNumber == 1 || attemptNumber == 5
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
        swpIndRow = cellfun(@(x) find(x & fromCell & toUse), ...
                            isRowCond, 'UniformOutput', false);

        % Check if there is enough data or not
        if any(cellfun(@isempty, swpIndRow))
            error('There is not enough data for %s!', cellName);
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, swpIndRow));
    case 4
        % Print message
        fprintf('Attempt #4: Using one "best trial" for each ');
        fprintf('condition) of %s ... \n', cellName);

        % Take only the most representative trace from each condition
        if size(rowConditions, 2) == 1
            for iRow = 1:nRows
                % Get the corresponding iPCond
                iPThis = iPEachRow(iRow);

                % Initialize
                swpIndThisRow = [];

                % Get three traces
                for iG = 1:nGCondToUse
                    % Determine whether each sweep is this condition
                    isThisCond = isGCond{iG} & isPCond{iPThis} & ...
                                    fromCell & toUse;

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
            swpIndRow = arrayfun(@(x) select_trace(x & fromCell & toUse, ...
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
            for iP = 1:nPCondToUse       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasLts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCondToUse   % for each pharmacological condition @ 200% g_incr
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
            for iP = 1:nPCondToUse       % for each pharmacological condition @ 200% g_incr
                swpIndTemp1{iP} = ...
                    intersect(swpIdxSCPGV(:, iC, iP, 4, :), swpIndHasBursts, 'sorted');
                if isempty(swpIndTemp1{iP})
                    toUse = false;
                    break;
                end
            end
            if toUse
                ct = ct + 1;
                for iP = 1:nPCondToUse   % for each pharmacological condition @ 200% g_incr
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
            for iP = 1:nPCondToUse       % for each pharmacological condition @ 200% g_incr
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
                for iP = 1:nPCondToUse   % for each pharmacological condition @ 200% g_incr
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
            for iP = 1:nPCondToUse       % for each pharmacological condition
                for iG = 1:nGCondToUse    % for each g_incr
                    swpIndTemp5{iG, iP} = ...
                        intersect(swpIdxSCPGV(:, iC, iP, iG + 2, :), ...
                                    swpIndToUse, 'sorted');
                    if isempty(swpIndTemp5{iG, iP})
                        toUse = false;
                        break;
                    end
                end
            end
            if toUse 
                ct = ct + 1;
                iRow = 0;
                for iP = 1:nPCondToUse   % for each pharmacological condition
                    for iG = 1:nGCondToUse    % for each g incr
                        iRow = iRow + 1;
                        swpIndRow{iRow}(ct) = swpIndTemp5{iG, iP}(1);
                    end
                end
            end
        end
        nColumns = ct;                % In our data set, nColumns == 36
    case 5
        % Attempt #5: Choose the "most representative" sweep from each cell in iCellToUse
        fprintf(['Attempt #5: Choose the "most representative" sweep for each cell ', ...
                    'in cells #%d to #%d ... \n'], iCellToUse(1), iCellToUse(end));
        % Get number of cells to fit
        ncellstofit = length(iCellToUse);

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
                for iG = 1:nGCondToUse   % for each g incr
                    for iFit = iCellToUse             % for each cell to fit

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
                for iFit = iCellToUse             % for each cell to fit
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

function swpIdxSelected = select_trace (isThisCond, hasLts, hasBursts, maxNoise)
%% Selects the most representative trace based on majority rule

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

% Make sure cell name is not a cell array
if iscell(cellName) && numel(cellName) == 1
    cellName = cellName{1};
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%