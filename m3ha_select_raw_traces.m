function [nColumns, g200p_ind, p_ind, row_ind] = m3ha_select_raw_traces (colmode, rowConditions, attempt_number, fiti, scpgv_ind, indtofit, gg_new, pp, fnrow, bursttime, ltspeaktime, cellstofit_ind, cellstofit_cellname, maxnoise)
%% Select raw traces to import
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Count the number of rows
nRows = size(rowConditions, 1);

% Initialize nColumns
nColumns = 0;

g200p_ind = cell(1, length(pp));
p_ind = cell(1, length(pp));
row_ind = cell(1, nRows);
temp_ind = cell(1, length(pp));
temp_ind1 = cell(1, length(pp));
temp_ind2 = cell(1, length(pp));
temp_ind3 = cell(1, length(pp));
temp2_ind = cell(length(gg_new), length(pp));
goodb_ind = find(bursttime > 0);
goodlts_ind = find(ltspeaktime > 0);
E091710_ind = find(~cellfun(@isempty, strfind(fnrow, 'E091710')));
switch colmode
case 1
    if attempt_number < 3
        fprintf('Fitting across trials of cell E091710 ... \n');
    else
        fprintf('Fitting across trials of cell %s ... \n', ...
                cellstofit_cellname{fiti});
    end
    switch attempt_number
    case 0
        % Attempt #0: Use 4 traces of E091710 @ 200% g_incr
        fprintf('Attempt #0: Use 4 traces of E091710 @ 200p g_incr ... \n');
        for hi = 1:length(pp)        % for each pharmacological condition @ 200% g_incr
            g200p_ind{hi} = ...
                intersect(scpgv_ind(:, :, hi, 4, :), E091710_ind, 'sorted');
            if isempty(g200p_ind{hi})
                error('Does not have enough data for E091710!')
            end
        end
        nColumns = 1;
    case 1
        % Attempt #1: Use all traces of E091710 @ 200% g_incr
        fprintf('Attempt #1: Use all traces of E091710 @ 200p g_incr ... \n');
        for hi = 1:length(pp)        % for each pharmacological condition @ 200% g_incr
            g200p_ind{hi} = ...
                intersect(scpgv_ind(:, :, hi, 4, :), E091710_ind, 'sorted');
            if isempty(g200p_ind{hi})
                error('Does not have enough data for E091710!')
            end
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, g200p_ind));
    case 2
        % Attempt #2: Use all traces of E091710
        fprintf('Attempt #2: Use all traces of E091710 ... \n');
        for hi = 1:length(pp)        % for each pharmacological condition
            p_ind{hi} = ...
                intersect(scpgv_ind(:, :, hi, :, :), E091710_ind, 'sorted');
            if isempty(p_ind{hi})
                error('Does not have enough data for E091710!')
            end
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, p_ind));
    case 3
        % Attempt #3: Use all traces of each cell to fit @ 100%, 200%, 400% g_incr
        fprintf('Attempt #3: Use all traces of ');
        fprintf('%s @ 100p, 200p, 400p g_incr ... \n', ...
                    cellstofit_cellname{fiti});
        if size(rowConditions, 2) == 1        % each row is one pharm condition
            for rowi = 1:nRows
                row_ind{rowi} = [];
                for gi = 1:length(gg_new)
                    row_ind{rowi} = [row_ind{rowi}; ...
                                cellstofit_ind{fiti}{gi, rowConditions(rowi, 1)}];
                end
            end
        elseif size(rowConditions, 2) == 2    % each row is a pharm condition paired with a g incr
            for rowi = 1:nRows
                row_ind{rowi} = ...
                    cellstofit_ind{fiti}{rowConditions(rowi, 2), rowConditions(rowi, 1)};
            end
        end

        % Make sure each row has the same number of traces
        nColumns = min(cellfun(@length, row_ind));
    case 4
        % Attempt #4: Use 12 traces of each cell to fit 
        %   (one "best trial" for each of 4 pharm conditions X 100%, 200%, 400% g_incr)
        fprintf('Attempt #4: Use 12 traces (one "best trial" for each condition)');
        fprintf(' of %s @ 100p, 200p, 400p g_incr ... \n', ...
                    cellstofit_cellname{fiti});
        if size(rowConditions, 2) == 1        % each row is one pharm condition
            % Take only the first trace from each condition
            for rowi = 1:nRows
                row_ind{rowi} = [];
                for gi = 1:length(gg_new)
                    % Select "most representative trace"
                    idxSelected = select_trace(cellstofit_ind{fiti}, ...
                                                gi, rowConditions(rowi, 1), ...
                                                goodlts_ind, goodb_ind, ...
                                                maxnoise);
                    if isempty(idxSelected)
                        error(['Most representative trace not found', ...
                                ' for gi == %d, hi == %d!'], ...
                                gi, rowConditions(rowi, 1));
                    else
                        row_ind{rowi} = [row_ind{rowi}; idxSelected];
                    end
                end
            end
            % Number of traces per row
            nColumns = 3;
        elseif size(rowConditions, 2) == 2    % each row is a pharm condition paired with a g incr
            % Take only the first trace from each condition
            for rowi = 1:nRows
                % Select "most representative trace"
                idxSelected = select_trace(cellstofit_ind{fiti}, ...
                                            rowConditions(rowi, 2), ...
                                            rowConditions(rowi, 1), ...
                                            goodlts_ind, goodb_ind, maxnoise);
                if isempty(idxSelected)
                    error(['Most representative trace not found', ...
                            ' for gi == %d, hi == %d!'], ...
                            rowConditions(rowi, 2), rowConditions(rowi, 1));
                else
                    row_ind{rowi} = idxSelected;
                end
            end
            % Number of traces per row
            nColumns = 1;
        end
    end
case 2        
    fprintf('Fitting across cells ... \n');
    switch attempt_number
    case 1
        % Attempt #1: Find cells with LTSs present for all pharm conditions @ 200% g_incr
        fprintf('Attempt #1: Find cells with LTSs present ');
        fprintf('for all pharm conditions @ 200p g_incr ... \n');
        ct = 0;                 % counts cells with an lts for all pharm conditions
        for ci = 1:length(cc)
            to_use_cell = true;
            for hi = 1:length(pp)       % for each pharmacological condition @ 200% g_incr
                temp_ind{hi} = ...
                    intersect(scpgv_ind(:, ci, hi, 4, :), goodlts_ind, 'sorted');
                if isempty(temp_ind{hi})
                    to_use_cell = false;
                    break;
                end
            end
            if to_use_cell == true 
                ct = ct + 1;
                for hi = 1:length(pp)   % for each pharmacological condition @ 200% g_incr
                    g200p_ind{hi}(ct) = temp_ind{hi}(1);
                end
            end
        end
        nColumns = ct;                     % In our data set, nColumns == 13
    case 2
        % Attempt #2: Find cells with bursts present for all pharm conditions @ 200% g_incr
        fprintf('Attempt #2: Find cells with bursts present ');
        fprintf('for all pharm conditions @ 200p g_incr ... \n');
        ct = 0;                % counts cells with an lts for all pharm conditions
        for ci = 1:length(cc)
            to_use_cell = true;
            for hi = 1:length(pp)       % for each pharmacological condition @ 200% g_incr
                temp_ind{hi} = ...
                    intersect(scpgv_ind(:, ci, hi, 4, :), goodb_ind, 'sorted');
                if isempty(temp_ind{hi})
                    to_use_cell = false;
                    break;
                end
            end
            if to_use_cell == true 
                ct = ct + 1;
                for hi = 1:length(pp)   % for each pharmacological condition @ 200% g_incr
                    g200p_ind{hi}(ct) = temp_ind{hi}(1);
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
        for ci = 1:length(cc)
            to_use_cell = true;
            ngoodb = 0;
            for hi = 1:length(pp)       % for each pharmacological condition @ 200% g_incr
                temp_ind1{hi} = ...
                    intersect(scpgv_ind(:, ci, hi, 4, :), ...
                                goodb_ind, 'sorted');
                temp_ind2{hi} = ...
                    intersect(scpgv_ind(:, ci, hi, 4, :), ...
                                goodlts_ind, 'sorted');
                temp_ind3{hi} = ...
                    setdiff(scpgv_ind(:, ci, hi, 4, :), 0);
                if ~isempty(temp_ind1{hi})
                    ngoodb = ngoodb + 1;
                end
                if isempty(temp_ind3{hi})
                    to_use_cell = false;
                end
            end
            if ngoodb < 3
                to_use_cell = false;
            end
            if to_use_cell == true
                ct = ct + 1;
                for hi = 1:length(pp)   % for each pharmacological condition @ 200% g_incr
                    if ~isempty(temp_ind1{hi})
                        g200p_ind{hi}(ct) = temp_ind1{hi}(1);
                    elseif ~isempty(temp_ind2{hi})
                        g200p_ind{hi}(ct) = temp_ind2{hi}(1);
                    else
                        g200p_ind{hi}(ct) = temp_ind3{hi}(1);
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
        for ci = 1:length(cc)
            to_use_cell = true;
            for hi = 1:length(pp)       % for each pharmacological condition
                for gi = 1:length(gg_new)    % for each g_incr
                    temp2_ind{gi, hi} = ...
                        intersect(scpgv_ind(:, ci, hi, gi + 2, :), ...
                                    indtofit, 'sorted');
                    if isempty(temp2_ind{gi, hi})
                        to_use_cell = false;
                        break;
                    end
                end
            end
            if to_use_cell == true 
                ct = ct + 1;
                rowi = 0;
                for hi = 1:length(pp)   % for each pharmacological condition
                    for gi = 1:length(gg_new)    % for each g incr
                        rowi = rowi + 1;
                        row_ind{rowi}(ct) = temp2_ind{gi, hi}(1);
                    end
                end
            end
        end
        nColumns = ct;                % In our data set, nColumns == 36
    case 5
        % Attempt #5: Choose the "most representative" sweep from each cell in fiti
        fprintf(['Attempt #5: Choose the "most representative" sweep for each cell ', ...
                    'in cells #%d to #%d ... \n'], fiti(1), fiti(end));
        % Get number of cells to fit
        ncellstofit = length(fiti);

        % Take a trace from each cell for each condition
        %   with priority given to a trace with bursts,
        %   then to a trace with LTSs
        if size(rowConditions, 2) == 1        % each row is one pharm condition              
            % Number of traces per row is 3 times number of cells to fit
            nColumns = 3*ncellstofit;

            for rowi = 1:nRows          % for each pharm condition
                % Initialize row_ind
                row_ind{rowi} = [];     

                % Build row_ind
                for gi = 1:length(gg_new)   % for each g incr
                    for iFit = fiti             % for each cell to fit
                        % Select "most representative trace"
                        idxSelected = select_trace(cellstofit_ind{iFit}, ...
                                                    gi, ...
                                                    rowConditions(rowi, 1), ...
                                                    goodlts_ind, goodb_ind, ...
                                                    maxnoise);
                        row_ind{rowi} = [row_ind{rowi}; idxSelected];
                    end
                end

                % Verify length of row_ind{rowi}
                if length(row_ind{rowi}) ~= nColumns
                    error(['Not enough traces selected for row #%d! ', ...
                            'Something is wrong!\n\n'], rowi);
                end
            end            
        elseif size(rowConditions, 2) == 2    % each row is a pharm condition paired with a g incr
            % Number of traces per row equals the number of cells to fit
            nColumns = ncellstofit;

            % Take only the first trace from each cell for each condition
            for rowi = 1:nRows          % for each pharm condition
                % Initialize row_ind
                row_ind{rowi} = [];     

                % Build row_ind
                for iFit = fiti             % for each cell to fit
                    % Select "most representative trace"
                    idxSelected = select_trace(cellstofit_ind{iFit}, ...
                                                rowConditions(rowi, 2), ...
                                                rowConditions(rowi, 1), ...
                                                goodlts_ind, goodb_ind, ...
                                                maxnoise);
                    row_ind{rowi} = [row_ind{rowi}; idxSelected];
                end

                % Verify length of row_ind{rowi}
                if length(row_ind{rowi}) ~= nColumns
                    error(['Not enough traces selected for row #%d! ', ...
                            'Something is wrong!\n\n'], rowi);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idxSelected = select_trace(allIndices, gi, hi, goodlts_ind, goodb_ind, maxnoise)

% Find all indices for this condition for this cell
indThisCond = allIndices{gi, hi};
nIndThisCond = length(indThisCond);
if nIndThisCond <= 0
    error('There must be traces for this cell & condition!\n');
end

% Find all indices of traces that do not have an LTS
indThisCondHasNoLts = setdiff(indThisCond, goodlts_ind);
nIndThisCondHasNoLts = length(indThisCondHasNoLts);

% Find all indices of traces that have an LTS
indThisCondHasLts = ...
    intersect(indThisCond, goodlts_ind, 'sorted');

% Find all indices of traces that have an LTS but not a burst
indThisCondHasLtsOnly = setdiff(indThisCondHasLts, goodb_ind);
nIndThisCondHasLtsOnly = length(indThisCondHasLtsOnly);

% Find all indices of traces that have a burst
indThisCondHasBurst = ...
    intersect(indThisCond, goodb_ind, 'sorted');
nIndThisCondHasBurst = length(indThisCondHasBurst);

% Select one "most representative" trace
if nIndThisCondHasNoLts > nIndThisCond / 2
    % Look for the trial with minimum noise in all trials with no LTS
    [~, iTemp] = min(maxnoise(indThisCondHasNoLts));
    idxSelected = indThisCondHasNoLts(iTemp);
elseif nIndThisCondHasLtsOnly > nIndThisCondHasBurst
    % Look for the trial with minimum noise in all trials with LTS but no burst
    [~, iTemp] = min(maxnoise(indThisCondHasLtsOnly));
    idxSelected = indThisCondHasLtsOnly(iTemp);
elseif nIndThisCondHasBurst >= nIndThisCondHasLtsOnly
    % Look for the trial with minimum noise in all trials with burst
    [~, iTemp] = min(maxnoise(indThisCondHasBurst));
    idxSelected = indThisCondHasBurst(iTemp);
else
    % Something is wrong 
    error('Error in LTS or burst data??\n');   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%        ng200p_ind(hi) = ct;
%        nColumns = min(min(ng200p_ind), 20);        
% number of cells to take per pharmacological condition
%                                
% TO FIX: doesn't work if too many cells are chosen

nColumns = min(nColumns, 5);

row_ind{rowi} = [row_ind{rowi}; ...
            cellstofit_ind{fiti}{gi, rowConditions(rowi, 1)}(1)];

% Select one trace with the priority given to 
%   bursts, then to LTSs
if ~isempty(indThisCondHasBurst)
    idxSelected = indThisCondHasBurst(1);
elseif ~isempty(indThisCondHasLts)
    idxSelected = indThisCondHasLts(1);
else
    idxSelected = indThisCond(1);
end

% Find all indices for this cell
indThisCell = cellstofit_ind{iFit};

% Find all indices for this condition for this cell
indThisCond = indThisCell{rowConditions(rowi, 2), rowConditions(rowi, 1)};

% Find all indices of traces that have an LTS
indThisCondHasLts = ...
    intersect(indThisCond, goodlts_ind, 'sorted');

% Find all indices of traces that have a burst
indThisCondHasBurst = ...
    intersect(indThisCond, goodb_ind, 'sorted');

% Select one trace with the priority given to 
%   bursts, then LTSs
if ~isempty(indThisCondHasBurst)
    idxSelected = indThisCondHasBurst(1);
elseif ~isempty(indThisCondHasLts)
    idxSelected = indThisCondHasLts(1);
else
    idxSelected = indThisCond(1);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%