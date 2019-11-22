function [numcells, numsets, numswps, cellnames, abffullfn, nswpsCPV, ...
                nswpsUsed, fnrow, cellidrow, prow, vrow, grow, swpnrow, ...
                gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, ...
                gabab_w, actIhold] = ...
                m3ha_CountSweeps (infolder, ljp, debugflag, debug_files, ...
                    preallocateflag, maxcells, maxsets, maxswps)
%% Counts total number of usable cells, sets and sweeps, generate filenames and record sweep properties
% Usage: [numcells, numsets, numswps, cellnames, abffullfn, nswpsCPV, ...
%               nswpsUsed, fnrow, cellidrow, prow, vrow, grow, swpnrow, ...
%               gabab_amp, gabab_Trise, gabab_TfallFast, gabab_TfallSlow, ...
%               gabab_w, actIhold] = ...
%               m3ha_CountSweeps (infolder, ljp, debugflag, debug_files, ...
%                   preallocateflag, maxcells, maxsets, maxswps)
% Arguments:
%       infolder    - 
%       ljp         - liquid junction potential used to correct voltage trace (mV)
%       preallocateflag    - (opt) whether maxcells, maxsets & maxswps are available
%       maxcells    - (opt) Total number of neurons recorded
%       maxsets     - (opt) Total number of cell-pharm-vhold conditions recorded
%       maxswps     - (opt) Total number of sweeps recorded
%
% Used by:
%       cd/m3ha_dclampDataExtractor.m
%

% File History:
% 2016-11-07 Moved from dclampDataExtractor.m
% 2016-11-30 E092810_0000 was added to broken_files and taken out of all analyses 
%        (otherwise the GAT1 block might be overrepresented statistically)
% 2016-12-13 Reversed sign of LJP
% 2018-08-06 Added actIhold
%

%% Fixed parameters used in the experiments
maxnsetspc = 12;                % A maximum of 12 sets per cell (4 pharm x 3 Vhold)
nswps_cpvg = 5;                 % # of sweeps per CPVG condition
pharm = [1; 2; 3; 4];           % Possible pharm conditions (1 - Control; 2 - GAT1 Block; 3 - GAT3 Block; 4 - Dual Block)
Vhold = [-50; -55; -60];        % Possible Vhold values (as shown on PClamp, not LJP-corrected)

%% Specify experiment names and conditions
% fprefix1 - dclamp experiments with "n-series" dyn files - 
% 25%, 50%, 100%, 200%, 400% gabab ipsc conductance increments
fprefix{1} = {'091710'; '091810'; '092110'; '092210'; '092710'; '092810'; '092910'};     % dates for the experimental series
ffol{1} = {'091710_dclampnew1'; '091810_dclampnew2'; '092110_dclampnew3'; '092210_dclampnew4'; '092710_dclampnew5'; '092810_dclampnew6'; '092910_dclampnew7'};                                     % directory (folder) names for the experimental series
fg{1} = [25, 50, 100, 200, 400];                                % G incr (%) values for the experimental series
% fprefix2 - dclamp experiments with "q-series" dyn files - 
% 50%, 100%, 200%, 400% gabab ipsc conductance increments
fprefix{2} = {'092910_7.5'; '100110'; '100810'};
ffol{2} = {'092910_dclampnew7.5'; '100110_dclampnew8'; '100810_dclampnew9'}; 
fg{2} = [50, 100, 200, 400];
% fprefix3 - dclamp experiments with "r-series" dyn files - 
% 100%, 200%, 400% gabab ipsc conductance increments
fprefix{3} = {'101210'};
ffol{3} = {'101210_dclampnew10'};
fg{3} = [100, 200, 400];
% fprefix4 - dclamp experiments with "t-series" dyn files - 
% 100%, 200%, 400%, 800% gabab ipsc conductance increments
fprefix{4} = {'101310'};
ffol{4} = {'101310_dclampnew11'};
fg{4} = [100, 200, 400, 800];

%% Problematic abf files
broken_files = {'B100110_0000', ...     % Cannot open with abf2load
            'E100110_0001', ...         % Only one sweep was recorded
            'G091810_0000', ...         % Channels were mixed up for at least one trace
            'E092810_0000', ...         % A GAT1 block curve was applied instead of a GAT3 block
            };                          %     but E092810_0002 is already the GAT1 block
            
%% Template parameters used for GABAB IPSC conductance waveforms
% Note: amp in Christine's thesis was actually for 200 % G incr
gtemp_amp = [16.00; 24.00; 8.88; 6.32];            % (nS)
gtemp_Trise = [52.00; 52.00; 38.63; 39.88];        % (ms)
gtemp_TfallFast = [90.10; 90.10; 273.40; 65.80];    % (ms)
gtemp_TfallSlow = [1073.20; 1073.20; 1022.00; 2600.00];    % (ms)
gtemp_w = [0.952; 0.952; 0.775; 0.629];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize counters
celln = 0;      % number of cells counted so far
setc = 0;       % counts number of sets
swpc = 0;       % counts number of sweeps

% Initialize variables if maxcells, maxsets & maxswps are known
if preallocateflag
    % Initialize variables different for each cell
    cellnames = cell(1, maxcells);

    % Initialize variables different for each set
    abffullfn = cell(1, maxsets);       % abf file name for this set
    nswpsCPV = zeros(1, maxsets);      % total number of sweeps for this set
    nswpsUsed = zeros(1, maxsets);     % total number of sweeps before this set
    holdingCurrents = cell(1, maxsets); % all holding currents for this set

    % Initialize variables different for each sweep
    fnrow = cell(1, maxswps);
    cellidrow = zeros(1, maxswps);
    prow = zeros(1, maxswps);
    vrow = zeros(1, maxswps);
    grow = zeros(1, maxswps);
    swpnrow = zeros(1, maxswps);
    gabab_amp = zeros(1, maxswps);
    gabab_Trise = zeros(1, maxswps);
    gabab_TfallFast = zeros(1, maxswps);
    gabab_TfallSlow = zeros(1, maxswps);
    gabab_w = zeros(1, maxswps);
    actIhold = zeros(1, maxswps);
end

% FOR each experimental series (each type of G incr (%) series)
for f = 1:4
    % Set up tag tables
    ptag = repmat(pharm, 1, length(Vhold));             % pharm conditions for each set
    vtag = repmat((Vhold' - ljp), length(pharm), 1);    % Vhold values for each set
    gtag = repmat(fg{f}, 1, nswps_cpvg);                % G incr (%) for each sweep
    swpntag = reshape(repmat(1:nswps_cpvg, numel(fg{f}), 1), 1, []);    % Sweep # for each sweep

    % FOR each day of experiment
    for q = 1:numel(fprefix{f})            
        % Load data from one day's experiments, 
        %     using "reorg_xxxxxx.xls" excel files to figure out which ones to load
        thisprefix = cell2mat(fprefix{f}(q));
        thisfolder = cell2mat(ffol{f}(q));
        excelfn = ['reorg_', thisprefix, '.xls'];
        dex0 = importdata(fullfile(infolder, excelfn));
        if isfield(dex0, 'textdata') == 1       % When does this happen?
            dex = dex0.textdata.Sheet1;
            dex(:,1) = [];
        elseif isfield(dex0, 'Sheet1') == 1     % This is what happens as far as I can tell
            dex = dex0.Sheet1;
        end
        col1 = dex(:,1);            % This contains the cell names

        % FOR each row in "reorg_xxxxxx.xls"
        for j = 1:length(col1)
            if col1{j} > 0                      % Found a new cell!
                % Increase cell count and record cell name
                celln = celln + 1;
                cellnames{celln} = col1{j}(6);  % This is a letter of the alphabet
                flist = dex((2:5) + j, 3:5);    % Contains the set # for each pharm-Vhold pair

                % FOR each possible set (each pharm-Vhold pair)
                for sn = 1:maxnsetspc            
                    % Get the content of this cell 
                    %   (a string in either the form '#### (## -> ##)' 
                    %                       or the form '#### (##)')
                    thisStr = flist{sn};

                    % Determine whether to load data
                    if isempty(thisStr)       % Skip this set if no data exist
                        continue;
                    end

                    % Get the set number
                    setn = thisStr(1:4);      % set #

                    if strcmp(thisprefix, '092910_7.5') == 1
                        datfn = [cellnames{celln}, '092910_', setn];
                    else
                        datfn = [cellnames{celln}, thisprefix, '_', setn];
                    end
                    if ismember(datfn, broken_files)    % Skip this set if abf file cannot be parsed
                        continue;
                    end
                    if debugflag && ~ismember(datfn, debug_files)
                        continue;
                    end

                    % Increase set count and create .abf file name
                    setc = setc + 1;
                    abffullfn{setc} = fullfile(infolder, thisfolder, [datfn, '.abf']);

                    % Get the holding currents in pA
                    [holdingCurrentsCPV, nHoldingCurrents] = ...
                        sscanf_full(thisStr(5:end), '%f');

                    % Store the holding currents
                    holdingCurrents{setc} = holdingCurrentsCPV;

                    % Check if there is a holding current
                    if nHoldingCurrents < 1
                        error('No holding current recorded for %s!!', datfn);
                    end

                    % Count the total number of sweeps for this set
                    nswpsCPV(setc) = length(gtag);

                    % Count the total number of sweeps before this set
                    nswpsUsed(setc) = swpc;

                    % Assign holding currents to each sweep in this set
                    actIholdThisSet = piecelinspace(holdingCurrentsCPV', ...
                                                    nswpsCPV(setc));

                    % Record sweep properties
                    fprintf('Recording sweep properties ...\n');
                    for iSwp = 1:nswpsCPV(setc)         % FOR each sweep
                        filebase = [datfn, '_', num2str(iSwp)];
                        swpc = swpc + 1;

                        % Record relevant sweep properties
                        fnrow{swpc} = [filebase, '.mat'];
                        cellidrow(swpc) = celln;        % Cell ID #
                        prow(swpc) = ptag(sn);          % Pharm condition
                        vrow(swpc) = vtag(sn);          % Vhold, ljp-corrected (mV)
                        grow(swpc) = gtag(iSwp);        % GABAB IPSC G incr (%)
                        swpnrow(swpc) = swpntag(iSwp);  % Within condition sweep #

                        % Record GABAB IPSC conductance waveform parameters for this sweep
                        gabab_amp(swpc) = (grow(swpc)/100) * gtemp_amp(prow(swpc));
                        gabab_Trise(swpc) = gtemp_Trise(prow(swpc));
                        gabab_TfallFast(swpc) = gtemp_TfallFast(prow(swpc));
                        gabab_TfallSlow(swpc) = gtemp_TfallSlow(prow(swpc));
                        gabab_w(swpc) = gtemp_w(prow(swpc));

                        % Average out the holding currents
                        %   based on the sweep number
                        actIhold(swpc) = actIholdThisSet(iSwp);
                    end
                end
            end
        end
    end
end

% Record total number of usable cells, sets and sweeps
numcells = celln;
numsets = setc;
numswps = swpc;

