%% m3ha_choose_heterogeneous_population.m
%% Chooses a heterogeneous population based on 2-cell network responses
%% Hard-coded parameters
oscDataDir = fullfile('/media', 'adamX', 'm3ha', 'manuscript', 'figures', 'Figure07');
oscDataFileName = '20200204T1042_using_bestparams_20200203_manual_singleneuronfitting0-102_vtraub_-65_rank2,4-5,7,9-10,12-13,16,20-21,23,25,29_oscillation_params.csv';

pCondDesired = 2;
gIncrDesired = 200;

%% Create full path
oscDataPath = fullfile(oscDataDir, oscDataFileName);

%% Read the table
oscDataTable = readtable(oscDataPath);

%% Extract from the table
cellName = oscDataTable.cellName;
rankNum = oscDataTable.rankNum;
pCond = oscDataTable.pCond;
gIncr = oscDataTable.gIncr;
oscPeriod2Ms = oscDataTable.oscPeriod2Ms;

%% Select rows
rowsSelected = pCond == pCondDesired & round(gIncr*12) == gIncrDesired;
oscPeriod2MsSelected = oscPeriod2Ms(rowsSelected);
cellNameSelected = cellName(rowsSelected);
rankNumSelected = rankNum(rowsSelected);

%% Reorder
[oscPeriod2MsSelected, origInd] = sort(oscPeriod2MsSelected);
cellNameSelected = cellNameSelected(origInd);
rankNumSelected = rankNumSelected(origInd);

%% Print results
cellfun(@(x, y, z) fprintf('%s (rank %d): %g\n', x, y, z), ...
        cellNameSelected, num2cell(rankNumSelected), num2cell(oscPeriod2MsSelected));