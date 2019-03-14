%% Tests the parse_multiunit function on all files in the present working directory

% File History:
% 2019-03-13 Created

% Parse all abfs
[allParams, allData] = parse_all_abfs('ChannelTypes', {'voltage', 'current'}, 'ChannelUnits', {'mV', 'arb'});

%%
siMs = allParams.siMs;
abfFullFileName = allParams.abfFullFileName;

tVecs = allData.tVec;
vVecs = allData.vVecs;
iVecs = allData.iVecs;

fileBases = extract_fileparts(abfFullFileName, 'base');

%% 
nFiles = numel(siMs);

%%
muParams = cell(nFiles, 1);
muData = cell(nFiles, 1);

%%
for iFile = 1:nFiles
%parfor iFile = 1:nFiles
%for iFile = 1    
    [muParams{iFile}, muData{iFile}] = parse_multiunit(vVecs{iFile}, siMs(iFile), 'PulseVectors', iVecs{iFile}, 'tVecs', tVecs{iFile}, 'PlotFlag', true, 'OutFolder', pwd, 'FileBase', fileBases{iFile});
    close all force hidden;
end
