% parse_all_multiunit.m
%% Tests the parse_multiunit function on all files in the present working directory
%
% Requires:
%       cd/argfun.m
%       cd/extract_fileparts.m
%       cd/parse_all_abfs.m
%       cd/parse_multiunit.m

% File History:
% 2019-03-13 Created
% 2019-03-14 Now combines all files from the same slice
% TODO: Make combining optional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hard-coded parameters
outFolder = pwd;
plotFlag = true; %false; %true;
nFilesPerSlice = 3;

%% Parse all abfs
[allParams, allData] = parse_all_abfs('ChannelTypes', {'voltage', 'current'}, 'ChannelUnits', {'mV', 'arb'});

%% Extract parameters and clear unused parameters
siMs = allParams.siMs;
abfFullFileName = allParams.abfFullFileName;
clear allParams;

% Count the number of files
nFiles = numel(abfFullFileName);

%% Extract data and clear unused data
vVecs = allData.vVecs;
iVecs = allData.iVecs;
clear allData;

% Count the number of vectors in each file
nVectorsEachFile = count_vectors(vVecs);

%% Combine all files from the same slice
% Count the number of slices
nSlices = ceil(nFiles / nFilesPerSlice);

% Create indices for the slices
indSlices = transpose(1:nSlices);

% Compute the index of the first file
idxFileFirst = arrayfun(@(x) nFilesPerSlice * (x - 1) + 1, indSlices);

% Compute the index of the last file
idxFileLast = arrayfun(@(x) min(nFilesPerSlice * x, nFiles), indSlices);

% Compute the new siMs
siMsSl = arrayfun(@(x, y) mean(siMs(x:y)), idxFileFirst, idxFileLast);

% Compute the slice bases
sliceBases = arrayfun(@(x, y) extract_fileparts(abfFullFileName(x:y), ...
                            'commonprefix'), idxFileFirst, idxFileLast, ...
                    'UniformOutput', false);

% Count the number of phase boundaries
nBoundaries = nFilesPerSlice - 1;

% Create phase boundaries if nFilesPerSlice is more than 1
if nBoundaries > 0
    phaseBoundaries = ...
        arrayfun(@(x, y) cumsum(nVectorsEachFile(x:(y - 1))) + 0.5, ...
                    idxFileFirst, idxFileLast, 'UniformOutput', false);
else
    phaseBoundaries = [];
end

% Concatenate vectors
horzcatcell = @(x) horzcat(x{:});
[vVecsSl, iVecsSl] = ...
    argfun(@(z) arrayfun(@(x, y) horzcatcell(z(x:y)), ...
                    idxFileFirst, idxFileLast, 'UniformOutput', false), ...
            vVecs, iVecs);

%% Parse all slices
% Preallocate parsed parameters and data
muParams = cell(nSlices, 1);
muData = cell(nSlices, 1);

% 
for iSlice = 1:nSlices
    [muParams{iSlice}, muData{iSlice}] = ...
        parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                        'PulseVectors', iVecsSl{iSlice}, ...
                        'PlotFlag', plotFlag, 'OutFolder', outFolder, ...
                        'FileBase', sliceBases{iSlice}, ...
                        'PhaseBoundaries', phaseBoundaries{iSlice});

    close all force hidden;
end

%% Plot tuning curves for all measures
clc2_plot_measures;
% for iSlice = 1:nSlices; [muParams{iSlice}, muData{iSlice}] = parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), 'PulseVectors', iVecsSl{iSlice}, 'PlotFlag', plotFlag, 'OutFolder', outFolder, 'FileBase', sliceBases{iSlice}, 'PhaseBoundaries', phaseBoundaries{iSlice}); close all force hidden; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

tVecs = allData.tVec;

% Extract file bases
fileBases = extract_fileparts(abfFullFileName, 'base');

muParams = cell(nFiles, 1);
muData = cell(nFiles, 1);
for iFile = 1:nFiles
    [muParams{iFile}, muData{iFile}] = ...
        parse_multiunit(vVecs{iFile}, siMs(iFile), ...
                        'PulseVectors', iVecs{iFile}, 'tVecs', tVecs{iFile}, ...
                        'PlotFlag', plotFlag, 'OutFolder', outFolder, ...
                        'FileBase', fileBases{iFile});

    close all force hidden;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%