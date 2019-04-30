function parse_all_multiunit(varargin)
%% Tests the parse_multiunit function on all files in the present working directory
%
% Requires:
%       cd/argfun.m
%       cd/count_vectors.m
%       cd/extract_fileparts.m
%       cd/parse_all_abfs.m
%       cd/parse_multiunit.m
%       cd/plot_measures.m

% File History:
% 2019-03-13 Created
% 2019-03-14 Now combines all files from the same slice
% 2019-04-29 Now saves combined data as a matfile
% TODO: Make outFolder, plotFlag optional parameters
% TODO: Make combining optional

%% Hard-coded parameters
inFolder = pwd;
outFolder = pwd;
plotFlag = false;
matFileSuffix = '_multiunit_data';
varsNeeded = {'vVecsSl', 'siMsSl', 'iVecsSl', 'sliceBases', 'phaseBoundaries'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate and save data vectors for each slice
inFolderName = extract_fileparts(inFolder, 'dirbase');

% Create a file name for all multi-unit data
matPath = fullfile(outFolder, [inFolderName, matFileSuffix, '.mat']);

% Load or process data for each slice
if isfile(matPath)
    % Load data for each slice
    load(matPath, varsNeeded{:});
else
    % Combine data from the same slice
    [vVecsSl, siMsSl, iVecsSl, sliceBases, phaseBoundaries] = ...
        combine_data_from_same_slice(inFolder);

    % Save data for each slice
    save(matPath, varsNeeded{:}, '-v7.3');
end

%% Parse all slices
% Preallocate parsed parameters and data
muParams = cell(nSlices, 1);
muData = cell(nSlices, 1);

% Parse and plot recordings from each slice
for iSlice = 1:nSlices
    % Parse and plot multi-unit recordings from this slice
    [muParams{iSlice}, muData{iSlice}] = ...
        parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                        'PulseVectors', iVecsSl{iSlice}, ...
                        'PlotFlag', plotFlag, 'OutFolder', outFolder, ...
                        'FileBase', sliceBases{iSlice}, ...
                        'PhaseBoundaries', phaseBoundaries{iSlice});

    % Close all figures
    close all force hidden;
end

% for iSlice = 1:nSlices; [muParams{iSlice}, muData{iSlice}] = parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), 'PulseVectors', iVecsSl{iSlice}, 'PlotFlag', plotFlag, 'OutFolder', outFolder, 'FileBase', sliceBases{iSlice}, 'PhaseBoundaries', phaseBoundaries{iSlice}); close all force hidden; end

%% Plot tuning curves for all measures
plot_measures;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vVecsSl, siMsSl, iVecsSl, sliceBases, phaseBoundaries] = ...
            combine_data_from_same_slice(inFolder, varargin)
%% Combines data from the same slice
% TODO: What if not all slices have the same number of files?

%% Hard-coded parameters
nFilesPerSlice = 3;

%% Parse all abfs
[allParams, allData] = ...
    parse_all_abfs('Directory', inFolder, ...
                    'ChannelTypes', {'voltage', 'current'}, ...
                    'ChannelUnits', {'mV', 'arb'});

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
