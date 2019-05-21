function parse_all_multiunit(varargin)
%% Tests the parse_multiunit function on all files in the present working directory
% Usage: parse_all_multiunit(varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
% Arguments:
%       varargin    - 'PlotAllFlag': whether to plot everything
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeDetectionFlag': whether to plot spike detection
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeDensityFlag': whether to plot spike density
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeHistogramFlag': whether to plot spike histograms
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotAutoCorrFlag': whether to plot autocorrelegrams
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotRawFlag': whether to plot raw traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotRasterFlag': whether to plot raster plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotMeasuresFlag': whether to plot time series 
%                                           of measures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/extract_fileparts.m
%       cd/force_matrix.m
%       cd/parse_all_abfs.m
%       cd/parse_multiunit.m
%       cd/plot_measures.m

% File History:
% 2019-03-13 Created
% 2019-03-14 Now combines all files from the same slice
% 2019-04-29 Now saves combined data as a matfile
% 2019-05-06 Added input parser and plot flags
% TODO: Make outFolder optional parameters
% TODO: Make combining optional

%% Hard-coded parameters
inFolder = pwd;
outFolder = pwd;
saveMatFlag = false;
matFileSuffix = '_multiunit_data';
varsNeeded = {'vVecsSl', 'siMsSl', 'iVecsSl', 'sliceBases', 'phaseBoundaries'};

%% Default values for optional arguments
plotAllFlagDefault = false;
plotSpikeDetectionFlagDefault = false;
plotSpikeDensityFlagDefault = false;
plotSpikeHistogramFlagDefault = false;
plotAutoCorrFlagDefault = false;
plotRawFlagDefault = false;
plotRasterFlagDefault = false;
plotMeasuresFlagDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeDetectionFlag', plotSpikeDetectionFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeDensityFlag', plotSpikeDensityFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotSpikeHistogramFlag', plotSpikeHistogramFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotAutoCorrFlag', plotAutoCorrFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotRawFlag', plotRawFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotRasterFlag', plotRasterFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotMeasuresFlag', plotMeasuresFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
plotAllFlag = iP.Results.PlotAllFlag;
plotSpikeDetectionFlag = iP.Results.PlotSpikeDetectionFlag;
plotSpikeDensityFlag = iP.Results.PlotSpikeDensityFlag;
plotSpikeHistogramFlag = iP.Results.PlotSpikeHistogramFlag;
plotAutoCorrFlag = iP.Results.PlotAutoCorrFlag;
plotRawFlag = iP.Results.PlotRawFlag;
plotRasterFlag = iP.Results.PlotRasterFlag;
plotMeasuresFlag = iP.Results.PlotMeasuresFlag;

%% Generate and save data vectors for each slice
inFolderName = extract_fileparts(inFolder, 'dirbase');

% Create a file name for all multi-unit data
matPath = fullfile(outFolder, [inFolderName, matFileSuffix, '.mat']);

% Load or process data for each slice
if isfile(matPath)
    % Load data for each slice
    fprintf("Loading data for each slice ...\n");
    load(matPath, varsNeeded{:});
else
    % Combine data from the same slice
    fprintf("Combining data for each slice ...\n");
    [vVecsSl, siMsSl, iVecsSl, sliceBases, phaseBoundaries] = ...
        combine_data_from_same_slice(inFolder);

    % Save data for each slice
    if saveMatFlag
        save(matPath, varsNeeded{:}, '-v7.3');
    end
end

%% Parse all slices
% Count the number of slices
nSlices = numel(vVecsSl);

% Preallocate parsed parameters and data
muParams = cell(nSlices, 1);
muData = cell(nSlices, 1);

% Parse and plot recordings from each slice
for iSlice = 1:nSlices
    % Print message
    fprintf("Parsing slice #%d ...\n", iSlice);

    % Parse and plot multi-unit recordings from this slice
    [muParams{iSlice}, muData{iSlice}] = ...
        parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                        'PulseVectors', iVecsSl{iSlice}, ...
                        'PlotAllFlag', plotAllFlag, ...
                        'PlotSpikeDetectionFlag', plotSpikeDetectionFlag, ...
                        'PlotSpikeDensityFlag', plotSpikeDensityFlag, ...
                        'PlotSpikeHistogramFlag', plotSpikeHistogramFlag, ...
                        'PlotAutoCorrFlag', plotAutoCorrFlag, ...
                        'PlotRawFlag', plotRawFlag, ...
                        'PlotRasterFlag', plotRasterFlag, ...
                        'PlotMeasuresFlag', plotMeasuresFlag, ...
                        'OutFolder', outFolder, ...
                        'FileBase', sliceBases{iSlice}, ...
                        'PhaseBoundaries', phaseBoundaries{iSlice});

    % Close all figures
    close all force hidden;
end

% for iSlice = 1:nSlices; [muParams{iSlice}, muData{iSlice}] = parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), 'PulseVectors', iVecsSl{iSlice}, 'PlotFlag', plotFlag, 'OutFolder', outFolder, 'FileBase', sliceBases{iSlice}, 'PhaseBoundaries', phaseBoundaries{iSlice}); close all force hidden; end

%% Plot tuning curves for all measures
% plot_measures;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vVecsSl, siMsSl, iVecsSl, sliceBases, phaseBoundaries] = ...
            combine_data_from_same_slice(inFolder, varargin)
%% Combines the data from the same slices
% TODO: What if not all slices have the same number of files?
% TODO: Use *phase* in the file name to detect phase boundaries

%% Hard-coded parameters
regexpSliceName = '.*slice[0-9]*';

%% Preparation
% Get all the abf file names
[~, allAbfPaths] = ...
    all_files('Directory', inFolder, 'Extension', 'abf');

% Extract all slice names
allSliceNames = extract_fileparts(allAbfPaths, 'base', ...
                                    'RegExp', regexpSliceName);

% Get unique slice names
sliceBases = unique(allSliceNames);

% Count the number of slices
nSlices = numel(sliceBases);

%% Do the job
vVecsSl = cell(nSlices, 1);
siMsSl = nan(nSlices, 1);
iVecsSl = cell(nSlices, 1);
phaseBoundaries = cell(nSlices, 1);
parfor iSlice = 1:nSlices
    [vVecsSl{iSlice}, siMsSl(iSlice), iVecsSl{iSlice}, phaseBoundaries{iSlice}] = ...
        combine_data_from_one_slice(inFolder, sliceBases{iSlice}, varargin{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vVecsSl, siMsSl, iVecsSl, phaseBoundaries] = ...
            combine_data_from_one_slice(inFolder, sliceBase, varargin)
%% Combines the data for one slice

%% Hard-coded parameters
regexpPhaseStr = 'phase[a-zA-Z0-9]*';

%% Count files and phases
% Get all .abf files for this slice
[~, allAbfPaths] = ...
    all_files('Directory', inFolder, 'Prefix', sliceBase, 'Extension', 'abf');

% Extract phase strings
allPhaseStrs = extract_fileparts(allAbfPaths, 'base', 'RegExp', regexpPhaseStr);

% Extract phase IDs
allPhaseIDs = extractAfter(allPhaseStrs, 'phase');

% Sort the unique phase IDs
%   Note: This assumes the phase IDs are in correct alphanumeric order
phaseIDs = unique(allPhaseIDs, 'sorted');

% Count the number of phases
nPhases = numel(phaseIDs);

%% Extract data to combine
% Parse all multi-unit recordings for this slice
[allParams, allData] = ...
    parse_all_abfs('FileNames', allAbfPaths, ...
                    'ChannelTypes', {'voltage', 'current'}, ...
                    'ChannelUnits', {'uV', 'arb'}, ...
                    'SaveSheetFlag', false);

% Extract parameters, then clear unused parameters
siMs = allParams.siMs;
clear allParams;

% Extract data, then clear unused data
vVecs = allData.vVecs;
iVecs = allData.iVecs;
clear allData;

%% Order the data correctedly
% Find the indices for each phase
indEachPhase = cellfun(@(x) find_in_strings(x, allAbfPaths), ...
                        phaseIDs, 'UniformOutput', false);

% Find the indices for each phase
sortOrder = vertcat(indEachPhase{:});

% Reorder data
[siMs, vVecs, iVecs] = argfun(@(x) x(sortOrder), siMs, vVecs, iVecs);

%% Combine the data
% Compute the new siMs
siMsSl = mean(siMs);

% Concatenate vectors
[vVecsSl, iVecsSl] = argfun(@force_matrix, vVecs, iVecs);

%% Create phase boundaries
% Count the number of phase boundaries
nBoundaries = nPhases - 1;

% Compute phase boundaries
if nBoundaries == 0
    phaseBoundaries = [];
else
    % Count the number of files for each phase
    nFilesEachPhase = count_samples(indEachPhase);

    % Get the index of the last file for each phase
    iFileLastEachPhase = cumsum(nFilesEachPhase);

    % Get the index of the first file for each phase
    iFileFirstEachPhase = iFileLastEachPhase - nFilesEachPhase + 1;

    % Count the number of sweeps in each file
    nSweepsEachFile = count_vectors(vVecs);

    % Count the number of sweeps in each phase
    nSweepsEachPhase = arrayfun(@(x, y) sum(nSweepsEachFile(x:y)), ...
                            iFileFirstEachPhase, iFileLastEachPhase);

    % Get the index of the last sweep for each phase
    iFileLastEachPhase = cumsum(nSweepsEachPhase);

    % Compute the phase boundaries
    phaseBoundaries = iFileLastEachPhase(1:nBoundaries) + 0.5;
end

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

plotFlag = false;

%% Hard-coded parameters
nFilesPerSlice = 3;

% Count the number of slices
nSlices = ceil(nFiles / nFilesPerSlice);

% Create indices for the slices
indSlices = transpose(1:nSlices);

% Compute the index of the first file
idxFileFirst = arrayfun(@(x) nFilesPerSlice * (x - 1) + 1, indSlices);

% Compute the index of the last file
idxFileLast = arrayfun(@(x) min(nFilesPerSlice * x, nFiles), indSlices);

% Compute the slice bases
sliceBases = arrayfun(@(x, y) extract_fileparts(abfFullFileName(x:y), ...
                            'commonprefix'), idxFileFirst, idxFileLast, ...
                    'UniformOutput', false);

abfFullFileName = allParams.abfFullFileName;

% Count the number of files
nFiles = numel(abfFullFileName);

% Compute the new siMs
siMsSl = arrayfun(@(x, y) mean(siMs(x:y)), idxFileFirst, idxFileLast);

% Concatenate vectors
horzcatcell = @(x) horzcat(x{:});
[vVecsSl, iVecsSl] = ...
    argfun(@(z) arrayfun(@(x, y) horzcatcell(z(x:y)), ...
                    idxFileFirst, idxFileLast, 'UniformOutput', false), ...
            vVecs, iVecs);

% Create phase boundaries if nFilesPerSlice is more than 1
if nBoundaries > 0
    phaseBoundaries = ...
        arrayfun(@(x, y) cumsum(nSweepsEachFile(x:(y - 1))) + 0.5, ...
                    idxFileFirst, idxFileLast, 'UniformOutput', false);
else
    phaseBoundaries = [];
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
