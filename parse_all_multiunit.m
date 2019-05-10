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
% 2019-05-06 Added input parser and plot flags
% TODO: Make outFolder, plotFlag optional parameters
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
plotSpikeHistogramFlag = iP.Results.PlotSpikeHistogramFlag;
plotAutoCorrFlag = iP.Results.PlotAutoCorrFlag;
plotRawFlag = iP.Results.PlotRawFlag;
plotRasterFlag = iP.Results.PlotRasterFlag;
plotMeasuresFlag = iP.Results.PlotMeasuresFlag;

%% Preparation
if plotAllFlag
    % TODO: Simplify with argfun.m

    if ~plotSpikeDetectionFlag
        plotSpikeDetectionFlag = true;
    end
    if ~plotSpikeHistogramFlag
        plotSpikeHistogramFlag = true;
    end
    if ~plotAutoCorrFlag
        plotAutoCorrFlag = true;
    end
    if ~plotRawFlag
        plotRawFlag = true;
    end
    if ~plotRasterFlag
        plotRasterFlag = true;
    end
    if ~plotMeasuresFlag
        plotMeasuresFlag = true;
    end
end

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
                        'PlotSpikeDetectionFlag', plotSpikeDetectionFlag, ...
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
%% Combines data from the same slice
% TODO: What if not all slices have the same number of files?

%% Hard-coded parameters
nFilesPerSlice = 3;

%% Parse all abfs
[allParams, allData] = ...
    parse_all_abfs('Directory', inFolder, ...
                    'ChannelTypes', {'voltage', 'current'}, ...
                    'ChannelUnits', {'uV', 'arb'});

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

plotFlag = false;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
