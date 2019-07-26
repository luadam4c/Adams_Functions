function [muParams, muData] = parse_all_multiunit(varargin)
%% Tests the parse_multiunit function on all files in the present working directory
% Usage: [muParams, muData] = parse_all_multiunit(varargin)
% Explanation:
%       TODO
% Example(s):
%       parse_all_multiunit;
%       apply_to_all_subdirs(@parse_all_multiunit);
%       [muParams, muData] = parse_all_multiunit;
%       parse_all_multiunit('PlotRaw', true);
%       parse_all_multiunit('PlotAll', true);
% Outputs:
%       muParams    - parsed multiunit parameters
%       muData      - parsed multiunit data
% Arguments:
%       varargin    - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'OutFolder': directory to place output files
%                   must be a string scalar or a character vector
%                   default == same as inFolder
%                   - 'PlotAllFlag': whether to plot everything
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotCombinedFlag': whether to plot raw data, 
%                           spike density and oscillation duration together
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'PlotSpikeDetectionFlag': whether to plot spike detection
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotSpikeDensityFlag': whether to plot spike density
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotSpikeHistogramFlag': whether to plot spike histograms
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotAutoCorrFlag': whether to plot autocorrelegrams
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotRawFlag': whether to plot raw traces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotRasterFlag': whether to plot raster plots
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'PlotMeasuresFlag': whether to plot time series 
%                                           of measures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == set in parse_multiunit.m
%                   - 'SaveMatFlag': whether to save combined data
%                                           as matfiles
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SaveResultsFlag': whether to save parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/all_files.m
%       cd/all_slice_bases.m
%       cd/combine_data_from_same_slice.m
%       cd/parse_multiunit.m
%       cd/plot_measures.m

% File History:
% 2019-03-13 Created
% 2019-03-14 Now combines all files from the same slice
% 2019-04-29 Now saves combined data as a matfile
% 2019-05-06 Added input parser and plot flags
% 2019-05-21 Now allows each slice to have different numbers of files
% 2019-05-21 Now uses 'slice' or 'phase' in the file name 
%               to detect sliceBase and phase boundaries
% 2019-05-30 Now saves combined vectors as mat files
% 2019-05-31 Updated plot flags
% 2019-06-10 Added 'PlotCombinedFlag'
% 2019-07-24 Added 'Directory', 'InFolder' & 'OutFolder'
% 2019-07-24 Added saveResultsFlag
% 2019-07-25 Now combines slice data if .abf files present 
%               but .mat file not present

%% Hard-coded parameters
matFileSuffix = '_multiunit_data';
varsNeeded = {'sliceBase', 'vVecsSl', 'siMsSl', 'iVecsSl', ...
                'phaseBoundaries', 'phaseStrs'};
regexpSliceMatFile = '.*slice[0-9]*.mat';
regexpSliceAbfFile = '.*slice[0-9]*.*.abf';

%% Default values for optional arguments
directoryDefault = pwd;
inFolderDefault = '';                   % set later
outFolderDefault = '';                  % set later
plotAllFlagDefault = false;
plotCombinedFlagDefault = false;
plotSpikeDetectionFlagDefault = [];     % set in parse_multiunit.m
plotSpikeDensityFlagDefault = [];       % set in parse_multiunit.m
plotSpikeHistogramFlagDefault = [];     % set in parse_multiunit.m
plotAutoCorrFlagDefault = [];           % set in parse_multiunit.m
plotRawFlagDefault = [];                % set in parse_multiunit.m
plotRasterFlagDefault = [];             % set in parse_multiunit.m
plotMeasuresFlagDefault = [];           % set in parse_multiunit.m
saveMatFlagDefault = true;
saveResultsFlagDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotAllFlag', plotAllFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PlotCombinedFlag', plotCombinedFlagDefault, ...
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
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveResultsFlag', saveResultsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
plotAllFlag = iP.Results.PlotAllFlag;
plotCombinedFlag = iP.Results.PlotCombinedFlag;
plotSpikeDetectionFlag = iP.Results.PlotSpikeDetectionFlag;
plotSpikeDensityFlag = iP.Results.PlotSpikeDensityFlag;
plotSpikeHistogramFlag = iP.Results.PlotSpikeHistogramFlag;
plotAutoCorrFlag = iP.Results.PlotAutoCorrFlag;
plotRawFlag = iP.Results.PlotRawFlag;
plotRasterFlag = iP.Results.PlotRasterFlag;
plotMeasuresFlag = iP.Results.PlotMeasuresFlag;
saveMatFlag = iP.Results.SaveMatFlag;
saveResultsFlag = iP.Results.SaveResultsFlag;

%% Preparation
% Decide on the input directory
if isempty(inFolder)
    inFolder = directory;
end

% Decide on the output directory
if isempty(outFolder)
    outFolder = inFolder;
end

% Get all the slice bases with .mat data
sliceBasesMat = all_slice_bases('Directory', inFolder, ...
                                'RegExpFile', regexpSliceMatFile, ...
                                'ForceCellOutput', true, 'SortBy', 'date', ...
                                'RegExpBase', '.*slice[0-9]*');

% Get all the slice bases with .abf data
sliceBasesAbf = all_slice_bases('Directory', inFolder, ...
                                'RegExpFile', regexpSliceAbfFile, ...
                                'ForceCellOutput', true, 'SortBy', 'date', ...
                                'RegExpBase', '.*slice[0-9]*');


% Either load .mat files or combine data from .abf files
if isempty(sliceBasesMat)
    % Combine data from the same slice
    fprintf('Combining data for each slice ...\n');
    allDataTable = ...
        combine_data_from_same_slice('Directory', inFolder, ...
                                    'SaveMatFlag', saveMatFlag, ...
                                    'VarsToSave', varsNeeded);
else
    % See if there are any slice data yet to be combined
    sliceToCombine = setdiff(sliceBasesAbf, sliceBasesMat);

    % Combine data from each slice that has not been combined
    cellfun(@(x) combine_data_from_same_slice('SliceBase', x, ...
                                            'SaveMatFlag', true, ...
                                            'VarsToSave', varsNeeded), ...
            sliceToCombine);

    % Get all the slice data .mat file names available
    [~, allMatPaths] = ...
        all_files('Directory', inFolder, 'RegExp', regexpSliceMatFile, ...
                    'SortBy', 'date', 'ForceCellOutput', true);

    % Load data for each slice as a structure array
    fprintf('Loading data for each slice ...\n');
    allDataStruct = cellfun(@(x) load(x, varsNeeded{:}), allMatPaths);

    % Convert to a table
    allDataTable = struct2table(allDataStruct, 'AsArray', true);
end

% Extract from the table
sliceBases = allDataTable.sliceBase;
vVecsSl = allDataTable.vVecsSl;
siMsSl = allDataTable.siMsSl;
iVecsSl = allDataTable.iVecsSl;
phaseBoundaries = allDataTable.phaseBoundaries;
phaseStrs = allDataTable.phaseStrs;

% Make sure phaseBoundaries is a cell array
if isnumeric(phaseBoundaries)
    phaseBoundaries = num2cell(phaseBoundaries);
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
                        'FileBase', sliceBases{iSlice}, ...
                        'PhaseBoundaries', phaseBoundaries{iSlice}, ...
                        'PlotAllFlag', plotAllFlag, ...
                        'PlotCombinedFlag', plotCombinedFlag, ...
                        'PlotSpikeDetectionFlag', plotSpikeDetectionFlag, ...
                        'PlotSpikeDensityFlag', plotSpikeDensityFlag, ...
                        'PlotSpikeHistogramFlag', plotSpikeHistogramFlag, ...
                        'PlotAutoCorrFlag', plotAutoCorrFlag, ...
                        'PlotRawFlag', plotRawFlag, ...
                        'PlotRasterFlag', plotRasterFlag, ...
                        'PlotMeasuresFlag', plotMeasuresFlag, ...
                        'SaveResultsFlag', saveResultsFlag, ...
                        'OutFolder', outFolder);

    % Close all figures
    close all force hidden;
end

% for iSlice = 1:nSlices; [muParams{iSlice}, muData{iSlice}] = parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), 'PulseVectors', iVecsSl{iSlice}, 'PlotFlag', plotFlag, 'OutFolder', outFolder, 'FileBase', sliceBases{iSlice}, 'PhaseBoundaries', phaseBoundaries{iSlice}); close all force hidden; end

%% Plot tuning curves for all measures
if plotMeasuresFlag
    plot_measures;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
