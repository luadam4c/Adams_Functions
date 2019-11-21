function varargout = parse_all_multiunit(varargin)
%% Runs the parse_multiunit function on all slices in the present working directory
% Usage: [muParams, muData] = parse_all_multiunit(varargin)
% Explanation:
%       TODO
%
% Example(s):
%       parse_all_multiunit;
%       apply_to_all_subdirs(@parse_all_multiunit);
%       [muParams, muData] = parse_all_multiunit;
%       parse_all_multiunit('PlotRaw', true);
%       parse_all_multiunit('PlotAll', true);
%
% Outputs:
%       muParams    - parsed multiunit parameters
%       muData      - parsed multiunit data
%
% Arguments:
%       varargin    - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SliceBases': names of slices to analyze
%                   must be a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from directory
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'OutFolder': directory to place output files
%                   must be a string scalar or a character vector
%                   default == same as inFolder
%                   - 'PlotMeasuresFlag': whether to plot time series 
%                                           of measures
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the parse_multiunit() function
%
% Requires:
%       cd/all_files.m
%       cd/all_slice_bases.m
%       cd/combine_multiunit_data.m
%       cd/parse_multiunit.m
%
% Used by:
%       cd/clc2_analyze.m

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
% 2019-08-06 Now accepts any parameter-value pair for parse_multiunit.m
% 2019-08-24 Now uses varargout
% TODO: Load data one file at a time

%% Hard-coded parameters
% TODO Make optional arguments
saveMatFlag = true;
% TODO: matFileSuffix = '_multiunit_data';?
varsNeeded = {'sliceBase', 'vVecsSl', 'siMsSl', 'iVecsSl', ...
                'phaseBoundaries', 'phaseStrs'};
regexpSliceMatFile = '.*slice[0-9]*.mat';
regexpSliceAbfFile = '.*slice[0-9]*.*.abf';
channelTypes = {'voltage', 'current'};
channelUnits = {'uV', 'arb'};

%% Default values for optional arguments
directoryDefault = pwd;
sliceBasesDefault = {};                 % detect from directory by default
inFolderDefault = '';                   % set later
outFolderDefault = '';                  % set later
plotMeasuresFlagDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SliceBases', sliceBasesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PlotMeasuresFlag', plotMeasuresFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
sliceBases = iP.Results.SliceBases;
inFolder = iP.Results.InFolder;
outFolder = iP.Results.OutFolder;
plotMeasuresFlag = iP.Results.PlotMeasuresFlag;

% Keep unmatched arguments for the parse_multiunit() function
otherArguments = struct2arglist(iP.Unmatched);

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


% Either combine data from .abf files or load .mat files
if isempty(sliceBasesMat)
    % Combine data from the same slice possibly without saving 
    %   the combined data as .mat files
    fprintf('Combining data for each slice ...\n');
    allDataTable = ...
        combine_multiunit_data('Directory', inFolder, ...
                                'SliceBase', sliceBases, ...
                                'ChannelTypes', channelTypes, ...
                                'ChannelUnits', channelUnits, ...
                                'SaveMatFlag', saveMatFlag, ...
                                'VarsToSave', varsNeeded);
else
    % See if there are any slice data yet to be combined
    sliceToCombine = setdiff(sliceBasesAbf, sliceBasesMat);

    % Combine data from each slice that has not been combined
    cellfun(@(x) combine_multiunit_data('SliceBase', x, ...
                                    'ChannelTypes', channelTypes, ...
                                    'ChannelUnits', channelUnits, ...
                                    'SaveMatFlag', true, ...
                                    'VarsToSave', varsNeeded), ...
            sliceToCombine);

    % Get all the slice data .mat file names available
    [~, allMatPaths] = ...
        all_files('Directory', inFolder, 'RegExp', regexpSliceMatFile, ...
                    'SortBy', 'date', 'ForceCellOutput', true);

    % Restricted to specific slices if provided
    if ~isempty(sliceBases)
        allMatPaths = allMatPaths(contains(allMatPaths, sliceBases));
    end

    % Load data for each slice as a structure array
    fprintf('Loading data for each slice ...\n');
    allDataStruct = cellfun(@(x) load(x, varsNeeded{:}), allMatPaths);

    % Convert to a table
    allDataTable = struct2table(allDataStruct, 'AsArray', true);
end

% Don't proceed if no files found
if isempty(allDataTable)
    return
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
    if nargout >= 2
        [muParams{iSlice}, muData{iSlice}] = ...
            parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                            'PulseVectors', iVecsSl{iSlice}, ...
                            'FileBase', sliceBases{iSlice}, ...
                            'PhaseBoundaries', phaseBoundaries{iSlice}, ...
                            'OutFolder', outFolder, ...
                            'PlotMeasuresFlag', plotMeasuresFlag, ...
                            otherArguments{:});
    elseif nargout >= 1
        muParams{iSlice} = ...
            parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                            'PulseVectors', iVecsSl{iSlice}, ...
                            'FileBase', sliceBases{iSlice}, ...
                            'PhaseBoundaries', phaseBoundaries{iSlice}, ...
                            'OutFolder', outFolder, ...
                            'PlotMeasuresFlag', plotMeasuresFlag, ...
                            otherArguments{:});
    else
        parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), ...
                        'PulseVectors', iVecsSl{iSlice}, ...
                        'FileBase', sliceBases{iSlice}, ...
                        'PhaseBoundaries', phaseBoundaries{iSlice}, ...
                        'OutFolder', outFolder, ...
                        'PlotMeasuresFlag', plotMeasuresFlag, ...
                        otherArguments{:});
    end

    % Close all figures
    close all force hidden;
end

%% Output results
varargout{1} = muParams;
if nargout >= 2
    varargout{2} = muData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%% Plot tuning curves for all measures
if plotMeasuresFlag
    plot_measures;
end

for iSlice = 1:nSlices; [muParams{iSlice}, muData{iSlice}] = parse_multiunit(vVecsSl{iSlice}, siMsSl(iSlice), 'PulseVectors', iVecsSl{iSlice}, 'PlotFlag', plotFlag, 'OutFolder', outFolder, 'FileBase', sliceBases{iSlice}, 'PhaseBoundaries', phaseBoundaries{iSlice}); close all force hidden; end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
