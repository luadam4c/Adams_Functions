function allData = combine_data_from_same_slice (varargin)
%% Combines data across multiple .abf files for each slice in the input folder (or for a particular slice)
% Usage: allData = combine_data_from_same_slice (varargin)
% Explanation:
%       TODO
%       Note: Current and voltage vectors are identified using 
%               identify_channels.m by default. If it's already labelled
%               correctly in the abf files, set 'UseOriginal' to be true.
%
% Example(s):
%       allData = combine_data_from_same_slice;
%       allData = combine_data_from_same_slice('SliceBase', 'slice4');
%       allData = combine_data_from_same_slice('ChannelTypes', {'voltage', 'current'});
%
% Outputs:
%       allData     - combined data for all slices (or for a particular slice)
%                   specified as a table (or struct)
%
% Arguments:
%       varargin    - 'Directory': working directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'InFolder': directory to read files from
%                   must be a string scalar or a character vector
%                   default == same as directory
%                   - 'SliceBases': name of slice(s) to combine
%                   must be a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == detected from inFolder
%                   - 'SaveMatFlag': whether to save combined data
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'VarsToSave': variables to save
%                   must be a cell array of strings
%                   default == {'sliceBase', 'vVecsSl', 'siMsSl', 'iVecsSl', ...
%                               'phaseBoundaries', 'phaseStrs'}
%                   - 'ForceTableOutput': whether to force output as a table
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for combine_abf_data()
%                           or parse_all_abfs()
%
% Requires:
%       cd/all_files.m
%       cd/all_slice_bases.m
%       cd/combine_abf_data.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/istext.m
%       cd/parse_all_abfs.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/combine_multiunit_data.m

% File History:
% 2019-07-24 Moved from parse_all_multiunit.m
% 2019-07-24 Made 'Directory', 'InFolder', 'SaveMatFlag', 'VarsToSave'
%               optional arguments
% 2019-07-24 Added 'SliceName' as an optional argument
% 2019-08-21 Now uses count_A_each_C.m
% 2019-08-21 Now uses compute_index_boundaries.m
% 2019-08-23 Fixed bug when there is only one phase
% 2019-08-23 Pull out code to function combine_abf_data.m
% TODO: Reorganize code to use structure arrays
% TODO: Allow combination of .mat files
% TODO: Combine gVecs as well

%% Hard-coded parameters
dataExt = 'abf';            % Currently only accepts abf files
sortBy = 'date';
regexpSliceBase = '.*slice[0-9]*';
regexpPhaseStr = 'phase[a-zA-Z0-9]*';

%% Default values for optional arguments
directoryDefault = pwd;
inFolderDefault = '';                   % set later
sliceBasesDefault = {};                 % set later
saveMatFlagDefault = true;
varsToSaveDefault = {'sliceBase', 'vVecsSl', 'siMsSl', 'iVecsSl', ...
                    'phaseBoundaries', 'phaseStrs'};
forceTableOutputDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'InFolder', inFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SliceBases', sliceBasesDefault, ...
    @istext);
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'VarsToSave', varsToSaveDefault, ...
    @iscellstr);
addParameter(iP, 'ForceTableOutput', forceTableOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
inFolder = iP.Results.InFolder;
sliceBases = iP.Results.SliceBases;
saveMatFlag = iP.Results.SaveMatFlag;
varsToSave = iP.Results.VarsToSave;
forceTableOutput = iP.Results.ForceTableOutput;

% Keep unmatched arguments for the combine_abf_data() or parse_all_abfs() 
%       or parse_abf() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the input directory
if isempty(inFolder)
    inFolder = directory;
end

% Decide on the unique slice bases
if isempty(sliceBases)
    % Look for unique slice bases that have .abf files
    sliceBases = all_slice_bases('Directory', inFolder, 'Extension', dataExt, ...
                                'ForceCellOutput', true, 'SortBy', sortBy, ...
                                'RegExpBase', regexpSliceBase);
end

% Make sure sliceBase is a cell array
% Note: The variable 'sliceBase' will be part of the table header
sliceBase = force_column_cell(sliceBases);

% Count the number of slices
nSlices = numel(sliceBase);

%% Do the job
if nSlices > 1 || forceTableOutput
    vVecsSl = cell(nSlices, 1);
    siMsSl = nan(nSlices, 1);
    iVecsSl = cell(nSlices, 1);
    phaseBoundaries = cell(nSlices, 1);
    phaseStrs = cell(nSlices, 1);
%    parfor iSlice = 1:nSlices
    for iSlice = 1:nSlices
        [vVecsSl{iSlice}, siMsSl(iSlice), iVecsSl{iSlice}, ...
            phaseBoundaries{iSlice}, phaseStrs{iSlice}] = ...
            combine_data_from_one_slice(inFolder, sliceBase{iSlice}, ...
                                        saveMatFlag, varsToSave, ...
                                        regexpPhaseStr, otherArguments);
    end

    % Return as a table
    allData = table(sliceBase, vVecsSl, siMsSl, ...
                        iVecsSl, phaseBoundaries, phaseStrs);
elseif nSlices == 1
    [vVecsSl, siMsSl, iVecsSl, phaseBoundaries, phaseStrs] = ...
        combine_data_from_one_slice(inFolder, sliceBase{1}, ...
                                    saveMatFlag, varsToSave, ...
                                    regexpPhaseStr, otherArguments);

    % Return as a struct
    allData.sliceBase = sliceBase{1};
    allData.vVecsSl = vVecsSl;
    allData.siMsSl = siMsSl;
    allData.iVecsSl = iVecsSl;
    allData.phaseBoundaries = phaseBoundaries;
    allData.phaseStrs = phaseStrs;
else
    % Return empty
    allData = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vVecsSl, siMsSl, iVecsSl, phaseBoundaries, phaseStrs] = ...
            combine_data_from_one_slice(inFolder, sliceBase, ...
                                        saveMatFlag, varsToSave, ...
                                        regexpPhaseStr, otherArguments)
%% Combines the data for one slice

%% Count files and phases
% Get all .abf files for this slice ordered by time stamp
[~, allAbfPaths] = all_files('Directory', inFolder, 'Prefix', sliceBase, ...
                            'Extension', 'abf', 'SortBy', 'date');

% If no file found, return
if isempty(allAbfPaths)
    vVecsSl = [];
    siMsSl = NaN;
    iVecsSl = [];
    phaseBoundaries = [];
    phaseStrs = {};
    return
end

%% Combine data
allData = combine_abf_data(allAbfPaths, 'RegExpPhaseStr', regexpPhaseStr, ...
                            'SaveMatFlag', false, 'SaveSheetFlag', true, ...
                            otherArguments{:});

vVecsSl = allData.vVecs;
siMsSl = allData.siMs;
iVecsSl = allData.iVecs;
phaseBoundaries = allData.phaseBoundaries;
phaseStrs = allData.phaseStrs;

%% Save to a matfile if requested
if saveMatFlag
    % Create a matfile path
    commonPrefix = extract_fileparts(allAbfPaths, 'commonprefix');
    commonDir = extract_fileparts(allAbfPaths, 'commondirectory');
    matPath = fullfile(commonDir, [commonPrefix, '.mat']);

    % Save data for this slice
    save(matPath, varsToSave{:}, '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%