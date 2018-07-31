function [ieisGrouped, ieisTable, ieisCellArray, ieisHeader, sheetPath, matfilePath] = ZG_extract_all_IEIs (varargin)
%% Extract all the inter-event intervals from an all_output directory
% Usage: [ieisGrouped, ieisTable, ieisCellArray, ieisHeader, sheetPath, matfilePath] = ZG_extract_all_IEIs (varargin)
% Outputs:
%       TODO
% Arguments:
%       varargin    - 'AllOutputDir': the all_output directory containing
%                                       many subdirectories named by slice_cell
%                   must be a valid directory
%                   default == pwd
%                   - 'OutFileBase': the output file base name
%                   must be a string scalar or a character vector
%                   default == 'ieisGrouped'
%                   - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
%                   - 'ClassesToInclude': classes of events to include
%                   must be a positive integer vector
%                   default == 1:5 (Paula's preference)
%                   - 'TimeWindow': time window to look for events in seconds
%                   must be a numeric vector
%                   default == [min(timeVector), max(timeVector)]
%
% Requires:
%       /home/barrettlab/holySheet/CaImagingExperiments_MasterList.xlsx
%       /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%       /home/Matlab/Adams_Functions/mat2sheet.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/minEASE/filter_minEASE_output.m
%
% Used by:
%       /home/Matlab/Adams_Functions/ZG_compute_IEI_thresholds.m
%
% File History:
% 2018-07-30 Created by Adam Lu
% 

%% Hard-coded parameters
masterListFile = '/home/barrettlab/holySheet/CaImagingExperiments_MasterList.xlsx';
masterListTab = '4mark';
dirPattern = 'data[0-9]*_cell[0-9]*';
dataDirNamesColNum = 3;     % dataDirNames column number in the master list tab
groupLabelsColNum = 5;      % groupLabels column number in the master list tab

%% Default values for optional arguments
allOutputDirDefault = pwd;
outFileBaseDefault = 'ieisGrouped';
sheetTypeDefault = 'xlsx';      % default spreadsheet type
classesToIncludeDefault = 1:5;
timeWindowDefault = [90, 600];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AllOutputDir', allOutputDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'OutFileBase', outFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'ClassesToInclude', classesToIncludeDefault, ...
    @(x) validateattributes(x, {'numeric', 'positive', 'integer'}, {'vector'}));
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, varargin{:});
allOutputDir = iP.Results.AllOutputDir;
outFileBase = iP.Results.OutFileBase;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
classesToInclude = iP.Results.ClassesToInclude;
timeWindow = iP.Results.TimeWindow;

% Check relationships between arguments
if ~isfolder(allOutputDir)
    fprint('%s does not exist or is not readable!\n', allOutputDir);
    return;
end

%% Extract inter-event intervals from each cell
% Get all files under the all output directory
files = dir(fullfile(allOutputDir));

% Filter for subdirectories
subDirs = files(cellfun(@(x) x, {files.isdir}) & ...
                cellfun(@(x) any(regexp(x, dirPattern)), {files.name}));

% Count the number of subdirectories
nSubdirs = numel(subDirs);

% Extract all inter-event intervals and slice and cell labels
ieisAllCells = cell(nSubdirs, 1);
sliceLabelAllCells = cell(nSubdirs, 1);
cellLabelAllCells = cell(nSubdirs, 1);
for iDir = 1:nSubdirs
    % Get the current parent directory
    parentDir = subDirs(iDir).folder;

    % Get the current subdirectory name (the slice_cell label)
    sliceCellLabel = subDirs(iDir).name;

    % Get the full path to the current subdirectory
    fullPath = fullfile(parentDir, sliceCellLabel);

    % Filter the output 
    output = filter_minEASE_output('OutputDir', fullPath, ...
                                    'ClassesToInclude', classesToInclude, ...
                                    'TimeWindow', timeWindow);

    % Extract the inter-event intervals
    ieisAllCells{iDir} = output.interEventIntervals;

    % Extract the slice and cell labels
    tempCell1 = strsplit(sliceCellLabel, '_');
    sliceLabelAllCells{iDir} = tempCell1{1};
    cellLabelAllCells{iDir} = tempCell1{2};
end

%% Combine inter-event intervals of the same slice and of the same group
% Read slice groupings
[~, ~, masterList] = xlsread(masterListFile, masterListTab);
dataDirNames = masterList(2:end, dataDirNamesColNum);
groupLabels = masterList(2:end, groupLabelsColNum);

% Convert the contents to strings if not already so
dataDirNames = cellfun(@num2str, dataDirNames, 'UniformOutput', false);
groupLabels = cellfun(@num2str, groupLabels, 'UniformOutput', false);

% Add an 'x' to the beginning of group labels so that they can be 
%   valid field names of structures
groupLabels = cellfun(@(x) ['x', x], groupLabels, 'UniformOutput', false);

% Find the group labels for all cells
groupLabelAllCells = cell(nSubdirs, 1);
parfor iDir = 1:nSubdirs
    % Get the slice label (dataXXX)
    sliceLabel = sliceLabelAllCells{iDir};

    % Find the matching slice label in the dataDirNames cell array
    %   Note: Must search for substrings and ignore case because 
    %           dataDirNames has  elements of the form DataXXX_caltracer_ORAMA
    idxSlice = find_ind_str_in_cell(sliceLabel, dataDirNames, ...
                                    'SearchMode', 'substrings', ...
                                    'IgnoreCase', true);

    if isempty(idxSlice)
        fprintf('Cannot find the slice label %s in dataDirNames!', sliceLabel);
        % Get the corresponding group label
        groupLabelAllCells{iDir} = 'Not found';
    else
        % Get the corresponding group label
        groupLabelAllCells{iDir} = groupLabels{idxSlice};
    end

end

% Find all unique groups
uniqueGroups = unique(groupLabelAllCells);

% Count the number of unique groups
nGroups = numel(uniqueGroups);

% Loop through all unique groups
for iGroup = 1:nGroups
    % Get the current group label
    groupLabel = uniqueGroups{iGroup};

    % Find all unique slices in this group
    indSlices = find_ind_str_in_cell(groupLabel, groupLabelAllCells);
    sliceLabels = sliceLabelAllCells(indSlices);
    uniqueSlices = unique(sliceLabels);

    % Count the number of unique slices in this group
    nSlices = numel(uniqueSlices);

    % Loop through all slices in this group
    ieisThisGroup = [];    
    for iSlice = 1:nSlices
        % Get the current slice label
        sliceLabel = uniqueSlices{iSlice};

        % Find all unique cells in this slice
        indCells = find_ind_str_in_cell(sliceLabel, sliceLabelAllCells);
        cellLabels = cellLabelAllCells(indCells);
        uniqueCells = unique(cellLabels);

        % Count the number of unique cells in this group
        nCells = numel(uniqueCells);

        % Loop through all cells in this slice
        ieisThisSlice = [];
        for iCell = 1:nCells
            % Get the current cell label
            cellLabel = uniqueCells{iCell};

            % Find the index of this cell in the cell arrays
            idxCell = find_ind_str_in_cell(cellLabel, cellLabelAllCells);

            % Get the ieis for this cell
            ieisThisCell = ieisAllCells{idxCell};

            % Save in structure
            ieisGrouped.(groupLabel).(sliceLabel).(cellLabel).interEventIntervals = ...
                ieisThisCell;

            % Append to ieisThisSlice
            ieisThisSlice = [ieisThisSlice; ieisThisCell];
        end

        % Save in structure
        ieisGrouped.(groupLabel).(sliceLabel).interEventIntervals = ...
            ieisThisSlice;        

        % Combine inter-event intervals from all slices for this group
        ieisThisGroup = [ieisThisGroup; ieisThisSlice];
    end

    % Save in structure
    ieisGrouped.(groupLabel).interEventIntervals = ...
        ieisThisGroup;
end

% Create matfile path
matfilePath = fullfile(allOutputDir, [outFileBase, '.mat']);

% Save in matfile
save(matfilePath, '-struct', 'ieisGrouped');

% Convert to spreadsheet file
[sheetPath, ieisTable, ieisCellArray, ieisHeader] = ...
    mat2sheet(matfilePath, 'SheetType', sheetType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

classesToInclude = 1:5;
timeWindow = [90, 600];

%}