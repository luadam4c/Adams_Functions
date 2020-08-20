function [ieisGrouped, ieisTable, ieisListed, ieisHeader, sheetPath, groupedFilePath, listedFilePath] = ZG_extract_all_IEIs (varargin)
%% Extract all the inter-event intervals from a directory containing multiple minEASE output subdirectories
% Usage: [ieisGrouped, ieisTable, ieisListed, ieisHeader, sheetPath, groupedFilePath, listedFilePath] = ZG_extract_all_IEIs (varargin)
% Outputs:
%       TODO
% Arguments:
%       varargin    - 'AllOutputDir': the all_output directory containing
%                                       many subdirectories named by slice_cell
%                   must be a valid directory
%                   default == pwd
%                   - 'GroupedFileBase': the base name of the output matfile 
%                                       containing inter-event interval vectors 
%                                       for each data source grouped as a structure
%                   must be a string scalar or a character vector
%                   default == ['ieisGrouped', '_' eventClassStr, ...
%                                   '_', timeWindowStr]
%                   - 'ListedFileBase': the base name of the matfile
%                                       containing inter-event interval vectors 
%                                       for each cell listed as a cell array
%                   must be a string scalar or a character vector
%                   default == ['ieisListed', '_' eventClassStr, ...
%                                   '_', timeWindowStr]
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
%       /home/Matlab/Adams_Functions/find_in_strings.m
%       /home/Matlab/Adams_Functions/mat2sheet.m
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/minEASE/minEASE_filter_output.m
%
% Used by:
%       /home/Matlab/Adams_Functions/ZG_compute_IEI_thresholds.m
%
% File History:
% 2018-07-30 Created by Adam Lu
% 2018-08-01 Updated the directory pattern
% 2018-08-01 Appended classesToInclude and timeWindow to the output file names
% 2018-08-02 Added ieisThisDataset
% 

%% Hard-coded parameters
dirRegexp = 'x[A-Za-z0-9]*_data[0-9]*_cell[0-9]*';
                                % regular expression pattern 
                                %   for each output subdirectory

%% Default values for optional arguments
allOutputDirDefault = pwd;
groupedFileBaseDefault = '';
listedFileBaseDefault = '';
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
addParameter(iP, 'GroupedFileBase', groupedFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ListedFileBase', listedFileBaseDefault, ...
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
groupedFileBase = iP.Results.GroupedFileBase;
listedFileBase = iP.Results.ListedFileBase;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
classesToInclude = iP.Results.ClassesToInclude;
timeWindow = iP.Results.TimeWindow;

% Check relationships between arguments
if ~isfolder(allOutputDir)
    fprint('%s does not exist or is not readable!\n', allOutputDir);
    return;
end

% Set defaults for dependent arguments
eventClassStr = ['eventClass', sscanf(num2str(classesToInclude), '%s')];
timeWindowStr = ['timeWindow', strjoin(strsplit(num2str(timeWindow)), 'to')];
if isempty(groupedFileBase)
    groupedFileBase = ['ieisGrouped', ...
                            '_', eventClassStr, '_', timeWindowStr];
end
if isempty(listedFileBase)
    listedFileBase = ['ieisListed', ...
                            '_', eventClassStr, '_', timeWindowStr];
end

%% Preparation
% Create matfile paths
groupedFilePath = fullfile(allOutputDir, [groupedFileBase, '.mat']);
listedFilePath = fullfile(allOutputDir, [listedFileBase, '.mat']);

%% Extract inter-event intervals from each cell
% Get all files under the all output directory
files = dir(fullfile(allOutputDir));

% Filter for subdirectories
subDirs = files(cellfun(@(x) x, {files.isdir}) & ...
                cellfun(@(x) any(regexp(x, dirRegexp)), {files.name}));

% Count the number of subdirectories
nSubdirs = numel(subDirs);

% Extract all inter-event intervals and slice and cell labels
ieisAllCells = cell(nSubdirs, 1);
groupLabelAllCells = cell(nSubdirs, 1);
sliceLabelAllCells = cell(nSubdirs, 1);
cellLabelAllCells = cell(nSubdirs, 1);
for iDir = 1:nSubdirs
    % Get the current parent directory
    parentDir = subDirs(iDir).folder;

    % Get the current subdirectory name (the slice_cell label)
    groupSliceCellLabel = subDirs(iDir).name;

    % Get the full path to the current subdirectory
    fullPath = fullfile(parentDir, groupSliceCellLabel);

    % Filter the output 
    output = minEASE_filter_output('OutputDir', fullPath, ...
                                    'ClassesToInclude', classesToInclude, ...
                                    'TimeWindow', timeWindow);

    % Extract the inter-event intervals
    ieisAllCells{iDir} = output.interEventIntervals;

    % Extract the group, slice and cell labels
    tempCell1 = strsplit(groupSliceCellLabel, '_');
    groupLabelAllCells{iDir} = tempCell1{1};
    sliceLabelAllCells{iDir} = tempCell1{2};
    cellLabelAllCells{iDir} = tempCell1{3};
end

% Save cell array and labels in a matfile
save(listedFilePath, 'groupLabelAllCells', 'cellLabelAllCells', ...
                    'sliceLabelAllCells', 'ieisAllCells', '-v7.3');

%% Combine inter-event intervals of the same slice and of the same group
% Find all unique groups
uniqueGroups = unique(groupLabelAllCells);

% Count the number of unique groups
nGroups = numel(uniqueGroups);

% Loop through all unique groups
ieisThisDataset = [];
for iGroup = 1:nGroups
    % Get the current group label
    groupLabel = uniqueGroups{iGroup};

    % Find all unique slices in this group
    indSlices = find_in_strings(groupLabel, groupLabelAllCells);
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
        indCells = find_in_strings(sliceLabel, sliceLabelAllCells);
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
            idxCell = find_in_strings(cellLabel, cellLabelAllCells);

            % Get the inter-event intervals for this cell
            ieisThisCell = ieisAllCells{idxCell};

            % Save inter-event intervals from this cell in the structure
            ieisGrouped.(groupLabel).(sliceLabel).(cellLabel).interEventIntervals = ...
                ieisThisCell;

            % Append inter-event intervals from this cell to this slice
            ieisThisSlice = [ieisThisSlice; ieisThisCell];
        end

        % Save inter-event intervals from this slice in the structure
        ieisGrouped.(groupLabel).(sliceLabel).interEventIntervals = ...
            ieisThisSlice;        

        % Append inter-event intervals from this slice to this group
        ieisThisGroup = [ieisThisGroup; ieisThisSlice];
    end

    % Save inter-event intervals from this group in the structure
    ieisGrouped.(groupLabel).interEventIntervals = ...
        ieisThisGroup;

    % Append inter-event intervals from this group to this dataset
    ieisThisDataset = [ieisThisDataset; ieisThisGroup];
end

% Save inter-event intervals from this dataset in the structure
ieisGrouped.interEventIntervals = ieisThisDataset;

% Save structure in a matfile
save(groupedFilePath, '-struct', 'ieisGrouped');

%% Convert structure to spreadsheet file
[sheetPath, ieisTable, ieisListed, ieisHeader] = ...
    mat2sheet(groupedFilePath, 'SheetType', sheetType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

classesToInclude = 1:5;
timeWindow = [90, 600];

% Find the group labels for all cells
groupLabelAllCells = cell(nSubdirs, 1);
parfor iDir = 1:nSubdirs
    % Get the slice label (dataXXX)
    sliceLabel = sliceLabelAllCells{iDir};

end

dirRegexp = 'data[0-9]*_cell[0-9]*';

% Get the current subdirectory name (the slice_cell label)
sliceCellLabel = subDirs(iDir).name;

% Get the full path to the current subdirectory
fullPath = fullfile(parentDir, sliceCellLabel);

% Extract the slice and cell labels
tempCell1 = strsplit(sliceCellLabel, '_');
sliceLabelAllCells{iDir} = tempCell1{1};
cellLabelAllCells{iDir} = tempCell1{2};

%}
