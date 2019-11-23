function thresholdsTable = ZG_extract_IEI_thresholds(fitsGrouped, varargin)
%% Extract/compute inter-event-interval distribution thresholds, separating events from spikes
% Usage: thresholdsTable = ZG_extract_IEI_thresholds(fitsGrouped, varargin)
% Outputs:
%       TODO
% Arguments:
%       fitsGrouped - fits grouped by experimental group, slice and cell
%                   must be a nonempty structure
%       varargin    - 'OutFolder': the output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SheetFileBase': the output spreadsheet file base name
%                   must be a string scalar or a character vector
%                   default == 'thresholdsTable'
%                   - 'MatFileBase': the output matfile base name
%                   must be a string scalar or a character vector
%                   default == 'thresholdsCellArray'
%                   - 'SheetType': sheet type; 
%                       e.g., 'xlsx', 'csv', etc.
%                   could be anything recognised by the readtable() function 
%                   (see issheettype.m under Adams_Functions)
%                   default == 'xlsx'
%
% Requires:
%       /home/Matlab/Adams_Functions/issheettype.m
%
% Used by:
%       /home/Matlab/Adams_Functions/ZG_compute_IEI_thresholds.m
%
% File History:
% 2018-07-30 Adapted code from zgThreshold_MEANS.m
% 2018-08-02 Added 'dataset' as a data source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hard-coded parameters
thresholdsHeader = {'Data Source', ...
            'IEI Kernel Threshold', ...
            'IEI Truncated GaussSingle Threshold', ...
            'IEI GaussOnly Threshold (intersection)', ...
            'IEI GaussOnly Threshold (minimum)', ...
            'IEI GaussOnly Threshold (lower)', ...
            'IEI GaussOnly Threshold (higher)', ...
            'IEI GaussExp Threshold (intersection)', ...
            'IEI GaussExp Threshold (minimum)', ...
            'IEI GaussExp Threshold (lower)', ...
            'IEI GaussExp Threshold (higher)', ...
            'IEI Bodova Threshold', ...
            'log(IEI) Kernel Threshold', ...
            'log(IEI) Truncated GaussSingle Threshold', ...
            'log(IEI) GaussOnly Threshold (intersection)', ...
            'log(IEI) GaussOnly Threshold (minimum)', ...
            'log(IEI) GaussOnly Threshold (lower)', ...
            'log(IEI) GaussOnly Threshold (higher)', ...
            'log(IEI) GaussExpExp Threshold (intersection)', ...
            'log(IEI) GaussExpExp Threshold (minimum)', ...
            'log(IEI) GaussExpExp Threshold (lower)', ...
            'log(IEI) GaussExpExp Threshold (higher)', ...
            };
varNames = {'dataSourceLabels',...
            'IEIKernelThreshold', ...
            'IEITruncatedGaussSingleThreshold', ...
            'IEIGaussOnlyThresholdIntersection', ...
            'IEIGaussOnlyThresholdMinimum', ...
            'IEIGaussOnlyThresholdLower', ...
            'IEIGaussOnlyThresholdHigher', ...
            'IEIGaussExpThresholdIntersection', ...
            'IEIGaussExpThresholdMinimum', ...
            'IEIGaussExpThresholdLower', ...
            'IEIGaussExpThresholdHigher', ...
            'IEIBodovaThreshold', ...
            'logIEIKernelThreshold', ...
            'logIEITruncatedGaussSingleThreshold', ...
            'logIEIGaussOnlyThresholdIntersection', ...
            'logIEIGaussOnlyThresholdMinimum', ...
            'logIEIGaussOnlyThresholdLower', ...
            'logIEIGaussOnlyThresholdHigher', ...
            'logIEIGaussExpExpThresholdIntersection', ...
            'logIEIGaussExpExpThresholdMinimum', ...
            'logIEIGaussExpExpThresholdLower', ...
            'logIEIGaussExpExpThresholdHigher', ...
            };
varSources = {'', ...
            'IEIparams.thresholdModel0', ...
            'IEIparams.thresholdModel1', ...
            'IEIparams.threshold1Model2', ...
            'IEIparams.threshold2Model2', ...
            'min(IEIparams.threshold1Model2, IEIparams.threshold2Model2)', ...
            'max(IEIparams.threshold1Model2, IEIparams.threshold2Model2)', ...
            'IEIparams.threshold1Model3', ...
            'IEIparams.threshold2Model3', ...
            'min(IEIparams.threshold1Model3, IEIparams.threshold2Model3)', ...
            'max(IEIparams.threshold1Model3, IEIparams.threshold2Model3)', ...
            'IEIparams.thresholdModel4', ...
            'exp(logIEIparams.thresholdModel0)', ...
            'exp(logIEIparams.thresholdModel1)', ...
            'exp(logIEIparams.threshold1Model2)', ...
            'exp(logIEIparams.threshold2Model2)', ...
            'exp(min(logIEIparams.threshold1Model2, logIEIparams.threshold2Model2))', ...
            'exp(max(logIEIparams.threshold1Model2, logIEIparams.threshold2Model2))', ...
            'exp(logIEIparams.threshold1Model3)', ...
            'exp(logIEIparams.threshold2Model3)', ...
            'exp(min(logIEIparams.threshold1Model3, logIEIparams.threshold2Model3))', ...
            'exp(max(logIEIparams.threshold1Model3, logIEIparams.threshold2Model3))', ...
            };

%% Default values for optional arguments
outFolderDefault = pwd;
sheetFileBaseDefault = 'thresholdsTable';
matFileBaseDefault = 'thresholdsCellArray';
sheetTypeDefault = 'xlsx';      % default spreadsheet type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'fitsGrouped', ...              % a structure
    @(x) validateattributes(x, {'struct'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'SheetFileBase', sheetFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'MatFileBase', matFileBaseDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, fitsGrouped, varargin{:});
outFolder = iP.Results.OutFolder;
sheetFileBase = iP.Results.SheetFileBase;
matFileBase = iP.Results.MatFileBase;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);

%% Preparation
% Create outFolder if it doesn't exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
end

% Construct the full file name for the output matfile
matFullFileName = fullfile(outFolder, [matFileBase, '.mat']);
                            
% Construct the full file name for the output spreadsheet file
sheetFullFileName = fullfile(outFolder, [sheetFileBase, '.', sheetType]);

%% Determine the number of data sources
% Start a counter
ct = 0;                        % Counts the data sources

% Removed pooled parameters from this data set
thisDataset = rmfield(fitsGrouped, {'IEIparams', 'logIEIparams'});
ct = ct + 1;

% Determine the number of experimental groups and add to counter
allGroups = fieldnames(thisDataset);
nGroups = numel(allGroups);
ct = ct + nGroups;

% Loop through each group
for iGroup = 1:nGroups
    % Extract the identifier for this experimental group
    idGroup = allGroups{iGroup};

    % Extract the structure for this group
    thisGroup = rmfield(fitsGrouped.(idGroup), ...
                        {'IEIparams', 'logIEIparams'});

    % Determine the number of slices in this group and add to counter
    allSlices = fieldnames(thisGroup);
    nSlices = numel(allSlices);
    ct = ct + nSlices;

    % Loop through each slice in this group
    for iSlice = 1:nSlices
        % Extract the identifier for this experimental group
        idSlice = allSlices{iSlice};

        % Extract the structure for this group
        thisSlice = rmfield(fitsGrouped.(idGroup).(idSlice), ...
                            {'IEIparams', 'logIEIparams'});

        % Determine the number of slices in this group and add to counter
        allCells = fieldnames(thisSlice);
        nCells = numel(allCells);
        ct = ct + nCells;
    end
end
nSources = ct;

%% Convert the data to arrays
% Preallocate vectors
dataSourceLabels = cell(nSources, 1);
allIEIparams = cell(nSources, 1);
allLogIEIparams = cell(nSources, 1);

% Start a counter for the row number
row = 0;

% Get the thresholds for pooled data from the entire dataset
row = row + 1;
dataSourceLabels{row} = 'dataset';
allIEIparams{row} = fitsGrouped.IEIparams;
allLogIEIparams{row} = fitsGrouped.logIEIparams;

% Loop through each group first
for iGroup = 1:nGroups
    % Extract the identifier for this experimental group
    idGroup = allGroups{iGroup};

    % Increment the row number
    row = row + 1;

    % Generate the data source label
    dataSourceLabels{row} = idGroup;

    % Add fits data to cell arrays
    allIEIparams{row} = fitsGrouped.(idGroup).IEIparams;
    allLogIEIparams{row} = fitsGrouped.(idGroup).logIEIparams;
end

% Loop through each slice
for iGroup = 1:nGroups
    % Extract the identifier for this experimental group
    idGroup = allGroups{iGroup};

    % Extract the structure for this group
    thisGroup = rmfield(fitsGrouped.(idGroup), ...
                        {'IEIparams', 'logIEIparams'});

    % Determine the slices of this group
    allSlices = fieldnames(thisGroup);
    nSlices = numel(allSlices);

    % Loop through each slice in this group
    for iSlice = 1:nSlices
        % Extract the identifier for this slice
        idSlice = allSlices{iSlice};

        % Increment the row number
        row = row + 1;

        % Generate the data source label
        dataSourceLabels{row} = [idGroup, '_', idSlice];

        % Add fits data to cell arrays
        allIEIparams{row} = ...
            fitsGrouped.(idGroup).(idSlice).IEIparams;
        allLogIEIparams{row} = ...
            fitsGrouped.(idGroup).(idSlice).logIEIparams;
    end
end

% Loop through each cell
for iGroup = 1:nGroups
    % Extract the identifier for this experimental group
    idGroup = allGroups{iGroup};

    % Extract the structure for this group
    thisGroup = rmfield(fitsGrouped.(idGroup), ...
                        {'IEIparams', 'logIEIparams'});

    % Determine the slices of this group
    allSlices = fieldnames(thisGroup);
    nSlices = numel(allSlices);

    % Loop through each slice in this group
    for iSlice = 1:nSlices
        % Extract the identifier for this slice
        idSlice = allSlices{iSlice};

        % Extract the structure for this slice
        thisSlice = rmfield(fitsGrouped.(idGroup).(idSlice), ...
                            {'IEIparams', 'logIEIparams'});

        % Determine the cells for this slice
        allCells = fieldnames(thisSlice);
        nCells = numel(allCells);

        % Loop through each cell in this slice
        for iCell = 1:nCells
            % Extract the identifier for this cell
            idCell = allCells{iCell};

            % Increment the row number
            row = row + 1;

            % Generate the data source label
            dataSourceLabels{row} = [idGroup, '_', idSlice, '_', idCell];

            % Add fits data to cell arrays
            allIEIparams{row} = ...
                fitsGrouped.(idGroup).(idSlice).(idCell).IEIparams;
            allLogIEIparams{row} = ...
                fitsGrouped.(idGroup).(idSlice).(idCell).logIEIparams;
        end
    end
end

%% Compute a threshold for each data source
% Count the number of columns
nCols = numel(thresholdsHeader);

% Create a thresholds cell array
thresholdsCellArray = cell(nSources, nCols);

% Place the data source labels in the first column
thresholdsCellArray(:, 1) = dataSourceLabels;

% Loop through all data sources
for iSource = 1:nSources
    % Extract the data for this source
    IEIparams = allIEIparams{iSource};
    logIEIparams = allLogIEIparams{iSource};

    % Loop through all thresholds
    for iCol = 2:nCols
        % Assign the variable
        eval(sprintf('thresholdsCellArray{iSource, iCol} = %s;\n', ...
                        varSources{iCol}));
    end
end

%% Output
% Convert the cell array into a table
thresholdsTable = cell2table(thresholdsCellArray, ...
                    'VariableNames', varNames);

% Write the table to a spreadsheet file
writetable(thresholdsTable, sheetFullFileName);

% Save the other variables in a matfile
save(matFullFileName, 'thresholdsCellArray', 'thresholdsHeader', ...
                    'varSources', 'varNames', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Preallocate variables
for iCol = 2:nCols
    eval(sprintf('%s = cell(nSources, 1);\n', varNames{iCol}));
end

% Loop through all data sources
parfor iSource = 1:nSources
    % Extract the data for this source
    IEIparams = allIEIparams{iSource};
    logIEIparams = allLogIEIparams{iSource};

    % Loop through all thresholds
    for iCol = 2:nCols
        % Assign the variable
        eval(sprintf('%s{iSource} = %s\n', ...
                        varNames{iCol}, varSources{iCol}));
    end
end

% Preallocate variables
for iCol = 2:nCols
    thresholdsCellArray{iCol} = cell(nSources, 1);
end

thresholdsTable = cell2table(thresholdsCellArray, ...
                    'VariableNames', thresholdsHeader);

%}

