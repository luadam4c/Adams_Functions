function thresholdsTable = ZG_compute_IEI_thresholds (varargin)
%% Compute all possible inter-event interval thresholds from the data within an all_output directory
% Usage: thresholdsTable = ZG_compute_IEI_thresholds (varargin)
% Outputs:
%       TODO
% Arguments:
%       varargin    - 'AllOutputDir': the all_output directory containing
%                                       many subdirectories named by slice_cell
%                   must be a valid directory
%                   default == pwd
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
%       /home/Matlab/Adams_Functions/issheettype.m
%       /home/Matlab/Adams_Functions/ZG_extract_all_IEIs.m
%       /home/Matlab/Adams_Functions/ZG_fit_IEI_distributions.m
%       /home/Matlab/Adams_Functions/ZG_extract_IEI_thresholds.m
%
% File History:
% 2018-07-30 Created by Adam Lu

%% Default values for optional arguments
allOutputDirDefault = pwd;
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
addParameter(iP, 'SheetType', sheetTypeDefault, ...
    @(x) all(issheettype(x, 'ValidateMode', true)));
addParameter(iP, 'ClassesToInclude', classesToIncludeDefault, ...
    @(x) validateattributes(x, {'numeric', 'positive', 'integer'}, {'vector'}));
addParameter(iP, 'TimeWindow', timeWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
                                                % introduced after R2016b
% Read from the Input Parser
parse(iP, varargin{:});
allOutputDir = iP.Results.AllOutputDir;
[~, sheetType] = issheettype(iP.Results.SheetType, 'ValidateMode', true);
classesToInclude = iP.Results.ClassesToInclude;
timeWindow = iP.Results.TimeWindow;

% Check relationships between arguments
if ~isfolder(allOutputDir)
    fprint('%s does not exist or is not readable!\n', allOutputDir);
    return;
end

%% Preparation
% Create strings for eventClass and timeWindow
eventClassStr = ['eventClass', sscanf(num2str(classesToInclude), '%s')];
timeWindowStr = ['timeWindow', strjoin(strsplit(num2str(timeWindow)), 'to')];

% Create output file names for ZG_extract_all_IEIs.m
ieisGroupedFileBase = ['ieisGrouped', '_', eventClassStr, '_', timeWindowStr];
ieisListedFileBase = ['ieisListed', '_', eventClassStr, '_', timeWindowStr];

% Create output directories and file names for ZG_fit_IEI_distributions.m
fitsOutFolderParent = allOutputDir;
fitsOutFileBase = ['fitsGrouped', '_', eventClassStr, '_', timeWindowStr];

% Create output directories and file names for ZG_extract_IEI_thresholds.m
thresholdsOutFolder = allOutputDir;
thresholdsSheetFileBase = ['thresholdsTable', ...
                            '_', eventClassStr, '_', timeWindowStr];
thresholdsMatFileBase = ['thresholdsCellArray', ...
                            '_', eventClassStr, '_', timeWindowStr];

%% Do the job
% Extract all the IEIs
ieisGrouped = ...
    ZG_extract_all_IEIs('AllOutputDir', allOutputDir, ...
                        'GroupedFileBase', ieisGroupedFileBase, ...
                        'ListedFileBase', ieisListedFileBase, ...
                        'SheetType', sheetType, ...
                        'ClassesToInclude', classesToInclude, ...
                        'TimeWindow', timeWindow);

% Fit inter-event-interval distributions and log distributions
fitsGrouped = ...
    ZG_fit_IEI_distributions(ieisGrouped, ...
                            'OutFolderParent', fitsOutFolderParent, ...
                            'OutFileBase', fitsOutFileBase);

% Extract/compute inter-event-interval distribution thresholds, 
%   separating events from spikes
thresholdsTable = ...
    ZG_extract_IEI_thresholds(fitsGrouped, ...
                            'OutFolder', thresholdsOutFolder, ...
                            'SheetFileBase', thresholdsSheetFileBase, ...
                            'MatFileBase', thresholdsMatFileBase, ...
                            'SheetType', sheetType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}