function [swdManualTable, swdManualCsvFile] = ...
                parse_atf_swd (originalEventFile, varargin)
%% Parse spike-wave-discharge (SWD) event info from .atf file or converted .csv file
% Usage: [swdManualTable, swdManualCsvFile] = ...
%               parse_atf_swd (originalEventFile, varargin)
% Explanation:
%       TODO
% Example(s):
%       swdManualTable = parse_atf_swd('WAGS04_30_2018_cage3_Manual_SWDs.atf');
% Outputs:
%       swdManualTable      - a table of spike-wave discharge event info
%                           specified as a 2D table
% Arguments:
%       originalEventFile   - original event file, could be .atf or 
%                               converted .csv
%                           must be a string scalar or a character vector
%       varargin    - 'AbfFileName': Name of the corresponding .abf file(s)
%                   must be empty, a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == [fileBase, '.abf']
%                   - 'OutFolder': directory to output swd table file, 
%                                   e.g. 'output'
%                   must be a string scalar or a character vector
%                   default == same as location of originalEventFile
%
% Requires:
%       cd/argfun.m
%       cd/atf2sheet.m
%       cd/construct_and_check_abfpath.m
%       cd/match_dimensions.m
%
% Used by:
%       cd/plot_traces_EEG.m

% File History:
% 2018-11-21 Created by Adam Lu
% 

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters
varNames = {'startTime', 'endTime', 'duration', 'abfPath', 'pathExists'};

%% Default values for optional arguments
abfFileNameDefault = '';        % set later
outFolderDefault = '';          % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'originalEventFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AbfFileName', abfFileNameDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, originalEventFile, varargin{:});
abfFileName = iP.Results.AbfFileName;
outFolder = iP.Results.OutFolder;

%% Preparation
% Decide what file type the first input is
if regexpi(originalEventFile, '.atf$')
    atfFile = originalEventFile;
    atfCsvFile = '';
elseif regexpi(originalEventFile, '.csv$')
    atfFile = '';
    atfCsvFile = originalEventFile;
else
    atfFile = '';
    atfCsvFile = '';
end

% Get the fileDir and fileBase
[fileDir, fileBase, ~] = fileparts(originalEventFile);

% Decide on the output folder
if isempty(outFolder)
    outFolder = fileDir;
end

% Construct manual SWD table csv file
swdManualCsvFile = fullfile(outFolder, [fileBase, '_Manual_SWDs.csv']);

%% Do the job
% Read the table from the file
if ~isfile(atfFile) && ~isfile(atfCsvFile)
    % Do nothing
    swdManualTable = [];
    swdManualCsvFile = '';
    return
elseif isfile(atfCsvFile)
    % Display warning if atf file also provided
    if isfile(atfFile)
        fprintf(['Table with be read from the csv file ', ...
                'instead of the atf file!\n']);
    end

    % Read in the SWD manual table from the converted csv file
    atfTable = readtable(atfCsvFile);
elseif isfile(atfFile)
    % Read in the SWD manual table and print to a csv file
    [atfTable, atfCsvFile] = atf2sheet(atfFile, 'SheetType', 'csv');
end

% Make sure there is an event recorded
if height(atfTable) == 0
    % Do nothing
    atfTable = [];
    return
end

% Get the first channel name
firstSignalName = atfTable{1, 'Signal'};

% Check whether each row is the same as firstSignalName
isFirstSignal = strcmp(atfTable.Signal, firstSignalName);

% Restrict to the entries for the first channel only
swdManualTableOfInterest = atfTable(isFirstSignal, :);

% Get the start and end times in ms
startTimesMs = swdManualTableOfInterest.Time1_ms_;
endTimesMs = swdManualTableOfInterest.Time2_ms_;

% Convert to seconds
startTime = startTimesMs / MS_PER_S;
endTime = endTimesMs / MS_PER_S;

% Compute duration
duration = endTime - startTime;

% If not provided, read in the .abf file names
if isempty(abfFileName)
    % Get the .abf file name for each SWD
    abfFileName = swdManualTableOfInterest.FileName;
end

% Construct full path to abf file
[abfPath, pathExists] = construct_and_check_abfpath(abfFileName);

% Make sure the dimensions match up
[abfPath, pathExists] = ...
    argfun(@(x) match_dimensions(x, size(startTime)), abfPath, pathExists);

%% Output results
% Create a table for the parsed SWDs
swdManualTable = table(startTime, endTime, duration, abfPath, pathExists, ...
                        'VariableNames', varNames);

% Write the table to a file
writetable(swdManualTable, swdManualCsvFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%