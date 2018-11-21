function SWDs = parse_swd_manual (swdManualFile, varargin)
%% Parse spike-wave-discharge (SWD) event info from manual file
% Usage: SWDs = parse_swd_manual (swdManualFile, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       SWDs        - a table of SWDs event info
%                   specified as a 2D table
% Arguments:
%       swdManualFile   - Manual SWD event file, could be .atf or .csv
%                       must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/atf2sheet.m
%       cd/construct_and_check_abfpath.m
%
% Used by:
%       cd/plot_EEG.m

% File History:
% 2018-11-21 Created by Adam Lu
% 

%% Hard-coded constants
MS_PER_S = 1000;

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

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
addRequired(iP, 'swdManualFile', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, swdManualFile, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Preparation
% Decide what file type the first input is
if regexpi(swdManualFile, '.atf$')
    swdManualAtfFile = swdManualFile;
    swdManualCsvFile = '';
elseif regexpi(swdManualFile, '.csv$')
    swdManualAtfFile = '';
    swdManualCsvFile = swdManualFile;
else
    swdManualAtfFile = '';
    swdManualCsvFile = '';
end

%% Do the job
% Read the table from the file
if ~isfile(swdManualAtfFile) && ~isfile(swdManualCsvFile)
    % Do nothing
    SWDs = [];
    return
elseif isfile(swdManualCsvFile)
    % Display warning if atf file also provided
    if isfile(swdManualAtfFile)
        fprintf(['Table with be read from the csv file ', ...
                'instead of the atf file!\n']);
    end

    % Read in the SWD manual table from the csv file
    swdManualTable = readtable(swdManualCsvFile);
elseif isfile(swdManualAtfFile)
    % Read in the SWD manual table and print to a csv file
    [swdManualTable, swdManualCsvFile] = ...
        atf2sheet(swdManualAtfFile, 'SheetType', 'csv');
end

% Make sure there is an event recorded
if height(swdManualTable) == 0
    % Do nothing
    SWDs = [];
    return
end

% Get the first channel name
firstSignalName = swdManualTable{1, 'Signal'};

% Check whether each row is the same as firstSignalName
isFirstSignal = strcmp(swdManualTable.Signal, firstSignalName);

% Restrict to the entries for the first channel only
swdManualTableOfInterest = swdManualTable(isFirstSignal, :);

% Get the start and end times in ms
startTimesMs = swdManualTableOfInterest.Time1_ms_;
endTimesMs = swdManualTableOfInterest.Time2_ms_;

% Convert to seconds
startTimes = startTimesMs / MS_PER_S;
endTimes = endTimesMs / MS_PER_S;

% Compute durations
durations = endTimes - startTimes;

% If not provided, read in the .abf file names
if isempty(abfFileName)
    % Get the .abf file name for each SWD
    abfFileName = swdManualTable.FileName;
end

% Construct full path to abf file
[abfPath, pathExists] = construct_and_check_abfpath(abfFileName);

% Make sure the dimensions match up
[abfPath, pathExists] = ...
    argfun(@(x), match_dimensions(x, size(startTimes)), abfPath, pathExists);

%% Output results
% Create a table for the parsed SWDs
SWDs = table(startTimes, endTimes, durations, abfPath, pathExists);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%